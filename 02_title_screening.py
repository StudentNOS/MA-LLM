# Input: Search terms
# Output: Pre-selected papers based on titles

from utils.dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='API-key upon request')

# Function to fetch all titles from the database
def get_all_titles():
    """
    Fetches all titles from the database.
    """
    query = "SELECT title FROM study"
    titles = execute_query(query, ENSURE)
    return [title[0] for title in titles]

# Function to create a prompt for OpenAI
def generate_prompt(search_term, titles):
    """
    Generates a prompt for OpenAI based on a search term and a list of titles.
    """
    prompt = f"Search Term: {search_term}\n\nList of titles:\n"
    prompt += "\n".join(f"{i}. {title}" for i, title in enumerate(titles, 1))
    prompt += "\n\nWhich of these titles are relevant to the search term?"
    return prompt

# Function to screen titles using OpenAI API
def screen_titles_with_openai(prompt):
    """
    Screens titles using the OpenAI API.
    """
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "system", "content": prompt}],
            max_tokens=2048
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return ""

def get_pubmed_ids_for_titles(titles):
    """
    Fetches PubMed IDs for the given list of titles.
    """
    pmids = []
    for title in titles:
        query = "SELECT pmid FROM study WHERE title = ?"
        result = execute_query(query, ENSURE, params=(title,))
        if result:
            pmids.append(result[0][0])
    return pmids

def save_pmids_to_file(pmids, file_path):
    """
    Saves a list of PubMed IDs to a text file.
    """
    with open(file_path, 'w') as file:
        for pmid in pmids:
            file.write(f"{pmid}\n")

def read_ids_from_file(file_path):
    """
    Reads IDs from a given file and returns them as a set.
    """
    with open(file_path, 'r') as file:
        return set(file.read().splitlines())

#Can potentially be adjusted to calculate per thousand instead of percent (when dealing with many PubMed-IDs)
def calculate_percentage(contained_set, total_set):
    """
    Calculates the percentage of IDs in total_set accounted for by IDs in contained_set.
    """
    if total_set:
        return (len(contained_set.intersection(total_set)) / len(total_set)) * 100
    else:
        return 0

def performance_indicator(gpt_file_path, gold_file_path, initial_search_file_path):
    """
    Compares IDs from three files and calculates the matching percentages.
    """
    gpt_ids = read_ids_from_file(gpt_file_path)
    gold_ids = read_ids_from_file(gold_file_path)
    initial_search_ids = read_ids_from_file(initial_search_file_path)

    # Calculate percentages
    gpt_in_gold_percentage = calculate_percentage(gpt_ids, gold_ids)  # Percentage of GPT_screening in gold_ref
    gold_in_initial_percentage = calculate_percentage(gold_ids, initial_search_ids)  # Percentage of gold_ref in initial_search
    gpt_in_initial_percentage = calculate_percentage(gpt_ids, initial_search_ids)  # Percentage of GPT_screening in initial_search

    return gpt_in_gold_percentage, gold_in_initial_percentage, gpt_in_initial_percentage

# Main function for GPT screening
def main():
    search_term = input("Enter search term: ")
    titles = get_all_titles()
    prompt = generate_prompt(search_term, titles)
    screened_title_responses = screen_titles_with_openai(prompt)

    # Process OpenAI's response to get the list of relevant titles
    screened_titles = [title for title in titles if title in screened_title_responses]

    print("Screened titles:")
    print(screened_titles)

    # Get PubMed IDs for screened titles
    pmids = get_pubmed_ids_for_titles(screened_titles)

    # Save PubMed IDs to a text file
    save_pmids_to_file(pmids, "GPT_screening_titles.txt")

    print("PubMed IDs saved to GPT_screening_titles.txt")
    
    gpt_file_path = "GPT_screening_titles.txt"
    gold_file_path = "gold_ref.txt"
    initial_search_file_path = "initial_search.txt"
    percentages = performance_indicator(gpt_file_path, gold_file_path, initial_search_file_path)
    
    print(f"Precision (GPT/Gold Standard): {percentages[0]:.2f}%")
    print(f"Proportion covered (Gold Standard): {percentages[1]:.2f}%")
    print(f"Proportion covered (GPT): {percentages[2]:.2f}%")

if __name__ == "__main__":
    main()
