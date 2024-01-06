# Input: Search terms
# Output: Pre-selected papers based on abstracts

from dbconnect import execute_query, ENSURE
from openai import OpenAI
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

client = OpenAI(api_key='API-Key upon request')

# Function to fetch all abstracts from the database
def get_all_abstracts():
    """
    Fetches all abstracts from the database.
    """
    query = "SELECT abstract FROM study"
    abstracts = execute_query(query, ENSURE)
    return [abstract[0] for abstract in abstracts]

# Function to create a prompt for OpenAI
def generate_prompt(search_term, abstracts):
    """
    Generates a prompt for OpenAI based on a search term and a list of abstracts.
    """
    prompt = f"Search Term: {search_term}\n\nList of abstracts:\n"
    prompt += "\n".join(f"{i}. {abstract}" for i, abstract in enumerate(abstracts, 1))
    prompt += "\n\nWhich of these abstracts are relevant to the search term?"
    return prompt

# Function to screen abstracts using OpenAI API
def screen_abstracts_with_openai(prompt):
    """
    Screens abstracts using the OpenAI API.
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

def get_pubmed_ids_for_abstracts(abstracts):
    """
    Fetches PubMed IDs for the given list of abstracts.
    """
    pmids = []
    for abstract in abstracts:
        query = "SELECT pmid FROM study WHERE abstract = ?"
        result = execute_query(query, ENSURE, params=(abstract,))
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
    abstracts = get_all_abstracts()
    prompt = generate_prompt(search_term, abstracts)
    screened_abstract_responses = screen_abstracts_with_openai(prompt)
    logging.info(f"OpenAI Response: {screened_abstract_responses}")

    # Process OpenAI's response to get the list of relevant abstracts
    screened_abstracts = [abstract for abstract in abstracts if abstract in screened_abstract_responses]

    print("Screened abstracts:")
    print(screened_abstracts)

    # Get PubMed IDs for screened abstracts
    pmids = get_pubmed_ids_for_abstracts(screened_abstracts)

    # Save PubMed IDs to a text file
    save_pmids_to_file(pmids, "GPT_screening_abstracts.txt")

    print("PubMed IDs saved to GPT_screening_abstracts.txt")
    
    gpt_file_path = "GPT_screening_abstracts.txt"
    gold_file_path = "gold_ref.txt"
    initial_search_file_path = "initial_search.txt"
    percentages = performance_indicator(gpt_file_path, gold_file_path, initial_search_file_path)
    
    print(f"GPT/Gold Standard (should be as high as possible): {percentages[0]:.2f}%")
    print(f"Gold Standard/Initial Search (reference value): {percentages[1]:.2f}%")
    print(f"GPT/Initial Search (should be as close to reference value as possible): {percentages[2]:.2f}%")

if __name__ == "__main__":
    main()
