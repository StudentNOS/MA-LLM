# Input: Search terms
# Output: Pre-selected papers based on titles

from dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='upon-request')

# Function to fetch all titles from the database
def get_all_titles():
    """
    Fetches all titles from the database.
    """
    query = "SELECT title FROM meta_analysis"  # Adjusted to the correct table
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

def get_pubmed_ids_for_titles(screened_titles):
    """
    Fetches PubMed IDs for the given list of screened titles.
    """
    pmids = []
    for title in screened_titles:
        # It's crucial that this matches exactly; consider potential variations in title formatting.
        query = "SELECT pmid FROM meta_analysis WHERE title = ?"
        result = execute_query(query, ENSURE, params=(title,))
        if result:
            pmids.extend([res[0] for res in result])
    return pmids

def save_pmids_to_file(pmids, file_path):
    """
    Saves a list of PubMed IDs to a text file.
    """
    with open(file_path, 'w') as file:
        for pmid in pmids:
            file.write(f"{pmid}\n")


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
    save_pmids_to_file(pmids, "C:/Users/tillj/Desktop/GPT_screening_titles.txt")

    print("PubMed IDs saved to GPT_screening_titles.txt")

if __name__ == "__main__":
    main()
