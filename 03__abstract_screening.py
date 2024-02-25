# Input: Search terms
# Output: Pre-selected papers based on abstracts

from dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='upon-request')

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
    # Filter out empty abstracts
    abstracts = [abstract for abstract in abstracts if abstract]

    # Create a prompt with the search term and list of abstracts
    prompt = f"Search Term: {search_term}\n\nList of abstracts:\n"
    prompt += "\n".join(abstract for abstract in abstracts)

    # Ask which abstracts are relevant to the search term
    prompt += "\n\nIdentify which of these abstracts are relevant to the search term."

    # Instruct to provide output in the specified format
    prompt += "\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes. "
    prompt += "Each entry should be the first 10 characters of each relevant abstract. "
    prompt += "For example: 'Example AB', 'Example CD', 'Example EF'."

    return prompt
# Function to screen abstracts using OpenAI API
def screen_abstracts_with_openai(prompt):
    """
    Screens abstracts using the OpenAI API and extracts the first 10 characters of relevant abstracts.
    """
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "system", "content": prompt}],
            max_tokens=2048
        )
        response_text = response.choices[0].message.content.strip()

        # Debugging: Print the complete response text
        print("Full Response Text:", response_text)

        # Split the response at each occurrence of ', '
        # Trim the single quotes from each item
        split_results = [item.strip("'") for item in response_text.split("', '")]

        # Debugging: Print the split results
        print("Split Results:", split_results)

        return split_results
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return []

def get_pubmed_ids_for_abstracts(abstract_initials, abstracts):
    pmids = []
    for initial in abstract_initials:
        initial = initial.strip("'")
        matching_abstracts = [abstract for abstract in abstracts if abstract.startswith(initial)]
        for abstract in matching_abstracts:
            query = "SELECT pmid FROM study WHERE abstract = ?"
            result = execute_query(query, ENSURE, params=(abstract,))
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
    abstracts = get_all_abstracts()
    prompt = generate_prompt(search_term, abstracts)
    screened_abstract_initials = screen_abstracts_with_openai(prompt)

    # Print OpenAI Response
    print("OpenAI Response with Abstract Initials:")
    print(screened_abstract_initials)

    # Get PubMed IDs for screened abstracts based on their initials
    pmids = get_pubmed_ids_for_abstracts(screened_abstract_initials, abstracts)

    # Save PubMed IDs to a text file
    save_pmids_to_file(pmids, "C:/Users/tillj/Desktop/GPT_screening_abstracts.txt")
    print("PubMed IDs saved to GPT_screening_abstracts.txt")

if __name__ == "__main__":
    main()
