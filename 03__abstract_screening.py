from dbconnect import execute_query, ENSURE
from openai import OpenAI

# Initialize the OpenAI client with a specified API key
client = OpenAI(api_key='upon request')

# Generator function to fetch abstracts and PubMed IDs in batches from the database
def get_abstracts_in_batches(batch_size=20):
    """
    Fetches batches of abstracts and their corresponding PubMed IDs from the database,
    specifically those that were selected in the previous title screening phase.
    This uses pagination to manage memory and API token limitations.
    """
    offset = 0
    while True:
        # SQL query to fetch a specific batch of abstracts based on the offset
        query = f"SELECT pmid, abstract FROM meta_analysis WHERE selected_title = TRUE LIMIT {batch_size} OFFSET {offset}"
        results = execute_query(query, ENSURE)
        if not results:
            break  # Break the loop if no more results are available
        yield [(result[0], result[1]) for result in results if result[1]]
        offset += batch_size  # Increase the offset to fetch the next batch

# Function to generate a prompt for OpenAI based on a search term and a list of abstracts
def generate_prompt(search_term, abstracts):
    """
    Constructs a prompt string for OpenAI API calls, formatted to handle abstracts. 
    It also includes instructions on how the API should format its output.
    """
    prompt = f"Search Term: {search_term}\n\nList of abstracts:\n"
    prompt += "\n".join(abstract for _, abstract in abstracts)
    prompt += "\n\nIdentify which of these abstracts are relevant to the search term."
    prompt += "\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes. "
    prompt += "Each entry should be the first 10 characters of each relevant abstract. "
    prompt += "For example: 'Example AB', 'Example CD', 'Example EF'."
    return prompt

# Function to screen abstracts using the OpenAI API
def screen_abstracts_with_openai(prompt):
    """
    Calls the OpenAI API with a specific prompt and extracts the relevant part of the API's response,
    focusing on abstracts that match the search criteria.
    """
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "system", "content": prompt}],
            max_tokens=2048
        )
        response_text = response.choices[0].message.content.strip()
        split_results = [item.strip("'") for item in response_text.split("', '")]
        return split_results
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return []

# Function to match the first 10 characters of abstracts to their corresponding PubMed IDs
def match_abstracts_to_ids(initials, abstracts_with_ids):
    """
    Matches the extracted initials from the OpenAI API response to their corresponding PubMed IDs,
    facilitating accurate tracking and selection of relevant studies.
    """
    matched_ids = []
    for initial in initials:
        for pmid, abstract in abstracts_with_ids:
            if abstract.startswith(initial):
                matched_ids.append(pmid)
                break
    return matched_ids

# Function to save PubMed IDs to a text file
def save_pmids_to_file(pmids, file_path):
    """
    Writes the list of PubMed IDs to a specified file, providing a permanent record
    of abstracts selected as relevant.
    """
    try:
        with open(file_path, 'w') as file:
            for pmid in pmids:
                file.write(f"{pmid}\n")
        print(f"Successfully saved {len(pmids)} PMIDs to {file_path}")
    except Exception as e:
        print(f"Failed to save PMIDs to file: {e}")

# Function to mark abstracts as selected in the database
def mark_abstract_selections(pmids, selected=True):
    """
    Updates the database to mark the selected abstracts, indicating they have been accepted
    based on their relevance to the search term.
    """
    placeholders = ', '.join('?' * len(pmids))
    query = f"UPDATE meta_analysis SET selected_abstract = ? WHERE pmid IN ({placeholders})"
    execute_query(query, params=[selected] + pmids)

# Main function to orchestrate the abstract screening process
def main():
    """
    Main function to handle user input, batch processing of abstracts,
    and saving of results. This is where the entire screening process is controlled.
    """
    search_term = input("Enter search term: ")
    all_screened_abstracts = []
    all_matched_pmids = []

    for abstracts_batch in get_abstracts_in_batches():
        prompt = generate_prompt(search_term, abstracts_batch)
        screened_abstract_initials = screen_abstracts_with_openai(prompt)
        matched_pmids = match_abstracts_to_ids(screened_abstract_initials, abstracts_batch)
        all_matched_pmids.extend(matched_pmids)
        all_screened_abstracts.extend(screened_abstract_initials)

    mark_abstract_selections(all_matched_pmids, selected=True)
    save_pmids_to_file(all_matched_pmids, "GPT_Selected.txt")
    print("PubMed IDs saved to GPT_Selected.txt")

if __name__ == "__main__":
    main()
