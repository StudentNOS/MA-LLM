# Input: Search terms
# Output: Pre-selected papers based on abstracts

from dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='upon-request')

# Function to fetch abstracts and PubMed IDs from the database (prescreened for titles)
def get_abstracts_and_ids():
    """
    Fetches all abstracts and their corresponding PubMed IDs from the database 
    for papers selected in title screening.
    """
    query = "SELECT pmid, abstract FROM meta_analysis WHERE selected_title = TRUE"
    results = execute_query(query, db=ENSURE)
    return [(result[0], result[1]) for result in results if result[1]]

#Alternatively: Fetch abstracts for ALL extracted PubMed IDs
#def get_all_abstracts():
#    query = "SELECT abstract FROM study"
#    abstracts = execute_query(query, ENSURE)
#    return [abstract[0] for abstract in abstracts]

# Function to create a prompt for OpenAI
def generate_prompt(search_term, abstracts):
    """
    Generates a prompt for OpenAI based on a search term and a list of abstracts.
    """
    prompt = f"Search Term: {search_term}\n\nList of abstracts:\n"
    prompt += "\n".join(abstract for _, abstract in abstracts)
    prompt += "\n\nIdentify which of these abstracts are relevant to the search term."
    
    #Formatting
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
        print("Full Response Text:", response_text)
        split_results = [item.strip("'") for item in response_text.split("', '")]
        print("Split Results:", split_results)
        return split_results
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return []

def match_abstracts_to_ids(initials, abstracts_with_ids):
    """
    Matches the first 10 characters of relevant abstracts to their PubMed IDs.
    """
    matched_ids = []
    for initial in initials:
        for pmid, abstract in abstracts_with_ids:
            if abstract.startswith(initial):
                matched_ids.append(pmid)
                break
    return matched_ids

def save_pmids_to_file(pmids, file_path):
    """
    Saves a list of PubMed IDs to a text file.
    """
    try:
        with open(file_path, 'w') as file:
            for pmid in pmids:
                file.write(f"{pmid}\n")
        print(f"Successfully saved {len(pmids)} PMIDs to {file_path}")
    except Exception as e:
        print(f"Failed to save PMIDs to file: {e}")

def mark_abstract_selections(pmids, selected=True):
    """
    Mark papers as selected based on abstract screening.
    """
    placeholders = ', '.join('?' * len(pmids))
    query = f"UPDATE meta_analysis SET selected_abstract = ? WHERE pmid IN ({placeholders})"
    execute_query(query, params=[selected] + pmids)

def main():
    search_term = input("Enter search term: ")
    abstracts_with_ids = get_abstracts_and_ids()
    prompt = generate_prompt(search_term, abstracts_with_ids)
    screened_abstract_initials = screen_abstracts_with_openai(prompt)
    matched_pmids = match_abstracts_to_ids(screened_abstract_initials, abstracts_with_ids)
    mark_abstract_selections(matched_pmids, selected=True)
    print("PMIDs to be saved:", matched_pmids)
    save_pmids_to_file(matched_pmids, "C:/Users/tillj/Desktop/GPT_Selected.txt")
    print("PubMed IDs saved to GPT_Selected.txt")

if __name__ == "__main__":
    main()
