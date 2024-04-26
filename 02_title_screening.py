from dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='upon request')

# Generator function to fetch titles in batches from the database
def get_titles_in_batches(batch_size=100):
    offset = 0  
    while True:
        query = f"SELECT title FROM meta_analysis LIMIT {batch_size} OFFSET {offset}" 
        titles = execute_query(query, ENSURE)  
        if not titles:
            break  
        yield [title[0] for title in titles]
        offset += batch_size  

# Function to create a prompt for OpenAI based on a given search term and list of titles
def generate_prompt(search_term, titles):
    """
    Constructs a detailed prompt for the OpenAI API call.
    """
    prompt = f"Search Term: {search_term}\n\nList of titles:\n"
    prompt += "\n".join(f"{i}. {title}" for i, title in enumerate(titles, 1))

    #Prompt Engineering Hier
    prompt += "\n\nWhich of these titles are relevant to the search term?"
    return prompt

# Function to screen titles using the OpenAI API
def screen_titles_with_openai(prompt):
    """
    Sends the constructed prompt to the OpenAI API and returns the response.
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

# Function to retrieve PubMed IDs for screened titles
def get_pubmed_ids_for_titles(screened_titles):
    """
    Fetches PubMed IDs for the titles that were identified as relevant by the GPT.
    """
    pmids = []
    for title in screened_titles:
        query = "SELECT pmid FROM meta_analysis WHERE title = ?"  
        result = execute_query(query, ENSURE, params=(title,))
        if result:
            pmids.extend([res[0] for res in result])
    return pmids

# Function to save PubMed IDs to a file
def save_pmids_to_file(pmids, file_path):
    """
    Saves the list of PubMed IDs to a text file at the specified path.
    """
    with open(file_path, 'w') as file:
        for pmid in pmids:
            file.write(f"{pmid}\n")  

# Function to mark titles as selected in the database
def mark_title_selections(pmids, selected=True):
    """
    Updates the database to mark the selected titles, indicating they have passed the screening process.
    """
    placeholders = ', '.join('?' * len(pmids))  
    query = f"UPDATE meta_analysis SET selected_title = ? WHERE pmid IN ({placeholders})"  
    execute_query(query, params=[selected] + pmids)  

# Main function for GPT screening
def main():
    search_term = input("Enter search term: ")  # Prompt user to enter a search term
    all_screened_titles = []  # List to hold all titles that pass the screening

    for titles_batch in get_titles_in_batches():  # Process each batch of titles
        prompt = generate_prompt(search_term, titles_batch)  # Generate the prompt for the current batch
        screened_title_responses = screen_titles_with_openai(prompt)  # Screen the titles
        screened_titles = [title for title in titles_batch if title in screened_title_responses]  # Filter relevant titles
        all_screened_titles.extend(screened_titles)  # Add relevant titles to the list

    pmids = get_pubmed_ids_for_titles(all_screened_titles)  # Fetch PubMed IDs for all screened titles
    mark_title_selections(pmids, selected=True)  # Mark the titles as selected in the database

    save_pmids_to_file(pmids, "GPT_Selected.txt")  # Save the PubMed IDs to a file

    # Print only the essential output
    print("PubMed IDs saved to GPT_Selected.txt")
    print(f"{len(pmids)} papers have been marked as selected based on title screening.")

if __name__ == "__main__":
    main()
