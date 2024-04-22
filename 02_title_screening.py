# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 19:03:55 2024 by Till Adam

This script is designed to screen PubMed titles using the OpenAI API to identify relevant papers based on a user-provided search term. 
The results are then saved and marked in a local database.

"""

from dbconnect import execute_query, ENSURE  # Import custom database functions and database path.
from openai import OpenAI  # Import the OpenAI library to interact with GPT models.

client = OpenAI(api_key='upon-request')  # Create an OpenAI client with the provided API key.

# Generator function to fetch titles in batches from the database
def get_titles_in_batches(batch_size=100):
    offset = 0  # Initialize offset for SQL pagination
    while True:
        query = f"SELECT title FROM meta_analysis LIMIT {batch_size} OFFSET {offset}"  # SQL query to fetch a batch of titles
        titles = execute_query(query, ENSURE)  # Execute the query using a custom function
        if not titles:
            break  # Exit the loop if no more titles are found
        yield [title[0] for title in titles]  # Yield the current batch of titles as a list
        offset += batch_size  # Increase the offset for the next batch

# Function to create a prompt for OpenAI based on a given search term and list of titles
def generate_prompt(search_term, titles):
    """
    Constructs a detailed prompt for the OpenAI API call.
    """
    prompt = f"Search Term: {search_term}\n\nList of titles:\n"
    prompt += "\n".join(f"{i}. {title}" for i, title in enumerate(titles, 1))  # Format titles as a numbered list
    prompt += "\n\nWhich of these titles are relevant to the search term?"  # Question to ask GPT
    return prompt

# Function to screen titles using the OpenAI API
def screen_titles_with_openai(prompt):
    """
    Sends the constructed prompt to the OpenAI API and returns the response.
    """
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",  # Specify the GPT model to use
            messages=[{"role": "system", "content": prompt}],
            max_tokens=2048  # Limit the maximum number of tokens in the response
        )
        return response.choices[0].message.content.strip()  # Extract and return the content from the response
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
        query = "SELECT pmid FROM meta_analysis WHERE title = ?"  # SQL query to fetch PubMed ID for a title
        result = execute_query(query, ENSURE, params=(title,))  # Execute the query with the title as a parameter
        if result:
            pmids.extend([res[0] for res in result])  # Collect all found PubMed IDs
    return pmids

# Function to save PubMed IDs to a file
def save_pmids_to_file(pmids, file_path):
    """
    Saves the list of PubMed IDs to a text file at the specified path.
    """
    with open(file_path, 'w') as file:
        for pmid in pmids:
            file.write(f"{pmid}\n")  # Write each PubMed ID to the file

# Function to mark titles as selected in the database
def mark_title_selections(pmids, selected=True):
    """
    Updates the database to mark the selected titles, indicating they have passed the screening process.
    """
    placeholders = ', '.join('?' * len(pmids))  # Create a placeholder for each PMID
    query = f"UPDATE meta_analysis SET selected_title = ? WHERE pmid IN ({placeholders})"  # SQL query to update selection status
    execute_query(query, params=[selected] + pmids)  # Execute the update query

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
