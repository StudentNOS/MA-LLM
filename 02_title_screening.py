from dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='upon request')

def get_titles_in_batches(batch_size=100):
    offset = 0
    while True:
        query = f"SELECT title FROM Initial LIMIT {batch_size} OFFSET {offset}"
        titles = execute_query(query, ENSURE)
        if not titles:
            break
        yield [title[0] for title in titles]
        offset += batch_size

def generate_prompt(search_term, titles):
    prompt = f"Search Term: {search_term}\n\nList of titles:\n"
    prompt += "\n".join(f"{i}. {title}" for i, title in enumerate(titles, 1))
    prompt += "\n\nWhich of these titles are relevant to the search term?"
    
    
    #Formatting
    prompt += "\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes. "
    prompt += "Each entry should be the first 10 characters of each relevant title. "
    prompt += "For example: 'Example AB', 'Example CD', 'Example EF'."
    return prompt

def screen_titles_with_openai(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "system", "content": prompt}],
            max_tokens=2048
        )
        response_text = response.choices[0].message.content.strip()
        # Extract the relevant parts from the formatted response
        split_results = [item.strip("'") for item in response_text.split("', '")]
        return split_results
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return []

def is_relevant(title, response_initials):
    return any(fuzz.partial_ratio(title, initial) > 80 for initial in response_initials)

def get_pubmed_ids_for_titles(screened_titles):
    pmids = []
    for title in screened_titles:
        query = "SELECT pmid FROM Initial WHERE title = ?"
        result = execute_query(query, ENSURE, params=(title,))
        if result:
            pmids.extend([res[0] for res in result])
    return pmids

def save_pmids_to_file(pmids, file_path):
    with open(file_path, 'w') as file:
        for pmid in pmids:
            file.write(f"{pmid}\n")

def move_titles_to_titles_table(pmids):
    placeholders = ', '.join('?' * len(pmids))
    select_query = f"SELECT * FROM Initial WHERE pmid IN ({placeholders})"
    selected_rows = execute_query(select_query, ENSURE, params=pmids)
    
    for row in selected_rows:
        columns = ["pmid", "title", "abstract", "fulltext_pdf", "short_summary", "search_terms", "inclusion_criteria", "exclusion_criteria", "n_studies_query", "n_studies_included", "patients", "interventions", "controls", "study_designs"]
        data = dict(zip(columns, row[1:]))
        insert("titles", data, ENSURE)

    delete_query = f"DELETE FROM Initial WHERE pmid IN ({placeholders})"
    execute_query(delete_query, ENSURE, params=pmids)

def main():
    search_term = input("Enter search term: ")
    all_screened_titles = []

    for titles_batch in get_titles_in_batches():
        prompt = generate_prompt(search_term, titles_batch)
        screened_title_initials = screen_titles_with_openai(prompt)
        screened_titles = [title for title in titles_batch if is_relevant(title, screened_title_initials)]
        all_screened_titles.extend(screened_titles)

    pmids = get_pubmed_ids_for_titles(all_screened_titles)
    move_titles_to_titles_table(pmids)

    save_pmids_to_file(pmids, "GPT_Selected.txt")
    print("PubMed IDs saved to GPT_Selected.txt")
    print(f"{len(pmids)} papers have been marked as selected based on title screening.")

if __name__ == "__main__":
    main()
