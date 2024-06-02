from dbconnect import execute_query, ENSURE
from openai import OpenAI

client = OpenAI(api_key='upon request')

# Funktion, um Titel in Chargen (Batches) zu holen
def get_titles_in_batches(batch_size=100):
    offset = 0
    while True:
        # SQL-Abfrage, um eine bestimmte Anzahl von Titeln zu holen
        query = f"SELECT title FROM Initial LIMIT {batch_size} OFFSET {offset}"
        titles = execute_query(query, ENSURE)
        if not titles:
            break
        # Liefert eine Liste von Titeln zurück
        yield [title[0] for title in titles]
        offset += batch_size

# Funktion, um ein Prompt für OpenAI zu erstellen
def generate_prompt(search_term, titles):
    prompt = f"Search Term: {search_term}\n\nList of titles:\n"
    prompt += "\n".join(f"{i}. {title}" for i, title in enumerate(titles, 1))
    
    #PROMPT ENGINEERING
    prompt += "\n\nWhich of these titles are relevant to the search term?"
    
    # Formatierung für GPT (NICHT ÄNDERN)
    prompt += "\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes. "
    prompt += "Each entry should be the first 10 characters of each relevant title. "
    prompt += "For example: 'Example AB', 'Example CD', 'Example EF'."
    return prompt

# Funktion, um die Titel mit OpenAI zu prüfen
def screen_titles_with_openai(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "system", "content": prompt}],
            max_tokens=2048
        )
        response_text = response.choices[0].message.content.strip()
        # Extrahieren der relevanten Teile aus der Antwort
        split_results = [item.strip("'") for item in response_text.split("', '")]
        return split_results
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return []

# Funktion, um zu prüfen, ob ein Titel relevant ist
def is_relevant(title, response_initials):
    return any(fuzz.partial_ratio(title, initial) > 80 for initial in response_initials)

# Funktion, um PubMed-IDs für die relevanten Titel zu bekommen
def get_pubmed_ids_for_titles(screened_titles):
    pmids = []
    for title in screened_titles:
        query = "SELECT pmid FROM Initial WHERE title = ?"
        result = execute_query(query, ENSURE, params=(title,))
        if result:
            pmids.extend([res[0] for res in result])
    return pmids

# Funktion, um die PubMed-IDs in eine Datei zu speichern
def save_pmids_to_file(pmids, file_path):
    with open(file_path, 'w') as file:
        for pmid in pmids:
            file.write(f"{pmid}\n")

# Funktion, um Titel von einer Tabelle in eine andere zu verschieben
def move_titles_to_titles_table(pmids):
    placeholders = ', '.join('?' * len(pmids))
    select_query = f"SELECT * FROM Initial WHERE pmid IN ({placeholders})"
    selected_rows = execute_query(select_query, ENSURE, params=pmids)
    
    for row in selected_rows:
        # Spalten der Datenbank und deren Werte
        columns = ["pmid", "title", "abstract", "fulltext_pdf", "short_summary", "search_terms", "inclusion_criteria", "exclusion_criteria", "n_studies_query", "n_studies_included", "patients", "interventions", "controls", "study_designs"]
        data = dict(zip(columns, row[1:]))
        insert("titles", data, ENSURE)

    delete_query = f"DELETE FROM Initial WHERE pmid IN ({placeholders})"
    execute_query(delete_query, ENSURE, params=pmids)

# Hauptfunktion, die den gesamten Ablauf steuert
def main():
    search_term = input("Enter search term: ")
    all_screened_titles = []

    # Iterieren über alle Titel in Chargen
    for titles_batch in get_titles_in_batches():
        prompt = generate_prompt(search_term, titles_batch)
        screened_title_initials = screen_titles_with_openai(prompt)
        screened_titles = [title for title in titles_batch if is_relevant(title, screened_title_initials)]
        all_screened_titles.extend(screened_titles)

    # PubMed-IDs für relevante Titel holen
    pmids = get_pubmed_ids_for_titles(all_screened_titles)
    # Relevante Titel in eine andere Tabelle verschieben
    move_titles_to_titles_table(pmids)
    # PubMed-IDs in eine Datei speichern
    save_pmids_to_file(pmids, "C:/Users/tillj/Desktop/ENSURE/GPT_Selected.txt")
    print("PubMed IDs saved to GPT_Selected.txt")
    print(f"{len(pmids)} papers have been marked as selected based on title screening.")

# Ausführung der Hauptfunktion
if __name__ == "__main__":
    main()
