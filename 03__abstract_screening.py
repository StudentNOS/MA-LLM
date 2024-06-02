from dbconnect import execute_query, ENSURE
from openai import OpenAI

# Initialize the OpenAI client with a specified API key
client = OpenAI(api_key='upon request')

# Funktion, um Abstracts in Chargen (Batches) zu holen
def get_abstracts_in_batches(batch_size=20):
    offset = 0
    while True:
        # SQL-Abfrage, um eine bestimmte Anzahl von Abstracts zu holen
        query = f"SELECT pmid, abstract FROM titles LIMIT {batch_size} OFFSET {offset}"
        results = execute_query(query, ENSURE)
        if not results:
            break
        # Liefert eine Liste von Abstracts zurück
        yield [(result[0], result[1]) for result in results if result[1]]
        offset += batch_size

# Funktion, um ein Prompt für OpenAI zu erstellen
def generate_prompt(search_term, abstracts):
    prompt = f"Search Term: {search_term}\n\nList of abstracts:\n"
    prompt += "\n".join(abstract for _, abstract in abstracts)
    
    # PROMPT ENGINEERING
    prompt += "\nThese abstracts will be screened to write a meta-analysis."
    prompt += "\nWhich of these abstracts would you select based on the search term?"
    prompt += "\nBe broad in your selection."
    
    # Formatierung für GPT (NICHT ÄNDERN)
    prompt += "\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes. "
    prompt += "Each entry should be the first 10 characters of each relevant abstract. "
    prompt += "For example: 'Example AB', 'Example CD', 'Example EF'."
    return prompt

# Funktion, um die Abstracts mit OpenAI zu prüfen
def screen_abstracts_with_openai(prompt):
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

# Funktion, um Abstracts mit ihren IDs zu verknüpfen
def match_abstracts_to_ids(initials, abstracts_with_ids):
    matched_ids = []
    for initial in initials:
        for pmid, abstract in abstracts_with_ids:
            if fuzz.partial_ratio(abstract[:10], initial) > 80:
                matched_ids.append(pmid)
                break
    return matched_ids

# Funktion, um die PubMed-IDs in eine Datei zu speichern
def save_pmids_to_file(pmids, file_path):
    try:
        with open(file_path, 'w') as file:
            for pmid in pmids:
                file.write(f"{pmid}\n")
        print(f"Successfully saved {len(pmids)} PMIDs to {file_path}")
    except Exception as e:
        print(f"Failed to save PMIDs to file: {e}")

# Funktion, um Abstracts von einer Tabelle in eine andere zu verschieben
def move_abstracts_to_abstracts_table(pmids):
    placeholders = ', '.join('?' * len(pmids))
    select_query = f"SELECT * FROM titles WHERE pmid IN ({placeholders})"
    selected_rows = execute_query(select_query, ENSURE, params=pmids)
    
    for row in selected_rows:
        # Spalten der Datenbank und deren Werte
        columns = ["pmid", "title", "abstract", "fulltext_pdf", "short_summary", "search_terms", "inclusion_criteria", "exclusion_criteria", "n_studies_query", "n_studies_included", "patients", "interventions", "controls", "study_designs"]
        data = dict(zip(columns, row[1:]))
        insert("abstracts", data, ENSURE)

    delete_query = f"DELETE FROM titles WHERE pmid IN ({placeholders})"
    execute_query(delete_query, ENSURE, params=pmids)

# Hauptfunktion, die den gesamten Ablauf steuert
def main():
    search_term = input("Enter search term: ")
    all_screened_abstracts = []
    all_matched_pmids = []

    # Iterieren über alle Abstracts in Chargen
    for abstracts_batch in get_abstracts_in_batches():
        prompt = generate_prompt(search_term, abstracts_batch)
        screened_abstract_initials = screen_abstracts_with_openai(prompt)
        matched_pmids = match_abstracts_to_ids(screened_abstract_initials, abstracts_batch)
        all_matched_pmids.extend(matched_pmids)
        all_screened_abstracts.extend(screened_abstract_initials)

    # Relevante Abstracts in eine andere Tabelle verschieben
    move_abstracts_to_abstracts_table(all_matched_pmids)
    # PubMed-IDs in eine Datei speichern
    save_pmids_to_file(all_matched_pmids, "C:/Users/tillj/Desktop/ENSURE/GPT_Selected.txt")
    print("PubMed IDs saved to GPT_Selected.txt")

# Ausführung der Hauptfunktion
if __name__ == "__main__":
    main()
