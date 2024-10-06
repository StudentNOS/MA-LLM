# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 16:43:00 2024

@author: tillj
"""


from dbconnect import execute_query, insert, ENSURE
from openai import OpenAI
from fuzzywuzzy import fuzz
from PATHS import api_key

client = OpenAI(api_key=api_key)
#%% Fetching
def get_data_in_batches(decision, batch_size=None):
    offset = 0
    
    if decision == "titles":
        if batch_size is None:
            batch_size = 100
        table = "Initial"
        fields = "pmid, title"
    elif decision == "abstracts":
        if batch_size is None:
            batch_size = 20
        table = "titles"
        fields = "pmid, abstract"
    else:
        raise ValueError("Invalid decision parameter")
    
    while True:
        query = f"SELECT {fields} FROM {table} LIMIT {batch_size} OFFSET {offset}"
        try:
            results = execute_query(query, ENSURE)
            #print(f"Query: {query}")  # Debugging line
            #print(f"Results: {results}")  # Debugging line
            if not results:
                break
            
            if decision == "titles":
                yield [(result[0], result[1]) for result in results]  # Yield PMIDs and titles
            elif decision == "abstracts":
                yield [(result[0], result[1]) for result in results if result[1]]  # Yield PMIDs and non-empty abstracts
            
        except Exception as e:
            print(f"Error executing query: {e}")
            break
        
        offset += batch_size

#%% Prompt
def generate_prompt(search_term, inclusion, exclusion, data, decision, manual):
    if decision == "titles":
        data_type = "titles"
        formatted_data = "\n".join(f"{i}. {title}" for i, title in enumerate(data, 1))
        format_instruction = "first 10 characters of each relevant title"
    elif decision == "abstracts":
        data_type = "abstracts"
        formatted_data = "\n".join(abstract for _, abstract in data)
        format_instruction = "first 10 characters of each relevant abstract"
    else:
        raise ValueError("Invalid decision parameter")

    prompt = f"Search Term: {search_term}\n\nList of {data_type}:\n"
    prompt += "\n\n Inclusion Criteria: {inclusion}"
    prompt += "\n Exclusion Criteria: {exclusion}"
    prompt += manual
    prompt += formatted_data
    
    # inflexible prompt characteristics -> innate part of screening
    prompt += "\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes. "
    prompt += f"Each entry should be the {format_instruction}. "
    prompt += "For example: 'Example AB', 'Example CD', 'Example EF'."
    prompt += "\n\nOutput this and only this. Format your response this way without exception."
    
    return prompt


#%% Screening
def screen_with_openai(prompt):
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

#%% Fulltext prompt
def screen_pdf_with_openai(pdf_text, search_term, inclusion, exclusion, manual_fulltext):
    # inflexible prompt characteristics -> innate part of screening
    prompt = f"Search Term: {search_term}\n"
    prompt += "\n\n Inclusion Criteria: {inclusion}"
    prompt += "\n Exclusion Criteria: {exclusion}"
    prompt += manual_fulltext
    prompt += "\n\n Full-text:\n{pdf_text}"
    
    prompt += "\n\nFormat the output as 'Relevant: first 10 characters of the text' or 'Irrelevant: -'."
    prompt += "\nFor example: \nRelevant: Example AB\nRelevant: Example CD\nIrrelevant: -\nRelevant: Example EF\nIrrelevant: -"
    prompt += "\nOutput this and only this. Format your response this way without exception."
    
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "system", "content": prompt}],
            max_tokens=2048
        )
        response_text = response.choices[0].message.content.strip()
        return response_text
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return "Irrelevant: -"

#%% Matching
def match_data_to_ids(initials, data, decision):
    matched_ids = []
    
    if decision == "titles":
        matched_ids = [
            pmid for i, (pmid, title) in enumerate(data)
            if any(fuzz.partial_ratio(title, initial) > 90 for initial in initials)
        ]
    elif decision == "abstracts":
        for initial in initials:
            for pmid, abstract in data:
                if fuzz.partial_ratio(abstract[:10], initial) > 90:
                    matched_ids.append(pmid)
                    break
    return matched_ids

#%% Save IDs
def save_pmids_to_file(pmids, file_path):
    try:
        with open(file_path, 'w') as file:
            for pmid in pmids:
                file.write(f"{pmid}\n")
        #print(f"Successfully saved {len(pmids)} PMIDs to {file_path}")
    except Exception as e:
        print(f"Failed to save PMIDs to file: {e}")

#%% Move selected to next table
def move_records(pmids, decision):
    placeholders = ', '.join('?' * len(pmids))
    
    if decision == "titles":
        source_table = "Initial"
        destination_table = "titles"
    elif decision == "abstracts":
        source_table = "titles"
        destination_table = "abstracts"
    elif decision == "full_texts":
        source_table = "abstracts"
        destination_table = "full_texts"
    else:
        print(f"Unknown decision: {decision}")
        return
    
    # Debugging
    #print(f"PMIDs to move: {pmids}")
    
    # select rows from the source table
    select_query = f"SELECT * FROM {source_table} WHERE pmid IN ({placeholders})"
    selected_rows = execute_query(select_query, ENSURE, params=pmids)
    #print(f"Selected Rows: {selected_rows}")  # Debugging
    
    # move each row to destination table
    for row in selected_rows:
        data = {
            "pmid": row[1],
            "title": row[2],
            "abstract": row[3],
            "fulltext_pdf": row[4],
            "short_summary": row[5],
            "search_terms": row[6],
            "inclusion_criteria": row[7],
            "exclusion_criteria": row[8],
            "n_studies_query": row[9],
            "n_studies_included": row[10],
            "patients": row[11],
            "interventions": row[12],
            "controls": row[13],
            "study_designs": row[14]
        }
        insert(destination_table, data, ENSURE, serialize_json=False)
    
    # delete rows from the source table
    delete_query = f"DELETE FROM {source_table} WHERE pmid IN ({placeholders})"
    execute_query(delete_query, ENSURE, params=pmids)

#%% Fulltext IDs to txt file 
def save_fulltext_pmids_to_file(database_path, output_file):
    try:
        import sqlite3
        # Connect to the SQLite database
        conn = sqlite3.connect(database_path)
        cursor = conn.cursor()

        # Query to select PMIDs from full_texts table
        cursor.execute("SELECT pmid FROM full_texts")
        pmids = cursor.fetchall()

        # Write PMIDs to the output file
        with open(output_file, 'w') as file:
            for pmid in pmids:
                file.write(f"{pmid[0]}\n")

        #print(f"Successfully saved {len(pmids)} PMIDs to {output_file}")
    
    except Exception as e:
        print(f"Failed to save PMIDs to file: {e}")
    
    finally:
        if conn:
            conn.close()
