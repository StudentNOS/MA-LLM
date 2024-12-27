# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 16:43:00 2024

@author: tillj
"""


from dbconnect import execute_query, insert, ENSURE
from openai import OpenAI
import json
import re
#from fuzzywuzzy import fuzz
#from PATHS import api_key

#client = OpenAI(api_key="OPENAI_API_KEY_NOS")

client = OpenAI(
    base_url= "https://llm1-compute.cms.hu-berlin.de/v1/",
    api_key="required-but-not-used",  # Placeholder as HU LLM does not validate keys
)

#%% Fetching
def get_data_in_batches(decision, batch_size=None):
    offset = 0
    
    if decision == "titles":
        if batch_size is None:
            batch_size = 10
        table = "Initial"
        fields = "pmid, title"
    elif decision == "abstracts":
        if batch_size is None:
            batch_size = 2
        # Check if the "titles" table has entries
        titles_count = execute_query("SELECT COUNT(*) FROM titles", ENSURE)[0][0]
        if titles_count > 0:
            table = "titles"
        else:
            table = "Initial"
        fields = "pmid, abstract"
    else:
        raise ValueError("Invalid decision parameter")
    
    while True:
        query = f"SELECT {fields} FROM {table} LIMIT {batch_size} OFFSET {offset}"
        try:
            results = execute_query(query, ENSURE)
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

def parse_json_llm_response(response: str):
    """
    Parses the response from the LLM to extract JSON content.

    The function first attempts to parse the entire response as JSON. If that fails,
    it looks for a JSON block enclosed in a ```json markdown block and attempts
    to parse it. If successful, it returns the parsed JSON object; otherwise, it returns None.

    :param response: The response string from the LLM.
    :type response: str
    :return: The parsed JSON object or None if parsing fails.
    :rtype: dict or None
    """
    # First, try to parse the entire response as JSON
    try:
        return json.loads(response)
    except json.JSONDecodeError:
        pass  # Proceed to regex if parsing fails

    # Use regex to find the JSON block within the markdown
    json_match = re.search(r'```json\n((?:.|\n)*)\n\s*```', response, re.DOTALL)

    if json_match:
        json_str = json_match.group(1).strip()  # Extract the JSON string
        return json.loads(json_str)  # Parse the JSON string into a Python dictionary

    return None  # Return None if no JSON block is found

#%% Prompt
def generate_prompt(data, decision, manual):
    if decision == "titles":
        data_type = "titles"
        formatted_data = "\n".join(f"{i}. {title}" for i, title in enumerate(data, 1))
        format_instruction = "title"
    elif decision == "abstracts":
        data_type = "abstracts"
        formatted_data = "\n".join(f"{i}. {abstract}" for i, abstract in enumerate(data, 1))
        format_instruction = "abstract"
    else:
        raise ValueError("Invalid decision parameter")

    prompt = (
        f"\n\nList of {data_type} for screening:\n"
        f"{formatted_data}"
        f"\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes."
        f"Each entry should be the unique PMID of the {format_instruction}."
        f"For example: '36823236', '35841079', '21608000'."
        f"\n\nOutput this and only this. Format your response this way without exception.")
    
    #print(f"Generated Prompt: {prompt}")
    
    return prompt


#%% Screening
def screen_with_openai(prompt):
    try:
        response = client.chat.completions.create(
            model='ModelCloud/Mistral-Large-Instruct-2407-gptq-4bit',
            messages=[{"role": "user", "content": prompt}],
            max_tokens=2048,
            response_format={"type": "json_object"}
        )
        response_text = response.choices[0].message.content #.strip()
        
        # Minimal addition: print the raw response for manual checking
        print(f"Raw LLM Response: {response_text}")

        response_dict = parse_json_llm_response(response_text)
        return response_dict
        
        #split_results = [item.strip(" '") for item in response_text.split(",")]
        
        #print(f"Split: {split_results}")
        
        return split_results
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return []

#%% Matching
def match_data_to_ids(initials, data, decision):
    # Extract PMIDs from the data and match against `initials` (IDs returned by LLM)
    matched_ids = [pmid for pmid, _ in data if pmid in initials]
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


#%% Fulltext prompt
def screen_pdf_with_openai(pdf_text, search_term, inclusion, exclusion, manual_fulltext):
    # inflexible prompt characteristics -> innate part of screening
    prompt = (
        f"\n\nSearch Term: {search_term}"
        f"\n\nInclusion Criteria: {inclusion}"
        f"\nExclusion Criteria: {exclusion}"
        f"\n\nFull-text:\n{pdf_text}"
        f"\n\nFormat the output as 'Relevant: first 10 characters of the text' or 'Irrelevant: -'."
        f"\nFor example: \nRelevant: Example AB\nRelevant: Example CD\nIrrelevant: -\nRelevant: Example EF\nIrrelevant: -"
        f"\nOutput this and only this. Format your response this way without exception.")
    
    try:
        response = client.chat.completions.create(
            messages=[{"role": "system", "content": prompt}],
            model="gpt-3.5-turbo",
            max_tokens=2048,
            response_format={"type": "json_object"}
        )
        response_text = response.choices[0].message.content #.strip()
        print(f"Raw response content: {response_text}")
        response_dict = parse_json_llm_response(response_text)
        return response_dict
    
    except Exception as e:
        print(f"Error in OpenAI API call: {e}")
        return "Irrelevant: -"

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
