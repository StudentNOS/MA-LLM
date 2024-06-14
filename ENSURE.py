# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 17:02:52 2024

@author: tillj
"""

from PATHS import ENSURE
from dbconnect import init_db, insert, delete_all_data
from PubMed import count_papers, create_excel_from_db, read_pmids_from_file, fetch_details 
from GPT import get_data_in_batches, generate_prompt, screen_with_openai, match_data_to_ids, save_pmids_to_file, move_records
from Performance import read_ids_from_file, calculate_performance_metrics

db_name = ENSURE
init_db(db_name)
print("Database initialized with the necessary tables.")

def main():
    # Step 1: Read PMIDs and fetch details
    pmid_file = "Initial.txt"
    print(f"Reading PMIDs from {pmid_file}...")
    pmids = read_pmids_from_file(pmid_file)
    print(f"Fetching details for {len(pmids)} PMIDs...")
    records = fetch_details(pmids)
    
    # Step 2: Insert fetched records into the 'Initial' table
    print("Inserting fetched records into the 'Initial' table...")
    for article in records.get("PubmedArticle", []):
        pmid = article["MedlineCitation"]["PMID"]
        title = article["MedlineCitation"]["Article"].get("ArticleTitle", "")
        abstract_text = ""
        if "Abstract" in article["MedlineCitation"]["Article"]:
            abstract_text = " ".join(str(abstract) for abstract in article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"])
        
        insert("Initial", {"pmid": pmid, "title": title, "abstract": abstract_text}, ENSURE, False)
        print(f"Inserted: {pmid}")
    
    paper_count = count_papers()
    print(f"Total number of papers in the database: {paper_count}")
    
    # Step 3: Generate Excel report from the database
    print("Generating Excel report from the database...")
    create_excel_from_db()

    # Step 4: Ask user if they want to screen titles
    screen_titles = input("Do you want to screen the titles (y/n)? ").strip().lower() == 'y'
    if screen_titles:
        search_term = input("Enter search term for titles: ").strip()
        manual = input("Prompt: ").strip()
        all_screened_titles = []
        all_matched_pmids_titles = []

        # Process titles
        print("Processing titles...")
        for titles_batch in get_data_in_batches("titles"):
            prompt = generate_prompt(search_term, titles_batch, "titles", manual)
            screened_title_initials = screen_with_openai(prompt)
            matched_pmids_titles = match_data_to_ids(screened_title_initials, titles_batch, "titles")
            
            all_matched_pmids_titles.extend(matched_pmids_titles)
            all_screened_titles.extend(screened_title_initials)

        # Move relevant titles to the 'titles' table
        print("Moving relevant titles to the 'titles' table...")
        move_records(all_matched_pmids_titles, "titles")

        # Save PubMed IDs for titles to a file
        output_file_titles = "GPT_Selected_Titles.txt"
        print(f"Saving PubMed IDs for titles to {output_file_titles}...")
        save_pmids_to_file(all_matched_pmids_titles, output_file_titles)
        
        print(f"PubMed IDs saved to {output_file_titles}")
        print(f"{len(all_matched_pmids_titles)} titles have been marked as selected based on title screening.")
        
        # Step 5: Calculate performance metrics for titles
        initial_ids = read_ids_from_file(pmid_file)
        gpt_selected_ids = read_ids_from_file(output_file_titles)
        goldstandard_selected_ids = read_ids_from_file("Goldstandard_Selected.txt")  # Assuming the gold standard file is located here
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial_ids, gpt_selected_ids, goldstandard_selected_ids)
        print("Performance metrics for titles:")
        print(f"Sensitivity: {sensitivity:.4f}")
        print(f"Specificity: {specificity:.4f}")
        print(f"Positive Predictive Value (PPV): {PPV:.4f}")
        print(f"Negative Predictive Value (NPV): {NPV:.4f}")
        print(f"True Positives: {tp}")
        print(f"True Negatives: {tn}")
        print(f"False Positives: {fp}")
        print(f"False Negatives: {fn}")

    # Step 6: Ask user if they want to screen abstracts
    screen_abstracts = input("Do you want to screen the abstracts (y/n)? ").strip().lower() == 'y'
    if screen_abstracts:
        search_term = input("Enter search term for abstracts: ").strip()
        manual = input("Prompt: ").strip()
        all_screened_abstracts = []
        all_matched_pmids_abstracts = []

        # Process abstracts
        print("Processing abstracts...")
        for abstracts_batch in get_data_in_batches("abstracts"):
            prompt = generate_prompt(search_term, abstracts_batch, "abstracts", manual)
            screened_abstract_initials = screen_with_openai(prompt)
            matched_pmids_abstracts = match_data_to_ids(screened_abstract_initials, abstracts_batch, "abstracts")
            
            all_matched_pmids_abstracts.extend(matched_pmids_abstracts)
            all_screened_abstracts.extend(screened_abstract_initials)

        # Move relevant abstracts to the 'abstracts' table
        print("Moving relevant abstracts to the 'abstracts' table...")
        move_records(all_matched_pmids_abstracts, "abstracts")

        # Save PubMed IDs for abstracts to a file
        output_file_abstracts = "GPT_Selected_Abstracts.txt"
        print(f"Saving PubMed IDs for abstracts to {output_file_abstracts}...")
        save_pmids_to_file(all_matched_pmids_abstracts, output_file_abstracts)
        
        print(f"PubMed IDs saved to {output_file_abstracts}")
        print(f"{len(all_matched_pmids_abstracts)} abstracts have been marked as selected based on abstract screening.")
        
        # Step 7: Calculate performance metrics for abstracts
        initial_ids = read_ids_from_file(pmid_file)
        gpt_selected_ids = read_ids_from_file(output_file_abstracts)
        goldstandard_selected_ids = read_ids_from_file("Goldstandard_Selected.txt")  # Assuming the gold standard file is located here
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial_ids, gpt_selected_ids, goldstandard_selected_ids)
        print("Performance metrics for abstracts:")
        print(f"Sensitivity: {sensitivity:.4f}")
        print(f"Specificity: {specificity:.4f}")
        print(f"Positive Predictive Value (PPV): {PPV:.4f}")
        print(f"Negative Predictive Value (NPV): {NPV:.4f}")
        print(f"True Positives: {tp}")
        print(f"True Negatives: {tn}")
        print(f"False Positives: {fp}")
        print(f"False Negatives: {fn}")

    delete_all_data(db_name)

if __name__ == "__main__":
    main()
