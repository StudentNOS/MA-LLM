# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 17:02:52 2024

@author: tillj
"""

import os
import fitz
from src.utils.dbconnect import insert, delete_all_data, ENSURE
from PubMed import count_papers, create_excel_from_db, read_pmids_from_file, fetch_details 
from GPT import get_data_in_batches, generate_prompt, screen_with_openai, screen_pdf_with_openai, match_data_to_ids, save_pmids_to_file, move_records, save_fulltext_pmids_to_file
from Performance import read_ids_from_file, calculate_performance_metrics, create_performance_table, fetch_data_for_report_and_chart, create_excel_report, draw_plot_chart

script_dir = os.path.dirname(os.path.abspath(__file__))

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
    create_excel_from_db() # For those who love spreadsheets
    
    # Load search term, inclusion, and exclusion criteria from text files
    search_term_file = os.path.join(script_dir, 'search_term.txt')
    inclusion_file = os.path.join(script_dir, 'inclusion.txt')
    exclusion_file = os.path.join(script_dir, 'exclusion.txt')
    
    with open(search_term_file, 'r') as f:
        search_term = f.read().strip()
    
    with open(inclusion_file, 'r') as f:
        inclusion = f.read().strip()
    
    with open(exclusion_file, 'r') as f:
        exclusion = f.read().strip()
    
    # Step 4: Ask user if they want to screen titles
    screen_titles = input("Do you want to screen the titles (y/n)? ").strip().lower() == 'y'
    if screen_titles:
        search_term = search_term
        inclusion = inclusion
        exclusion = exclusion
        manual = input("Prompt: ").strip() # because computers can't read minds (yet...)
        all_screened_titles = []
        all_matched_pmids_titles = []

        # Process titles
        print("Processing titles... (Please wait)")
        for titles_batch in get_data_in_batches("titles"):
            prompt = generate_prompt(search_term, inclusion, exclusion, titles_batch, "titles", manual)
            screened_title_initials = screen_with_openai(prompt)
            matched_pmids_titles = match_data_to_ids(screened_title_initials, titles_batch, "titles")
            
            all_matched_pmids_titles.extend(matched_pmids_titles)
            all_screened_titles.extend(screened_title_initials)

        # Move relevant titles to the 'titles' table
        print("Moving relevant titles to the 'titles' table...")
        move_records(all_matched_pmids_titles, "titles")

        # Save PubMed IDs for titles to a file
        output_file_titles = os.path.join(script_dir, 'GPT_Selected_Titles.txt')
        print(f"Saving PubMed IDs for titles to {output_file_titles}...")
        save_pmids_to_file(all_matched_pmids_titles, output_file_titles)
        
        print(f"PubMed IDs saved to {output_file_titles}")
        print(f"{len(all_matched_pmids_titles)} titles have been marked as selected based on title screening.")
        
        # Step 5: Calculate performance metrics for titles
        initial_ids = read_ids_from_file(pmid_file)
        gpt_selected_ids = read_ids_from_file(output_file_titles)
        goldstandard_selected_ids = read_ids_from_file(os.path.join(script_dir, 'Goldstandard_Selected.txt'))
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial_ids, gpt_selected_ids, goldstandard_selected_ids)
        print("Performance metrics for titles:")
        print("Sensitivity & Specificity:")
        print(f"{sensitivity:.2f}")
        print(f"{specificity:.2f}")
        print(f"Positive Predictive Value (PPV): {PPV:.2f}")
        print(f"Negative Predictive Value (NPV): {NPV:.2f}")
        print(f"True Positives: {tp}")
        print(f"True Negatives: {tn}")
        print(f"False Positives: {fp}")
        print(f"False Negatives: {fn}")

    # Step 6: Ask user if they want to screen abstracts
    screen_abstracts = input("Do you want to screen the abstracts (y/n)? ").strip().lower() == 'y'
    if screen_abstracts:
        search_term = search_term
        inclusion = inclusion
        exclusion = exclusion
        manual = input("Prompt: ").strip()
        all_screened_abstracts = []
        all_matched_pmids_abstracts = []

        # Process abstracts
        print("Processing abstracts... (Please wait)")
        for abstracts_batch in get_data_in_batches("abstracts"):
            prompt = generate_prompt(search_term, inclusion, exclusion, abstracts_batch, "abstracts", manual)
            screened_abstract_initials = screen_with_openai(prompt)
            matched_pmids_abstracts = match_data_to_ids(screened_abstract_initials, abstracts_batch, "abstracts")
            
            all_matched_pmids_abstracts.extend(matched_pmids_abstracts)
            all_screened_abstracts.extend(screened_abstract_initials)

        # Move relevant abstracts to the 'abstracts' table
        print("Moving relevant abstracts to the 'abstracts' table...")
        move_records(all_matched_pmids_abstracts, "abstracts")

        # Save PubMed IDs for abstracts to a file
        output_file_abstracts = os.path.join(script_dir, 'GPT_Selected_Abstracts.txt')
        print(f"Saving PubMed IDs for abstracts to {output_file_abstracts}...")
        save_pmids_to_file(all_matched_pmids_abstracts, output_file_abstracts)
        
        print(f"PubMed IDs saved to {output_file_abstracts}")
        print(f"{len(all_matched_pmids_abstracts)} abstracts have been marked as selected based on abstract screening.")
        
        # Step 7: Calculate performance metrics for abstracts
        initial_ids = read_ids_from_file(pmid_file)
        gpt_selected_ids = read_ids_from_file(output_file_abstracts)
        goldstandard_selected_ids = read_ids_from_file(os.path.join(script_dir, 'Goldstandard_Selected.txt'))
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial_ids, gpt_selected_ids, goldstandard_selected_ids)
        print("Performance metrics for abstracts:")
        print("Sensitivity & Specificity:")
        print(f"{sensitivity:.2f}")
        print(f"{specificity:.2f}")
        print(f"Positive Predictive Value (PPV): {PPV:.2f}")
        print(f"Negative Predictive Value (NPV): {NPV:.2f}")
        print(f"True Positives: {tp}")
        print(f"True Negatives: {tn}")
        print(f"False Positives: {fp}")
        print(f"False Negatives: {fn}")

    # Step 8: Ask user if they want to create a visual overview of screening - Almost there ;)
    create_visual_overview = input("Would you like to create a visual overview for your Screening? (y/n)").strip().lower() == 'y'
    if create_visual_overview:
        # Create performance table
        create_performance_table(tp, tn, fp, fn, sensitivity, specificity, PPV, NPV, "performance_table.png")
        
        # Fetch data for report and chart
        data = fetch_data_for_report_and_chart()
        
        # Create Excel report
        create_excel_report(data)
        
        # Draw plot chart
        draw_plot_chart(data)
        
        print("Plots and overview table have been saved to folder. Either open in folder or click 'Plots' in the console to view them.") # Reward: Visuals!
    
    # Step 9: Fulltext Screening
    fulltexts_folder = "FullTexts"
    if not os.path.exists(fulltexts_folder):
        os.makedirs(fulltexts_folder)
    
    print("Thank you! The title and abstract screenings have been successful. To proceed to the full text screening, please follow these steps:")
    print("1. Open the 'GPT_Selected_Abstracts.txt' file in your folder. This file contains the PubMed IDs of the papers selected up to this point.")
    print("2. Look for these IDs on PubMed and download the papers as PDF files.")
    print('3. Please give each PDF file its PubMed ID as a name (e.g., "27339398.pdf").')
    print('4. Save these PDF files to the folder "FullTexts", which should have been created in the folder where you keep all the scripts.')
    
    proceed_fulltext_screening = input("Have you completed these steps and would like to proceed to the screening? (y/n)").strip().lower() == 'y'
    
    if proceed_fulltext_screening:
        search_term = search_term
        inclusion = inclusion
        exclusion = exclusion
        manual_fulltext = input("Prompt: ").strip()
        relevant_pmids = []
    
        for file_name in os.listdir(fulltexts_folder):
            if file_name.endswith(".pdf"):
                pmid = os.path.splitext(file_name)[0]
                
                with fitz.open(os.path.join(fulltexts_folder, file_name)) as doc:
                    num_pages = doc.page_count
                    chunk_size = 5  # Set the chunk size
                    
                    for start_page in range(0, num_pages, chunk_size):
                        text_chunk = ""
                        end_page = min(start_page + chunk_size, num_pages)
                        
                        for page_num in range(start_page, end_page):
                            page = doc.load_page(page_num)
                            text_chunk += page.get_text()
                        
                        gpt_output = screen_pdf_with_openai(text_chunk, search_term, inclusion, exclusion, manual_fulltext)
                        
                        if "Relevant: " in gpt_output:
                            relevant_pmids.append(pmid)
                            break  # Once a chunk is relevant, mark the full PDF as relevant and skip the rest
    
        if relevant_pmids:
            move_records(relevant_pmids, "full_texts")
            print(f"Moved {len(relevant_pmids)} records to full_texts.")
            
            output_file_fulltexts = os.path.join(script_dir, 'GPT_Selected_Fulltexts.txt')
            save_fulltext_pmids_to_file(ENSURE, output_file_fulltexts)
            print(f"PubMed IDs saved to {output_file_fulltexts}")
        else:
            print("No relevant full texts found.")
    else:
        print("Fulltext screening process halted. Complete the steps and try again.")

    # Step 10: Performance Evaluation for Fulltext Screening
    # Read the IDs from the initial file, fulltext screening output, and gold standard file
    initial_ids = read_ids_from_file(pmid_file)  # Original PMIDs
    gpt_selected_ids = read_ids_from_file(os.path.join(script_dir, 'GPT_Selected_Fulltexts.txt'))
    goldstandard_selected_ids = read_ids_from_file(os.path.join(script_dir, 'Goldstandard_Selected.txt'))  # Gold standard PMIDs
    
    # Calculate performance metrics
    sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
        initial_ids, gpt_selected_ids, goldstandard_selected_ids)
    
    # Print performance metrics
    print("Performance metrics for fulltext screening:")
    print("Sensitivity & Specificity:")
    print(f"{sensitivity:.2f}")
    print(f"{specificity:.2f}")
    print(f"Positive Predictive Value (PPV): {PPV:.2f}")
    print(f"Negative Predictive Value (NPV): {NPV:.2f}")
    print(f"True Positives: {tp}")
    print(f"True Negatives: {tn}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")

    delete_all_data(ENSURE)

if __name__ == "__main__":
    main()
