# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 17:02:52 2024

@author: tillj
"""

import pandas as pd
from numpy import nan
import os
#import fitz
from dbconnect import insert, delete_all_data, ENSURE
from PubMed import count_papers, create_excel_from_db, read_pmids_from_file, fetch_details 
from GPT import get_data_in_batches, generate_prompt, screen_with_openai, screen_pdf_with_openai, match_data_to_ids, save_pmids_to_file, move_records, save_fulltext_pmids_to_file # for LLama models import GPT_forLLama
from Performance import read_ids_from_file, calculate_performance_metrics, create_performance_table, fetch_data_for_report_and_chart, create_excel_report, draw_plot_chart
import ssl

ssl._create_default_https_context = ssl._create_unverified_context
script_dir = os.path.dirname(os.path.abspath(__file__))

def main(screen_titles, TitlePrompt, screen_abstracts, AbstractPrompt, row_index, Prompts):
    # Step 1: Read PMIDs and fetch details
    pmid_file = "Initial.txt"
    pmids = read_pmids_from_file(pmid_file)
    records = fetch_details(pmids)
    create_excel_from_db()
    
    # Step 2: Insert fetched records into the 'Initial' table
    for article in records.get("PubmedArticle", []):
        pmid = article["MedlineCitation"]["PMID"]
        title = article["MedlineCitation"]["Article"].get("ArticleTitle", "")
        abstract_text = ""
        if "Abstract" in article["MedlineCitation"]["Article"]:
            abstract_text = " ".join(str(abstract) for abstract in article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"])
        
        insert("Initial", {"pmid": pmid, "title": title, "abstract": abstract_text}, ENSURE, False)
    
    # Step 3: Screen titles
    if screen_titles == 1:
        all_screened_titles = []
        all_matched_pmids_titles = []

        for titles_batch in get_data_in_batches("titles"):
            prompt = generate_prompt(titles_batch, "titles", TitlePrompt)
            screened_title_initials = screen_with_openai(prompt)
            matched_pmids_titles = match_data_to_ids(screened_title_initials, titles_batch, "titles")
            
            all_matched_pmids_titles.extend(matched_pmids_titles)
            all_screened_titles.extend(screened_title_initials)

        move_records(all_matched_pmids_titles, "titles")
        output_file_titles = os.path.join(script_dir, 'GPT_Selected_Titles.txt')
        save_pmids_to_file(all_matched_pmids_titles, output_file_titles)
        
        # Calculate performance metrics for titles
        initial_ids = read_ids_from_file(pmid_file)
        gpt_selected_ids = read_ids_from_file(output_file_titles)
        goldstandard_selected_ids = read_ids_from_file(os.path.join(script_dir, 'Goldstandard_Selected.txt'))
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial_ids, gpt_selected_ids, goldstandard_selected_ids)
        
        # Update Prompts DataFrame
        Prompts.at[row_index, "sensitivity_titles"] = sensitivity
        Prompts.at[row_index, "specificity_titles"] = specificity
        Prompts.at[row_index, "PPV_titles"] = PPV
        Prompts.at[row_index, "NPV_titles"] = NPV
        Prompts.at[row_index, "tp_titles"] = tp
        Prompts.at[row_index, "tn_titles"] = tn
        Prompts.at[row_index, "fp_titles"] = fp
        Prompts.at[row_index, "fn_titles"] = fn
    else:
        # Insert NaN for skipped title screening
        for col in ["sensitivity_titles", "specificity_titles", "PPV_titles", "NPV_titles", "tp_titles", "tn_titles", "fp_titles", "fn_titles"]:
            Prompts.at[row_index, col] = nan
    
    # Step 5: Screen abstracts
    if screen_abstracts == 1:
        all_screened_abstracts = []
        all_matched_pmids_abstracts = []

        for abstracts_batch in get_data_in_batches("abstracts"):
            prompt = generate_prompt(abstracts_batch, "abstracts", AbstractPrompt)
            screened_abstract_initials = screen_with_openai(prompt)
            matched_pmids_abstracts = match_data_to_ids(screened_abstract_initials, abstracts_batch, "abstracts")
            
            all_matched_pmids_abstracts.extend(matched_pmids_abstracts)
            all_screened_abstracts.extend(screened_abstract_initials)

        move_records(all_matched_pmids_abstracts, "abstracts")
        output_file_abstracts = os.path.join(script_dir, 'GPT_Selected_Abstracts.txt')
        save_pmids_to_file(all_matched_pmids_abstracts, output_file_abstracts)
        
        # Calculate performance metrics for abstracts
        initial_ids = read_ids_from_file(pmid_file)
        gpt_selected_ids = read_ids_from_file(output_file_abstracts)
        goldstandard_selected_ids = read_ids_from_file(os.path.join(script_dir, 'Goldstandard_Selected.txt'))
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial_ids, gpt_selected_ids, goldstandard_selected_ids)
        
        # Update Prompts DataFrame
        Prompts.at[row_index, "sensitivity_abstracts"] = sensitivity
        Prompts.at[row_index, "specificity_abstracts"] = specificity
        Prompts.at[row_index, "PPV_abstracts"] = PPV
        Prompts.at[row_index, "NPV_abstracts"] = NPV
        Prompts.at[row_index, "tp_abstracts"] = tp
        Prompts.at[row_index, "tn_abstracts"] = tn
        Prompts.at[row_index, "fp_abstracts"] = fp
        Prompts.at[row_index, "fn_abstracts"] = fn
    else:
        # Insert NaN for skipped abstract screening
        for col in ["sensitivity_abstracts", "specificity_abstracts", "PPV_abstracts", "NPV_abstracts", "tp_abstracts", "tn_abstracts", "fp_abstracts", "fn_abstracts"]:
            Prompts.at[row_index, col] = nan
    
    delete_all_data(ENSURE)

if __name__ == "__main__":
    Prompts = pd.read_excel("Prompts.xlsx")
    
    for Index, CurrentPrompt in Prompts.iterrows():
        print(Index)
        main(CurrentPrompt["screen_titles"], CurrentPrompt["TitlePrompt"], 
             CurrentPrompt["screen_abstracts"], CurrentPrompt["AbstractPrompt"], Index, Prompts)
    
    # Save updated DataFrame back to the same Excel file
    Prompts.to_excel("Prompts.xlsx", index=False)
