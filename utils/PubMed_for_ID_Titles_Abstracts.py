# Input: Search terms
# Output: Pre-selected papers based on titles

import sqlite3
import requests
import xml.etree.ElementTree as ET
import pandas as pd
from db import insert, execute_query, ENSURE

# Insert data into the 'study' table
def insert_data(ENSURE, data):
    with sqlite3.connect(ENSURE) as conn:
        cursor = conn.cursor()
        cursor.executemany("INSERT OR IGNORE INTO study (pmid, title, abstract) VALUES (?, ?, ?)", data)
        conn.commit()

# Get PubMed IDs based on a search term and a maximum number of results
def get_pubmed_ids(search_term, max_results=5):
    """
    Fetches PubMed IDs based on a search term and a maximum number of results.
    """
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        search_url = f"{base_url}esearch.fcgi?db=pubmed&term={search_term}&retmax={max_results}&retmode=json"
        response = requests.get(search_url)
        data = response.json()
        id_list = data['esearchresult']['idlist']
        return id_list
    except requests.RequestException as e:
        print(f"An error occurred while fetching PubMed IDs: {e}")
        return []

# Fetch papers from PubMed based on a list of PubMed IDs
def fetch_papers(id_list):
    """
    Fetches papers from PubMed based on a list of PubMed IDs.
    """
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        ids = ','.join(id_list)
        fetch_url = f"{base_url}efetch.fcgi?db=pubmed&id={ids}&retmode=xml"
        response = requests.get(fetch_url)
        root = ET.fromstring(response.content)

        papers_data = []
        for article in root.findall('.//PubmedArticle'):
            pubmed_id = article.find('.//PMID').text
            title = article.find('.//ArticleTitle').text
            abstract = article.find('.//Abstract/AbstractText')
            abstract_text = "No abstract available" if abstract is None else abstract.text
            papers_data.append((pubmed_id, title, abstract_text))

        return papers_data
    except requests.RequestException as e:
        print(f"An error occurred while fetching papers: {e}")
        return []

# Function to count the number of papers
def count_papers():
    """
    Counts the number of papers in the database.
    """
    query = "SELECT COUNT(*) FROM study"
    try:
        count = execute_query(query, ENSURE)[0][0]
        return count
    except IndexError:
        return 0

# Main function to fetch, insert, and count paper data
def main():
    """
    Main function to fetch paper data based on a search term, insert it into the database, 
    and then count the number of papers in the database.
    """
    search_term = input("Enter search term: ")
    id_list = get_pubmed_ids(search_term)
    papers_data = fetch_papers(id_list)

    for data in papers_data:
        insert('study', {'pmid': data[0], 'title': data[1], 'abstract': data[2]}, ENSURE)

    paper_count = count_papers()
    print(f"Total number of papers in the database: {paper_count}")

if __name__ == "__main__":
    main()
    

#Create Excel table from SQLite database
#with sqlite3.connect(ENSURE) as conn:
    #df = pd.read_sql_query("SELECT * FROM study", conn)
    #df.to_excel("C:/Users/tillj/Desktop/ensure.xlsx", index=False)

#Function for deletion
#def delete_all_abstracts(ENSURE='ensure.sqlite'):
    #with sqlite3.connect(ENSURE) as conn:
        #cursor = conn.cursor()
        #cursor.execute("DELETE FROM study")
        #conn.commit()
        #print("Study table has been cleared.")
#def main():
    #ENSURE = 'ensure.sqlite'
    #delete_all_abstracts(ENSURE)
#if __name__ == "__main__":
    #main()
