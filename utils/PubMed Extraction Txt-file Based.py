import sqlite3
from Bio import Entrez
from dbconnect import insert, execute_query, ENSURE  # Importing from your dbconnect script
import pandas as pd  # Import pandas for data manipulation

# Replace 'your_email' with your actual email address
Entrez.email = "Mail"

def count_papers():
    """
    Counts the number of papers in the database.
    """
    query = "SELECT COUNT(*) FROM Initial"
    try:
        count = execute_query(query, ENSURE)[0][0]
        return count
    except IndexError:
        return 0

def create_excel_from_db():
    with sqlite3.connect(ENSURE) as conn:
        df = pd.read_sql_query("SELECT pmid, title, abstract FROM Initial", conn)
        df.to_excel("PubMed_Data.xlsx", index=False)
        print("Excel file 'PubMed_Data.xlsx' created.")

def read_pmids_from_file(file_path):
    with open(file_path, "r") as file:
        return [line.strip() for line in file if line.strip().isdigit()]

def fetch_details(pmids):
    ids = ','.join(pmids)
    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

def main():
    pmids = read_pmids_from_file("Initial.txt")
    records = fetch_details(pmids)
    
    for article in records["PubmedArticle"]:
        pmid = article["MedlineCitation"]["PMID"]
        title = article["MedlineCitation"]["Article"].get("ArticleTitle", "")
        abstract_text = ""
        if "Abstract" in article["MedlineCitation"]["Article"]:
            abstract_text = " ".join([str(abstract) for abstract in article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]])
        
        insert("Initial", {"pmid": pmid, "title": title, "abstract": abstract_text}, ENSURE, False)
        print(f"Inserted: {pmid}")
    
    paper_count = count_papers()
    print(f"Total number of papers in the database: {paper_count}")
    
    create_excel_from_db()

if __name__ == "__main__":
    main()
