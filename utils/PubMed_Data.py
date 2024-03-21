import sqlite3
from Bio import Entrez
from dbconnect import insert, execute_query, ENSURE  # Importing from your dbconnect script
import pandas as pd  # Import pandas for data manipulation

# Replace 'your_email' with your actual email address
Entrez.email = "Entrez.Email"

def search_pubmed(search_term, start_year, end_year):
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=1000, mindate=start_year, maxdate=end_year, datetype="pdat")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

# Function to count the number of papers in the database
def count_papers():
    """
    Counts the number of papers in the database.
    """
    query = "SELECT COUNT(*) FROM meta_analysis"  # Ensure this matches your table name
    try:
        count = execute_query(query, ENSURE)[0][0]
        return count
    except IndexError:
        return 0

# Function to create an Excel table from the SQLite database
def create_excel_from_db():
    with sqlite3.connect(ENSURE) as conn:
        df = pd.read_sql_query("SELECT pmid, title, abstract FROM meta_analysis", conn)  # Adjust query as needed
        df.to_excel("PubMed_Data.xlsx", index=False)  # Writing to an Excel file
        print("Excel file 'PubMed_Data.xlsx' created.")

def main():
    search_term = input("Enter PubMed search term: ")
    start_year = input("Enter start year (Lower Bound): ")
    end_year = input("Enter end year (Upper Bound): ")
    
    pmids = search_pubmed(search_term, start_year, end_year)
    records = fetch_details(pmids)
    
    for article in records["PubmedArticle"]:
        pmid = article["MedlineCitation"]["PMID"]
        article_data = article["MedlineCitation"]["Article"]
        title = article_data.get("ArticleTitle", "")
        
        abstract_text = ""
        if "Abstract" in article_data:
            abstract_text = " ".join([str(abstract) for abstract in article_data["Abstract"]["AbstractText"]])
        
        insert("meta_analysis", {"pmid": pmid, "title": title, "abstract": abstract_text}, ENSURE, False)
        print(f"Inserted: {pmid}")

    paper_count = count_papers()
    print(f"Total number of papers in the database: {paper_count}")

    create_excel_from_db()

if __name__ == "__main__":
    main()
