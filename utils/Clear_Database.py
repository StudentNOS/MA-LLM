import sqlite3
from Bio import Entrez
from dbconnect import insert, execute_query, ENSURE  # Importing from your dbconnect script
import pandas as pd  # Import pandas for data manipulation


def delete_all_abstracts(ENSURE):
    with sqlite3.connect(ENSURE) as conn:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM meta_analysis")
        conn.commit()
        print("Study table has been cleared.")
def main():
    delete_all_abstracts(ENSURE)
if __name__ == "__main__":
    main()
