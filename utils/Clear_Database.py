import sqlite3
from Bio import Entrez
from dbconnect import insert, execute_query, ENSURE  # Importing from your dbconnect script
import pandas as pd  # Import pandas for data manipulation


def delete_all_data(db_path):
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        # List all tables you want to clear
        tables_to_clear = ['Initial', 'titles', 'abstracts', 'full_texts']
        for table in tables_to_clear:
            cursor.execute(f"DELETE FROM {table}")
        print("All tables have been cleared.")
        conn.commit()

def main():
    db_path = ENSURE  # Ensure this is your database path
    delete_all_data(db_path)

if __name__ == "__main__":
    main()
