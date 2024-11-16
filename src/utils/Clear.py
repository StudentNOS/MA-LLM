import sqlite3
from src.utils.dbconnect import ENSURE


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
    delete_all_data(ENSURE)

if __name__ == "__main__":
    main()
