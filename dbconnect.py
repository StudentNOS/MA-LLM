from PATHS import ENSURE

import sqlite3
import json

def init_db(db=ENSURE):
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        
        # Drop existing tables if they exist
        cursor.execute("DROP TABLE IF EXISTS Initial;")
        cursor.execute("DROP TABLE IF EXISTS titles;")
        cursor.execute("DROP TABLE IF EXISTS abstracts;")
        cursor.execute("DROP TABLE IF EXISTS full_texts;")
        
        # Create the Initial table
        cursor.execute("""
            CREATE TABLE Initial (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT,
                title TEXT,
                abstract TEXT,
                fulltext_pdf TEXT DEFAULT '',
                short_summary TEXT DEFAULT '',
                search_terms TEXT DEFAULT '',
                inclusion_criteria TEXT DEFAULT '',
                exclusion_criteria TEXT DEFAULT '',
                n_studies_query INTEGER DEFAULT 0,
                n_studies_included INTEGER DEFAULT 0,
                patients TEXT DEFAULT '',
                interventions TEXT DEFAULT '',
                controls TEXT DEFAULT '',
                study_designs TEXT DEFAULT ''
            );
        """)

        # Create the titles table
        cursor.execute("""
            CREATE TABLE titles (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT,
                title TEXT,
                abstract TEXT,
                fulltext_pdf TEXT DEFAULT '',
                short_summary TEXT DEFAULT '',
                search_terms TEXT DEFAULT '',
                inclusion_criteria TEXT DEFAULT '',
                exclusion_criteria TEXT DEFAULT '',
                n_studies_query INTEGER DEFAULT 0,
                n_studies_included INTEGER DEFAULT 0,
                patients TEXT DEFAULT '',
                interventions TEXT DEFAULT '',
                controls TEXT DEFAULT '',
                study_designs TEXT DEFAULT ''
            );
        """)

        # Create the abstracts table
        cursor.execute("""
            CREATE TABLE abstracts (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT,
                title TEXT,
                abstract TEXT,
                fulltext_pdf TEXT DEFAULT '',
                short_summary TEXT DEFAULT '',
                search_terms TEXT DEFAULT '',
                inclusion_criteria TEXT DEFAULT '',
                exclusion_criteria TEXT DEFAULT '',
                n_studies_query INTEGER DEFAULT 0,
                n_studies_included INTEGER DEFAULT 0,
                patients TEXT DEFAULT '',
                interventions TEXT DEFAULT '',
                controls TEXT DEFAULT '',
                study_designs TEXT DEFAULT ''
            );
        """)

        # Create the full_texts table
        cursor.execute("""
            CREATE TABLE full_texts (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT,
                title TEXT,
                abstract TEXT,
                fulltext_pdf TEXT DEFAULT '',
                short_summary TEXT DEFAULT '',
                search_terms TEXT DEFAULT '',
                inclusion_criteria TEXT DEFAULT '',
                exclusion_criteria TEXT DEFAULT '',
                n_studies_query INTEGER DEFAULT 0,
                n_studies_included INTEGER DEFAULT 0,
                patients TEXT DEFAULT '',
                interventions TEXT DEFAULT '',
                controls TEXT DEFAULT '',
                study_designs TEXT DEFAULT ''
            );
        """)

        # Ensure the tables are created
        conn.commit()



def execute_query(query, db=ENSURE, params=None):
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        if params:
            cursor.execute(query, params)
        else:
            cursor.execute(query)
        conn.commit()
        return cursor.fetchall()

    
def insert(table, key_dict, db=ENSURE, serialize_json=False):
    """
    Inserts data into a specified table. Can serialize complex data types if needed.
    """
    try:
        with sqlite3.connect(db) as conn:
            cursor = conn.cursor()
            serialized_values = []
            for value in key_dict.values():
                if serialize_json and isinstance(value, (dict, list)):
                    serialized_values.append(json.dumps(value))
                else:
                    serialized_values.append(value)

            keys = ", ".join(key_dict.keys())
            placeholders = ", ".join(["?"] * len(key_dict))
            query = f"INSERT INTO {table} ({keys}) VALUES ({placeholders})"
            cursor.execute(query, serialized_values)
            conn.commit()
    except sqlite3.Error as e:
        print(f"An error occurred: {e}")
        

def insert2(table, data, db=ENSURE):
    """
    Alternative insert function using different syntax.
    """
    keys = ", ".join(data.keys())
    values = ", ".join(f"'{value}'" for value in data.values())
    query = f"INSERT INTO {table} ({keys}) VALUES ({values})"
    execute_query(query, db)

def select(table, columns, where=None, db=ENSURE):
    """
    Selects data from a specified table with optional where clause.
    """
    columns = ", ".join(columns)
    query = f"SELECT {columns} FROM {table}"
    if where:
        query += f" WHERE {where}"
    return execute_query(query, db)

def create(db=ENSURE):
    """
    Drops existing tables and recreates the database.
    """
    input("WARNING: This will delete the existing database and create a new one. Press enter to continue.")
    query = """
        DROP TABLE IF EXISTS Initial; 
        DROP TABLE IF EXISTS titles; 
        DROP TABLE IF EXISTS abstracts; 
        DROP TABLE IF EXISTS full_texts;
    """
    execute_query(query, db)
    init_db(db)

def delete_all_data(db_name):
    with sqlite3.connect(db_name) as conn:
        cursor = conn.cursor()
        # List all tables you want to clear
        tables_to_clear = ['Initial', 'titles', 'abstracts', 'full_texts']
        for table in tables_to_clear:
            cursor.execute(f"DELETE FROM {table}")
        print("All tables have been cleared.")
        conn.commit()
