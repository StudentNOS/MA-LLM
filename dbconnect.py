import os
import sqlite3
import json

# set database path
script_dir = os.path.dirname(os.path.abspath(__file__))
ENSURE = os.path.join(script_dir, 'ensure.sqlite')

# initialize database
def init_db(db=ENSURE):
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()

        # drop tables
        cursor.execute("DROP TABLE IF EXISTS Initial;")
        cursor.execute("DROP TABLE IF EXISTS titles;")
        cursor.execute("DROP TABLE IF EXISTS abstracts;")
        cursor.execute("DROP TABLE IF EXISTS full_texts;")
        
        # create tables
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

        # commit changes
        conn.commit()

# actually create ensure database
init_db(ENSURE)
print("Database initialized with the necessary tables.")

# execute SQL query
def execute_query(query, db=ENSURE, params=None):
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        if params:
            cursor.execute(query, params)
        else:
            cursor.execute(query)
        conn.commit()
        return cursor.fetchall()

# insert data into table
def insert(table, key_dict, db=ENSURE, serialize_json=False):
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

# select data from table
def select(table, columns, where=None, db=ENSURE):
    columns = ", ".join(columns)
    query = f"SELECT {columns} FROM {table}"
    if where:
        query += f" WHERE {where}"
    return execute_query(query, db)

# delete all data from tables
def delete_all_data(db_name):
    with sqlite3.connect(db_name) as conn:
        cursor = conn.cursor()
        tables_to_clear = ['Initial', 'titles', 'abstracts', 'full_texts']
        for table in tables_to_clear:
            cursor.execute(f"DELETE FROM {table}")
        print("All tables have been cleared.")
        conn.commit()
