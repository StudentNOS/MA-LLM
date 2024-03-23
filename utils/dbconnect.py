import sqlite3
import json

ENSURE="ensure.sqlite"

def init_db(db=ENSURE):
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS meta_analysis (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT,
                title TEXT,
                abstract TEXT,
                fulltext_pdf TEXT,
                short_summary TEXT,
                search_terms TEXT,
                inclusion_criteria TEXT,
                exclusion_criteria TEXT,
                n_studies_query INTEGER,
                n_studies_included INTEGER,
                patients TEXT,
                interventions TEXT,
                controls TEXT,
                study_designs TEXT
            );
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS study (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT,
                meta_analysis_id INTEGER FOREIGNKEY REFERENCES meta_analysis(id),
                title TEXT,
                abstract TEXT,
                fulltext_pdf TEXT,
                total_patients INTEGER,
                int_n_patients INTEGER,
                int_age_central_tendency REAL,
                int_age_central_tendency_type TEXT,
                int_age_dispersion REAL,
                int_age_dispersion_type TEXT,
                int_perc_male REAL,
                int_outcome_metric REAL,
                con_n_patients INTEGER,
                con_age_central_tendency REAL,
                con_age_central_tendency_type TEXT,
                con_age_dispersion REAL,
                con_age_dispersion_type TEXT,
                con_perc_male REAL,
                con_outcome_metric REAL,
                rob2_overall INTEGER,
                rob2_random_sequence INTEGER,
                rob2_deviation INTEGER,
                rob2_missing_outcome INTEGER,
                rob2_measure_out INTEGER,
                rob2_selection INTEGER,
                grade_overall INTEGER
        );
    """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS llm_study_extraction (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                study_id INTEGER FOREIGNKEY REFERENCES study(id),
                total_patients INTEGER,
                int_n_patients INTEGER,
                int_age_central_tendency REAL,
                int_age_central_tendency_type TEXT,
                int_age_dispersion REAL,
                int_age_dispersion_type TEXT,
                int_perc_male REAL,
                int_outcome_metric REAL,
                con_n_patients INTEGER,
                con_age_central_tendency REAL,
                con_age_central_tendency_type TEXT,
                con_age_dispersion REAL,
                con_age_dispersion_type TEXT,
                con_perc_male REAL,
                con_outcome_metric REAL,
                rob2_overall INTEGER,
                rob2_random_sequence INTEGER,
                rob2_deviation INTEGER,
                rob2_missing_outcome INTEGER,
                rob2_measure_out INTEGER,
                rob2_selection INTEGER,
                grade_overall INTEGER
        );
    """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS llm_ma_generation (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                meta_analysis_id INTEGER FOREIGNKEY REFERENCES meta_analysis(id),
                title text,
                fulltext_pdf TEXT,
                abstract TEXT,
                table_1 TEXT
        );
    """)
        # New commands to add selection columns
    try:
        cursor.execute("""
            ALTER TABLE meta_analysis ADD COLUMN selected_title BOOLEAN DEFAULT FALSE;
        """)
    except sqlite3.OperationalError:
        pass  # Ignore error if column already exists

    try:
        cursor.execute("""
            ALTER TABLE meta_analysis ADD COLUMN selected_abstract BOOLEAN DEFAULT FALSE;
        """)
    except sqlite3.OperationalError:
        pass  # Ignore error if column already exists
        conn.commit()

def execute_query(query, db=ENSURE, params=None):
    """
    Executes a given SQL query with optional parameters and returns the fetched results.
    """
    try:
        with sqlite3.connect(db) as conn:
            cursor = conn.cursor()
            if params:
                cursor.execute(query, params)
            else:
                # Splitting and executing each statement if multiple statements are passed without params
                for statement in query.split(";"):
                    statement = statement.strip()
                    if statement:
                        cursor.execute(statement)
            conn.commit()
            return cursor.fetchall()
    except sqlite3.Error as e:
        print(f"An error occurred: {e}")
        return []
    
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
    keys = ", ".join(data.keys())
    values = ", ".join(f"'{value}'" for value in data.values())
    query = f"INSERT INTO {table} ({keys}) VALUES ({values})"
    exec(query, db)

def select(table, columns, where=None, db=ENSURE):
    columns = ", ".join(columns)
    query = f"SELECT {columns} FROM {table}"
    if where:
        query += f" WHERE {where}"
    return exec(query, db)

def create(db=ENSURE):
    input("WARNING: This will delete the existing database and create a new one. Press enter to continue.")
    query = """
        DROP TABLE IF EXISTS meta_analysis; 
        DROP TABLE IF EXISTS study; 
        DROP TABLE IF EXISTS llm_study_extraction; 
        DROP TABLE IF EXISTS llm_ma_generation;
    """
    exec(query, db)
    init_db(db)

# Main function to initialize the database
def main():
    db_name = ENSURE
    init_db(db_name)
    print("Database initialized with the necessary tables.")

if __name__ == "__main__":
    main()
