import sqlite3
import json

def exec(query, db="ensure.sqlite"):
    with sqlite3.connect(db) as conn:
        c = conn.cursor()
        for statement in query.split(";"):
            statement = statement.strip()
            if not statement:
                continue
            c.execute(statement)
        conn.commit()
        return c.fetchall()
    
def insert(table, key_dict, db="ensure.sqlite"):
    # Serialize complex data types to JSON
    serialized_values = []
    for value in key_dict.values():
        if isinstance(value, (dict, list)):
            # Serialize dictionaries and lists to JSON
            serialized_values.append(json.dumps(value))
        else:
            # Use other values directly
            serialized_values.append(value)

    # Prepare the keys and placeholders for the query
    keys = ", ".join(key_dict.keys())
    placeholders = ", ".join(["?"] * len(key_dict))

    # Create the INSERT query with placeholders
    query = f"INSERT INTO {table} ({keys}) VALUES ({placeholders})"

    # Execute the query with the serialized values
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        cursor.execute(query, serialized_values)
        conn.commit()
        

def insert2(table, data, db="ensure.sqlite"):
    keys = ", ".join(data.keys())
    values = ", ".join(f"'{value}'" for value in data.values())
    query = f"INSERT INTO {table} ({keys}) VALUES ({values})"
    exec(query, db)

def select(table, columns, where=None, db="ensure.sqlite"):
    columns = ", ".join(columns)
    query = f"SELECT {columns} FROM {table}"
    if where:
        query += f" WHERE {where}"
    return exec(query, db)

def create(db="ensure.sqlite"):
    input("WARNING: This will delete the existing database and create a new one. Press enter to continue.")
    query = "DROP TABLE IF EXISTS meta_analysis; DROP TABLE IF EXISTS study; DROP TABLE IF EXISTS llm_study_extraction; DROP TABLE IF EXISTS llm_ma_generation;"
    query += """
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
        grade_overall INTEGER);
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
    CREATE TABLE IF NOT EXISTS llm_ma_generation (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        meta_analysis_id INTEGER FOREIGNKEY REFERENCES meta_analysis(id),
        title text,
        fulltext_pdf TEXT,
        abstract TEXT,
        table_1 TEXT
        );
        """
    exec(query, db)