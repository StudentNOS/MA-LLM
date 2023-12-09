import sqlite3
import json

def db_insert(table, key_dict, db="aied.sqlite"):
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
        
        
def db_exec(query, db="aied.sqlite"):
    with sqlite3.connect(db) as conn:
        c = conn.cursor()
        c.execute(query)
        conn.commit()
        return c.fetchall()

def db_insert2(table, data, db="aied.sqlite"):
    keys = ", ".join(data.keys())
    values = ", ".join(f"'{value}'" for value in data.values())
    query = f"INSERT INTO {table} ({keys}) VALUES ({values})"
    db_exec(query, db)
