import sqlite3
from dbconnect import ENSURE
DB_PATH = "ensure.sqlite"

def fetch_pmids_from_db(db_path):
    """
    Fetches all PMIDs from the meta_analysis table in the database.
    """
    query = "SELECT pmid FROM meta_analysis"
    pmids = []
    try:
        with sqlite3.connect(db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            pmids = [row[0] for row in cursor.fetchall()]
    except sqlite3.Error as e:
        print(f"An error occurred: {e}")
    return pmids

def create_initial_gpt_file(pmids, filename="GPT_Initial.txt"):
    """
    Creates a text file with one PMID per line.
    """
    with open(filename, "w") as file:
        for pmid in pmids:
            file.write(f"{pmid}\n")
    print(f"File {filename} created with PMIDs.")

def read_ids_from_database():
    with sqlite3.connect(ENSURE) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT pmid FROM meta_analysis")
        rows = cursor.fetchall()
    return {row[0] for row in rows}

def read_ids_from_file(file_path):
    with open(file_path, 'r') as file:
        ids = {line.strip() for line in file if line.strip()}
    return ids

def calculate_overlap_percentage(set1, set2):
    intersection = set1.intersection(set2)
    if not set1:
        return 0 
    percentage_overlap = (len(intersection) / len(set1)) * 100
    return percentage_overlap

def calculate_performance_metrics(goldstandard_initial, goldstandard_selected, gpt_initial, gpt_selected):
    # Convert lists to sets for efficient intersection and difference operations
    gs_initial_set = set(goldstandard_initial)
    gs_selected_set = set(goldstandard_selected)
    gpt_initial_set = set(gpt_initial)
    gpt_selected_set = set(gpt_selected)

    tp = len(gs_selected_set.intersection(gpt_selected_set))
    fp = len(gpt_selected_set.difference(gs_selected_set))
    fn = len(gs_selected_set.difference(gpt_selected_set))
    all_ids = gs_initial_set.union(gpt_initial_set)  # Combine all IDs for a complete set
    tn = len(all_ids) - tp - fp - fn

    sensitivity = tp / (tp + fn) if tp + fn else 0
    specificity = tn / (tn + fp) if tn + fp else 0

    return sensitivity, specificity, tp, tn, fp, fn

def main():
    pmids = fetch_pmids_from_db(DB_PATH)
    create_initial_gpt_file(pmids)

    # Part 1: Calculate overlap percentages
    database_ids = read_ids_from_database()
    file_ids = read_ids_from_file("Goldstandard_Initial.txt")  # Adjust file path as necessary
    db_in_file_percentage = calculate_overlap_percentage(database_ids, file_ids)
    file_in_db_percentage = calculate_overlap_percentage(file_ids, database_ids)
    print(f"Percentage of database IDs in file: {db_in_file_percentage:.2f}%")
    print(f"Percentage of file IDs in database: {file_in_db_percentage:.2f}%")

    # Part 2: Calculate performance metrics
    goldstandard_initial_path = "Goldstandard_Initial.txt"
    goldstandard_selected_path = "Goldstandard_Selected.txt"
    gpt_initial_path = "GPT_Initial.txt"
    gpt_selected_path = "GPT_Selected.txt"
    
    goldstandard_initial = read_ids_from_file(goldstandard_initial_path)
    goldstandard_selected = read_ids_from_file(goldstandard_selected_path)
    gpt_initial = read_ids_from_file(gpt_initial_path)
    gpt_selected = read_ids_from_file(gpt_selected_path)
    
    sensitivity, specificity, tp, tn, fp, fn = calculate_performance_metrics(
        goldstandard_initial, goldstandard_selected, gpt_initial, gpt_selected)
    
    print(f"Sensitivity: {sensitivity:.4f}")
    print(f"Specificity: {specificity:.4f}")
    print(f"True Positives: {tp}")
    print(f"True Negatives: {tn}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")

if __name__ == "__main__":
    main()
