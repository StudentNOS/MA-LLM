import sqlite3
import pandas as pd
import matplotlib.pyplot as plt

DB_PATH = "ensure.sqlite"  # Update this to your actual database path

def fetch_data_for_report_and_chart():
    """
    Fetches data from the database.
    """
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        
        # Fetch initial records
        query_initial = "SELECT pmid FROM Initial"
        cursor.execute(query_initial)
        initial_data = cursor.fetchall()
        
        # Fetch titles records
        query_titles = "SELECT pmid FROM titles"
        cursor.execute(query_titles)
        titles_data = cursor.fetchall()
        
        # Fetch abstracts records
        query_abstracts = "SELECT pmid FROM abstracts"
        cursor.execute(query_abstracts)
        abstracts_data = cursor.fetchall()
        
        # Create a combined dataset
        data = set()
        for pmid in initial_data:
            data.add((pmid[0], None, None))
        
        for pmid in titles_data:
            data.add((pmid[0], pmid[0], None))
        
        for pmid in abstracts_data:
            data.add((pmid[0], pmid[0], pmid[0]))
        
        # Convert set to list and return
        return list(data)

def create_excel_report(data, file_path="Screening_Progress_Report.xlsx"):
    """
    Creates an Excel report from the provided data.
    """
    df = pd.DataFrame(data, columns=["Initial", "titles", "abstracts"])
    df.to_excel(file_path, index=False)
    print(f"Excel report generated: {file_path}")

def draw_plot_chart(data):
    """
    Draws a plot chart based on the provided data.
    """
    # Count the number of papers at each stage
    initial_count = len(data)
    title_screened_count = sum(1 for _, title, _ in data if title is not None)
    abstract_screened_count = sum(1 for _, _, abstract in data if abstract is not None)

    # Setup for drawing the chart
    stages = ['Abstract Screening', 'Title Screening', 'Initial']  # Reversed order
    counts = [abstract_screened_count, title_screened_count, initial_count]  # Reversed order

    fig, ax = plt.subplots()
    ax.barh(stages, counts, color='skyblue')
    ax.set_xlabel('Number of Articles')
    ax.set_title('Screening Progression')

    # Add the text on the bars
    for i, count in enumerate(counts):
        ax.text(count, i, str(count), va='center')

    plt.tight_layout()
    plt.savefig("Screening_Progression.png")
    plt.show()

def main():
    data = fetch_data_for_report_and_chart()
    create_excel_report(data)
    draw_plot_chart(data)

if __name__ == "__main__":
    main()
