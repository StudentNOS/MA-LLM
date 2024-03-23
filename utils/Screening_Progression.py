import sqlite3
import pandas as pd
import matplotlib.pyplot as plt

DB_PATH = "C:/Users/tillj/.spyder-py3/ensure.sqlite"  # Update this to your actual database path

def fetch_data_for_report_and_chart():
    """
    Fetches data from the database.
    """
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        query = """
            SELECT pmid,
                   CASE WHEN selected_title THEN pmid END as title_screened,
                   CASE WHEN selected_abstract THEN pmid END as abstract_screened
            FROM meta_analysis
        """
        cursor.execute(query)
        data = cursor.fetchall()
        return data

def create_excel_report(data, file_path="Screening_Progress_Report.xlsx"):
    """
    Creates an Excel report from the provided data.
    """
    df = pd.DataFrame(data, columns=["Initial", "Title Screening", "Abstract Screening"])
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
