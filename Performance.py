import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
from dbconnect import ENSURE

def read_ids_from_file(file_path):
    with open(file_path, 'r') as file:
        ids = {line.strip() for line in file if line.strip()}
    return ids

def calculate_performance_metrics(initial, gpt_selected, goldstandard_selected):
    initial_set = set(initial)
    gs_selected_set = set(goldstandard_selected)
    gpt_selected_set = set(gpt_selected)

    tp = len(gs_selected_set.intersection(gpt_selected_set))
    fn = len(gs_selected_set.difference(gpt_selected_set))
    fp = len(gpt_selected_set.difference(gs_selected_set))
    tn = len(initial_set.difference(gs_selected_set.union(gpt_selected_set)))

    sensitivity = tp / (tp + fn) if tp + fn else 0
    specificity = tn / (tn + fp) if tn + fp else 0
    PPV = tp / (tp + fp) if tp + fp else 0
    NPV = tn / (fn + tn) if fn + tn else 0

    return sensitivity, specificity, PPV, NPV, tp, tn, fp, fn

def create_performance_table(tp, tn, fp, fn, sensitivity, specificity, PPV, NPV, filename):
    # Table data
    cell_text = [
        [f"True Positive\n(TP)\n{tp}", f"False Positive\n(FP)\n{fp}", f"Positive\n Predictive\n Value\n{PPV:.2f}"],
        [f"False Negative\n(FN)\n{fn}", f"True Negative\n(TN)\n{tn}", f"Negative\n Predictive\n Value\n{NPV:.2f}"],
        [f"Sensitivity\n{sensitivity:.2f}", f"Specificity\n{specificity:.2f}", ""]
    ]
    
    # Cell colors
    cell_colors = [
        ["#FFDDDD", "#DDFFDD", "#FFFFDD"],
        ["#DDFFDD", "#FFDDDD", "#DDFFFF"],
        ["#DDFFFF", "#FFFFDD", "#FFFFFF"]
    ]

    # Create a figure and a single subplot
    fig, ax = plt.subplots()
    
    # Hide the axes
    ax.axis('off')

    # Add a table at the bottom of the axes
    the_table = ax.table(cellText=cell_text,
                         cellColours=cell_colors,
                         colWidths=[0.3]*3,
                         loc='center',
                         cellLoc='center')
    
    # Scaling the table
    the_table.scale(1, 10)

    # Save the figure
    plt.savefig(filename, bbox_inches='tight', pad_inches=0.05)
    # Table data
    cell_text = [
        [f"True Positive\n(TP)\n{tp}", f"False Positive\n(FP)\n{fp}", f"Positive\n Predictive\n Value\n{PPV:.2f}"],
        [f"False Negative\n(FN)\n{fn}", f"True Negative\n(TN)\n{tn}", f"Negative\n Predictive\n Value\n{NPV:.2f}"],
        [f"Sensitivity\n{sensitivity:.2f}", f"Specificity\n{specificity:.2f}", ""]
    ]
    
    # Cell colors
    cell_colors = [
        ["#FFDDDD", "#DDFFDD", "#FFFFDD"],
        ["#DDFFDD", "#FFDDDD", "#DDFFFF"],
        ["#DDFFFF", "#FFFFDD", "#FFFFFF"]
    ]

    # Create a figure and a single subplot
    fig, ax = plt.subplots()
    
    # Hide the axes
    ax.axis('off')

    # Add a table at the bottom of the axes
    the_table = ax.table(cellText=cell_text,
                         cellColours=cell_colors,
                         colWidths=[0.3]*3,
                         loc='center',
                         cellLoc='center')
    
    # Scaling the table
    the_table.scale(1, 10)

    # Save the figure
    plt.savefig(filename, bbox_inches='tight', pad_inches=0.05)

def fetch_data_for_report_and_chart():
    """
    Fetches data from the database.
    """
    with sqlite3.connect(ENSURE) as conn:
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
