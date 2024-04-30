import matplotlib.pyplot as plt

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
    """
    Generates a table with performance metrics, arranged similarly to the uploaded image, and saves it as an image.
    """
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

def main():
    # Part 1: Calculate overlap percentages
    initial_ids = read_ids_from_file("Initial.txt")
    gpt_selected_ids = read_ids_from_file("GPT_Selected.txt")
    goldstandard_selected_ids = read_ids_from_file("Goldstandard_Selected.txt")
    
    # Calculate performance metrics
    sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
        initial_ids, gpt_selected_ids, goldstandard_selected_ids)
    
    create_performance_table(tp, tn, fp, fn, sensitivity, specificity, PPV, NPV, "performance_table.png")
    
    print(f"Sensitivity: {sensitivity:.4f}")
    print(f"Specificity: {specificity:.4f}")
    print(f"Positive Predictive Value (PPV): {PPV:.4f}")
    print(f"Negative Predictive Value (NPV): {NPV:.4f}")
    print(f"True Positives: {tp}")
    print(f"True Negatives: {tn}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")

if __name__ == "__main__":
    main()
