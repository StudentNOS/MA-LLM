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
    # True Negatives (TN) are those IDs in the initial set that are neither in GPT selected nor in Goldstandard selected
    tn = len(initial_set.difference(gs_selected_set.union(gpt_selected_set)))

    sensitivity = tp / (tp + fn) if tp + fn else 0
    specificity = tn / (tn + fp) if tn + fp else 0

    return sensitivity, specificity, tp, tn, fp, fn

def main():
    # Part 1: Calculate overlap percentages
    initial_ids = read_ids_from_file("Initial.txt")
    gpt_selected_ids = read_ids_from_file("GPT_Selected.txt")
    goldstandard_selected_ids = read_ids_from_file("Goldstandard_Selected.txt")
    
    # Calculate performance metrics
    sensitivity, specificity, tp, tn, fp, fn = calculate_performance_metrics(
        initial_ids, gpt_selected_ids, goldstandard_selected_ids)
    
    print(f"Sensitivity: {sensitivity:.4f}")
    print(f"Specificity: {specificity:.4f}")
    print(f"True Positives: {tp}")
    print(f"True Negatives: {tn}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")

if __name__ == "__main__":
    main()
