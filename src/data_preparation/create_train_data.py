import csv

# Read the initial numbers
with open('./data/Initial.txt', 'r') as file:
    initial_numbers = set(file.read().splitlines())

# Read the gold standard numbers
with open('./data/Goldstandard_Selected.txt', 'r') as file:
    gold_standard_numbers = set(file.read().splitlines())

# Create the CSV file
with open('output.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Number', 'Label'])  # Write header

    for number in initial_numbers:
        label = 'yes' if number in gold_standard_numbers else 'no'
        csvwriter.writerow([number, label])

# Count the 'yes' labels in the output file
yes_count = 0
with open('output.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    next(csvreader)  # Skip header
    for row in csvreader:
        if row[1] == 'yes':
            yes_count += 1

# Count the entries in Goldstandard_Selected.txt
gold_standard_count = len(gold_standard_numbers)

# Print if they are the same
if yes_count == gold_standard_count:
    print("The counts are the same.")
else:
    print("The counts are different.")