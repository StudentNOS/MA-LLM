# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 10:57:15 2025

@author: Till
"""

# Read IDs from a file into a set
def read_ids(file_path):
    with open(file_path, 'r') as file:
        return set(line.strip() for line in file)

# Main function
def check_ids(initial_file, goldstandard_file):
    initial_ids = read_ids(initial_file)
    goldstandard_ids = read_ids(goldstandard_file)

    # Find missing IDs
    missing_ids = goldstandard_ids - initial_ids

    if not missing_ids:
        print("All IDs match")
    else:
        print("Missing IDs:", *missing_ids)

# Specify file paths
initial_file = "Initial.txt"
goldstandard_file = "Goldstandard_Selected.txt"

# Run the script
check_ids(initial_file, goldstandard_file)