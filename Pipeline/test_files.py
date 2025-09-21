#!/usr/bin/env python3
"""
Simple test script to debug file loading issues.
This script only loads the input files and prints debug information.
"""

import os
import pandas as pd

def test_file_loading():
    print("=== Testing File Loading ===")

    # Test Initial.txt
    print("\n1. Testing Initial.txt")
    try:
        with open('Initial.txt', 'r', encoding='utf-8') as f:
            initial_content = f.read()
        print("Initial.txt content length:", len(initial_content))
        initial_lines = [line.strip() for line in initial_content.splitlines() if line.strip()]
        print("Number of PMIDs in Initial.txt:", len(initial_lines))
        print("First 5 PMIDs:", initial_lines[:5])
        print("Last 5 PMIDs:", initial_lines[-5:])
        print("Any empty lines?", any(line == '' for line in initial_lines))
    except Exception as e:
        print("ERROR loading Initial.txt:", str(e))

    # Test Goldstandard_Selected.txt
    print("\n2. Testing Goldstandard_Selected.txt")
    try:
        with open('Goldstandard_Selected.txt', 'r', encoding='utf-8') as f:
            goldstandard_content = f.read()
        print("Goldstandard_Selected.txt content length:", len(goldstandard_content))
        goldstandard_lines = [line.strip() for line in goldstandard_content.splitlines() if line.strip()]
        print("Number of PMIDs in Goldstandard_Selected.txt:", len(goldstandard_lines))
        print("First 5 PMIDs:", goldstandard_lines[:5])
        print("Last 5 PMIDs:", goldstandard_lines[-5:])
        print("Any empty lines?", any(line == '' for line in goldstandard_lines))
    except Exception as e:
        print("ERROR loading Goldstandard_Selected.txt:", str(e))

    # Test Prompts.xlsx
    print("\n3. Testing Prompts.xlsx")
    try:
        prompts_df = pd.read_excel('Prompts.xlsx')
        print("Prompts.xlsx shape:", prompts_df.shape)
        print("Columns:", list(prompts_df.columns))
        print("Data types:")
        print(prompts_df.dtypes)
        print("\nFirst few rows:")
        print(prompts_df.head().to_string())

        # Check specific columns
        if 'screen_titles' in prompts_df.columns:
            print("\nscreen_titles unique values:", prompts_df['screen_titles'].unique())
            print("screen_titles value counts:", prompts_df['screen_titles'].value_counts().to_dict())
        else:
            print("ERROR: screen_titles column not found!")

        if 'screen_abstracts' in prompts_df.columns:
            print("screen_abstracts unique values:", prompts_df['screen_abstracts'].unique())
            print("screen_abstracts value counts:", prompts_df['screen_abstracts'].value_counts().to_dict())
        else:
            print("ERROR: screen_abstracts column not found!")

    except Exception as e:
        print("ERROR loading Prompts.xlsx:", str(e))
        import traceback
        traceback.print_exc()

    # Test if files exist
    print("\n4. File existence check:")
    files_to_check = ['Initial.txt', 'Goldstandard_Selected.txt', 'Prompts.xlsx']
    for filename in files_to_check:
        exists = os.path.exists(filename)
        size = os.path.getsize(filename) if exists else 0
        print(f"{filename}: {'EXISTS' if exists else 'MISSING'} (size: {size} bytes)")

    print("\n=== File Loading Test Complete ===")

if __name__ == '__main__':
    test_file_loading()
