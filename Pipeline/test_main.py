#!/usr/bin/env python3
"""
Test script to run the main screening process directly.
"""

import os
import pandas as pd
import sys
sys.path.append('.')

# Import the main functions
from screening_updated import process_prompts, AIModelClient

def test_main_process():
    print("=== Testing Main Screening Process ===")

    # Load test data
    try:
        with open('Initial.txt', 'r', encoding='utf-8') as f:
            initial_pmids = [line.strip() for line in f.read().splitlines() if line.strip()]

        with open('Goldstandard_Selected.txt', 'r', encoding='utf-8') as f:
            goldstandard_pmids = [line.strip() for line in f.read().splitlines() if line.strip()]

        prompts_df = pd.read_excel('Prompts.xlsx')

        print(f"Loaded {len(initial_pmids)} initial PMIDs")
        print(f"Loaded {len(goldstandard_pmids)} goldstandard PMIDs")
        print(f"Loaded prompts DataFrame with shape: {prompts_df.shape}")

        # Only use first 2 prompts for testing
        test_prompts_df = prompts_df.head(2)
        print(f"Using first 2 prompts for testing")

        # Create a mock AI client (this will fail, but we'll see the debug output)
        try:
            ai_client = AIModelClient('groq', 'llama3-8b-8192', 'test_key')
            print("AI client created successfully")
        except Exception as e:
            print(f"AI client creation failed (expected): {str(e)}")
            return

        # Run the process (this will show all debug output)
        print("\n=== Starting Process ===")
        process_prompts(test_prompts_df, initial_pmids[:10], goldstandard_pmids, ai_client)  # Only first 10 PMIDs for testing

    except Exception as e:
        print(f"ERROR in main process: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    test_main_process()
