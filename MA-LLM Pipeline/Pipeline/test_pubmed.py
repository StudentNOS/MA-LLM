#!/usr/bin/env python3
"""
Test script to debug PubMed API and article fetching issues.
"""

import os
import sqlite3
import pandas as pd
from Bio import Entrez

# Test PubMed API connection
def test_pubmed_connection():
    print("=== Testing PubMed Connection ===")

    # Set email (required by NCBI)
    Entrez.email = "test@example.com"
    print("Entrez email set to: test@example.com")

    # Test with a simple query
    try:
        print("\n1. Testing simple PubMed search...")
        handle = Entrez.esearch(db="pubmed", term="cancer", retmax=5)
        results = Entrez.read(handle)
        handle.close()
        print(f"Search successful! Found {results['Count']} articles")
        print(f"Sample PMIDs: {results['IdList'][:3]}")
    except Exception as e:
        print(f"ERROR in PubMed search: {str(e)}")
        return False

    return True

# Test article fetching
def test_article_fetching():
    print("\n=== Testing Article Fetching ===")

    # Use PMIDs from Initial.txt
    try:
        with open('Initial.txt', 'r', encoding='utf-8') as f:
            initial_pmids = [line.strip() for line in f.read().splitlines() if line.strip()]

        print(f"Loaded {len(initial_pmids)} PMIDs from Initial.txt")
        print(f"First 5 PMIDs: {initial_pmids[:5]}")

        # Test fetching first 3 articles
        test_pmids = initial_pmids[:3]
        print(f"\nFetching details for PMIDs: {test_pmids}")

        Entrez.email = "test@example.com"
        handle = Entrez.efetch(db="pubmed", id=",".join(test_pmids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        print(f"Successfully fetched {len(records.get('PubmedArticle', []))} articles")

        # Print details of first article
        if records.get('PubmedArticle'):
            article = records['PubmedArticle'][0]
            pmid = article["MedlineCitation"]["PMID"]
            title = article["MedlineCitation"]["Article"].get("ArticleTitle", "No title")

            abstract = ""
            if "Abstract" in article["MedlineCitation"]["Article"]:
                abstract = " ".join(
                    str(t) for t in article["MedlineCitation"]["Article"]["Abstract"].get("AbstractText", [])
                )

            print("
First article details:")
            print(f"PMID: {pmid}")
            print(f"Title: {title[:100]}...")
            print(f"Abstract length: {len(abstract)} characters")
            print(f"Abstract preview: {abstract[:200]}..." if abstract else "No abstract")

    except Exception as e:
        print(f"ERROR in article fetching: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

    return True

# Test database operations
def test_database():
    print("\n=== Testing Database Operations ===")

    db_file = 'test_debug.db'

    try:
        # Create test database
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()

        # Create table
        cursor.execute("""
            CREATE TABLE Articles (
                id INTEGER PRIMARY
