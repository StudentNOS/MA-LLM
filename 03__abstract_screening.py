# Input: Search terms
# Output: Pre-selected papers based on abstracts
from utils import dbconnect as db, ncbi
import pandas as pd
PROMPT = """
"""

def get_abstracts(search_terms):
    """Get abstracts from papers based on search terms
    Returns:
        list: List of abstracts
    """
    abstr_list = ncbi.search(search_terms)
    column_names = ncbi.FIELDS
    df = pd.DataFrame(abstr_list, columns=column_names)
    df.rename(columns={'pubmedid': 'pmid'}, inplace=True)
    # Append multiple fields into the abstract column
    df['abstract'] = df["abstract"] + df['objective'] + df['background'] + df['methods'] + df['results'] + df['conclusions']
    
    db.insert('study', df[["pmid","title", "abstract"]])
    
    

def exclude_abstracts(abstracts):
    """Exclude abstracts based on exclusion criteria
    Returns:
        list: List of abstracts
    """
    

