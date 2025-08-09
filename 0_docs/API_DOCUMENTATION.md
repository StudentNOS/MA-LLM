# ENSURE API Documentation

This document provides comprehensive API documentation for all modules in the ENSURE systematic review automation pipeline.

## Module Overview

The ENSURE pipeline consists of the following main modules:

- **ENSURE.py**: Main pipeline orchestrator
- **PubMed.py**: PubMed data retrieval and management
- **GPT.py**: Large Language Model interface for screening
- **Performance.py**: Performance metrics calculation and visualization
- **dbconnect.py**: Database connection and management
- **Scramble.py**: Data randomization utilities

## Core API Reference

### ENSURE.py

Main pipeline orchestrator that coordinates the entire screening workflow.

#### `main(screen_titles, TitlePrompt, screen_abstracts, AbstractPrompt, row_index, Prompts)`

Executes the complete screening pipeline for a single prompt configuration.

**Parameters:**
- `screen_titles` (int): Flag to enable title screening (1=yes, 0=no)
- `TitlePrompt` (str): Prompt text for title screening
- `screen_abstracts` (int): Flag to enable abstract screening (1=yes, 0=no)
- `AbstractPrompt` (str): Prompt text for abstract screening
- `row_index` (int): Index of current prompt in DataFrame
- `Prompts` (DataFrame): Pandas DataFrame containing prompt configurations

**Returns:**
- None (modifies Prompts DataFrame in-place)

**Side Effects:**
- Updates database with screening results
- Creates output files with selected PMIDs
- Calculates and stores performance metrics

**Example:**
```python
import pandas as pd
from ENSURE import main

# Load prompt configuration
prompts_df = pd.read_excel("Prompts.xlsx")

# Execute pipeline for first prompt
main(
    screen_titles=1,
    TitlePrompt="Screen titles for diabetes research...",
    screen_abstracts=1, 
    AbstractPrompt="Evaluate abstracts for diabetes interventions...",
    row_index=0,
    Prompts=prompts_df
)
```

---

### PubMed.py

Handles all PubMed data retrieval and database export functionality.

#### `read_pmids_from_file(file_path)`

Reads PubMed IDs from a text file.

**Parameters:**
- `file_path` (str): Path to file containing PMIDs (one per line)

**Returns:**
- `list`: List of PMID strings

**Raises:**
- `FileNotFoundError`: If specified file doesn't exist
- `ValueError`: If file contains invalid PMID formats

**Example:**
```python
from PubMed import read_pmids_from_file

pmids = read_pmids_from_file("Initial.txt")
print(f"Loaded {len(pmids)} PMIDs")
```

#### `fetch_details(pmids)`

Fetches detailed metadata for a list of PMIDs from PubMed.

**Parameters:**
- `pmids` (list): List of PMID strings

**Returns:**
- `dict`: Entrez record dictionary containing paper metadata

**Raises:**
- `HTTPError`: If PubMed API request fails
- `ValueError`: If invalid PMIDs provided

**Example:**
```python
from PubMed import fetch_details

pmids = ["32420649", "22901346"]
records = fetch_details(pmids)

for article in records.get("PubmedArticle", []):
    pmid = article["MedlineCitation"]["PMID"]
    title = article["MedlineCitation"]["Article"].get("ArticleTitle", "")
    print(f"{pmid}: {title}")
```

#### `count_papers()`

Returns the count of papers currently in the Initial database table.

**Returns:**
- `int`: Number of papers in database

**Example:**
```python
from PubMed import count_papers

paper_count = count_papers()
print(f"Database contains {paper_count} papers")
```

#### `create_excel_from_db()`

Exports all papers from the Initial database table to an Excel file.

**Side Effects:**
- Creates "PubMed_Data.xlsx" in current directory

**Example:**
```python
from PubMed import create_excel_from_db

create_excel_from_db()
print("Excel export completed")
```

---

### GPT.py

Interfaces with Large Language Models for automated paper screening.

#### `get_data_in_batches(decision, batch_size=None)`

Generator that yields batches of papers from the database for processing.

**Parameters:**
- `decision` (str): Type of screening ("titles" or "abstracts")
- `batch_size` (int, optional): Number of papers per batch (defaults: 50 for titles, 6 for abstracts)

**Yields:**
- `list`: List of tuples containing (pmid, content) pairs

**Raises:**
- `ValueError`: If invalid decision parameter provided

**Example:**
```python
from GPT import get_data_in_batches

# Process titles in batches
for batch in get_data_in_batches("titles", batch_size=25):
    print(f"Processing batch of {len(batch)} titles")
    for pmid, title in batch:
        print(f"{pmid}: {title[:50]}...")
```

#### `generate_prompt(data, decision, manual)`

Generates a complete prompt for LLM screening including data and formatting instructions.

**Parameters:**
- `data` (list): List of tuples containing (pmid, content) pairs
- `decision` (str): Type of screening ("titles" or "abstracts")
- `manual` (str): Custom prompt text provided by user

**Returns:**
- `str`: Complete formatted prompt for LLM

**Example:**
```python
from GPT import generate_prompt

data = [("12345", "Example paper title"), ("67890", "Another paper title")]
manual_prompt = "Select papers about diabetes treatment"

prompt = generate_prompt(data, "titles", manual_prompt)
print(prompt)
```

#### `screen_with_openai(prompt)`

Sends a prompt to the LLM and returns the screening results.

**Parameters:**
- `prompt` (str): Complete prompt text for screening

**Returns:**
- `list`: List of selected PMID strings

**Raises:**
- `APIError`: If LLM API request fails
- `ValueError`: If response cannot be parsed

**Example:**
```python
from GPT import screen_with_openai

prompt = "Screen the following titles for diabetes research..."
selected_pmids = screen_with_openai(prompt)
print(f"LLM selected {len(selected_pmids)} papers")
```

#### `match_data_to_ids(screened_results, original_data, decision)`

Matches LLM screening results back to original PMIDs.

**Parameters:**
- `screened_results` (list): Results from LLM screening
- `original_data` (list): Original data sent to LLM
- `decision` (str): Type of screening performed

**Returns:**
- `list`: List of matched PMID strings

**Example:**
```python
from GPT import match_data_to_ids

screened = ["12345", "67890"]
original = [("12345", "Title 1"), ("67890", "Title 2"), ("11111", "Title 3")]

matched_pmids = match_data_to_ids(screened, original, "titles")
print(f"Matched {len(matched_pmids)} PMIDs")
```

#### `save_pmids_to_file(pmids, filepath)`

Saves a list of PMIDs to a text file.

**Parameters:**
- `pmids` (list): List of PMID strings
- `filepath` (str): Output file path

**Side Effects:**
- Creates or overwrites specified file

**Example:**
```python
from GPT import save_pmids_to_file

selected_pmids = ["12345", "67890", "11111"]
save_pmids_to_file(selected_pmids, "GPT_Selected_Titles.txt")
```

#### `move_records(pmids, table_name)`

Moves database records from Initial table to specified target table.

**Parameters:**
- `pmids` (list): List of PMID strings to move
- `table_name` (str): Target database table ("titles" or "abstracts")

**Side Effects:**
- Inserts records into target table
- Records remain in Initial table

---

### Performance.py

Calculates and visualizes screening performance metrics.

#### `read_ids_from_file(file_path)`

Reads PMID set from a text file.

**Parameters:**
- `file_path` (str): Path to file containing PMIDs

**Returns:**
- `set`: Set of PMID strings

**Example:**
```python
from Performance import read_ids_from_file

initial_ids = read_ids_from_file("Initial.txt")
selected_ids = read_ids_from_file("GPT_Selected_Titles.txt")
print(f"Selected {len(selected_ids)} out of {len(initial_ids)} papers")
```

#### `calculate_performance_metrics(initial, gpt_selected, goldstandard_selected)`

Calculates comprehensive performance metrics for screening results.

**Parameters:**
- `initial` (set): Set of all initial PMIDs
- `gpt_selected` (set): Set of PMIDs selected by LLM
- `goldstandard_selected` (set): Set of PMIDs in gold standard

**Returns:**
- `tuple`: (sensitivity, specificity, PPV, NPV, tp, tn, fp, fn)
  - `sensitivity` (float): True positive rate
  - `specificity` (float): True negative rate  
  - `PPV` (float): Positive predictive value
  - `NPV` (float): Negative predictive value
  - `tp` (int): True positives count
  - `tn` (int): True negatives count
  - `fp` (int): False positives count
  - `fn` (int): False negatives count

**Example:**
```python
from Performance import calculate_performance_metrics, read_ids_from_file

initial = read_ids_from_file("Initial.txt")
gpt_selected = read_ids_from_file("GPT_Selected_Titles.txt")
goldstandard = read_ids_from_file("Goldstandard_Selected.txt")

sensitivity, specificity, ppv, npv, tp, tn, fp, fn = calculate_performance_metrics(
    initial, gpt_selected, goldstandard
)

print(f"Sensitivity: {sensitivity:.3f}")
print(f"Specificity: {specificity:.3f}")
print(f"PPV: {ppv:.3f}")
print(f"NPV: {npv:.3f}")
```

#### `create_performance_table(tp, tn, fp, fn, sensitivity, specificity, PPV, NPV, filename)`

Creates and saves a visual confusion matrix table.

**Parameters:**
- `tp`, `tn`, `fp`, `fn` (int): Confusion matrix counts
- `sensitivity`, `specificity`, `PPV`, `NPV` (float): Performance metrics
- `filename` (str): Output image filename

**Side Effects:**
- Saves confusion matrix visualization as image file

#### `create_excel_report(data, file_path)`

Creates Excel report from screening progress data.

**Parameters:**
- `data` (list): List of tuples containing screening progress
- `file_path` (str): Output Excel file path

**Side Effects:**
- Creates Excel file with screening progress report

#### `draw_plot_chart(data)`

Creates and displays screening progression visualization.

**Parameters:**
- `data` (list): List of tuples containing screening progress

**Side Effects:**
- Displays plot and saves as "Screening_Progression.png"

---

### dbconnect.py

Handles all database operations and connectivity.

#### `insert(table, data, database, create_if_missing=True)`

Inserts data into specified database table.

**Parameters:**
- `table` (str): Target table name
- `data` (dict): Dictionary of column-value pairs
- `database` (str): Database file path
- `create_if_missing` (bool): Whether to create table if it doesn't exist

**Side Effects:**
- Inserts record into database
- Creates table if necessary and create_if_missing=True

**Example:**
```python
from dbconnect import insert, ENSURE

data = {
    "pmid": "12345678",
    "title": "Example paper title",
    "abstract": "Example abstract text"
}

insert("Initial", data, ENSURE)
```

#### `execute_query(query, database)`

Executes SQL query and returns results.

**Parameters:**
- `query` (str): SQL query string
- `database` (str): Database file path

**Returns:**
- `list`: List of tuples containing query results

**Example:**
```python
from dbconnect import execute_query, ENSURE

results = execute_query("SELECT COUNT(*) FROM Initial", ENSURE)
count = results[0][0]
print(f"Database contains {count} records")
```

#### `delete_all_data(database)`

Deletes all data from all tables in the database.

**Parameters:**
- `database` (str): Database file path

**Side Effects:**
- Removes all records from all tables
- Preserves table structure

**Warning:**
This operation is irreversible. Use with caution.

---

### Scramble.py

Utilities for data randomization and anonymization.

#### `scramble()`

Randomly shuffles the order of PMIDs in Initial.txt file.

**Side Effects:**
- Modifies Initial.txt file in current directory
- Maintains same PMIDs but in random order

**Example:**
```python
from Scramble import scramble

# Randomize order of papers
scramble()
print("PMID order randomized")
```

## Configuration

### Environment Variables

The pipeline uses the following configuration:

**PATHS.py**:
```python
api_key = "your_openai_api_key"      # OpenAI API key
email = "your_email@domain.com"       # Email for PubMed access
```

### Database Configuration

**Database File**: `ensure.sqlite` (SQLite database)

**Tables**:
- `Initial`: All papers from initial PMID list
- `titles`: Papers selected during title screening
- `abstracts`: Papers selected during abstract screening

## Error Handling

### Common Exceptions

- **APIError**: LLM API request failures
- **HTTPError**: PubMed API connection issues
- **FileNotFoundError**: Missing input files
- **SQLiteError**: Database operation failures
- **ValueError**: Invalid input parameters

### Error Recovery

```python
try:
    results = screen_with_openai(prompt)
except APIError as e:
    print(f"API error: {e}")
    # Implement retry logic
except ValueError as e:
    print(f"Invalid prompt format: {e}")
    # Handle prompt formatting issues
```

## Performance Considerations

### Batch Processing
- Title screening: 50 papers per batch (recommended)
- Abstract screening: 6 papers per batch (due to length)
- Adjust batch sizes based on API limits and paper length

### Rate Limiting
```python
import time

# Add delays between API calls
time.sleep(1.0)  # 1 second between requests

# For large batches
time.sleep(300)  # 5 minutes between large operations
```

### Memory Management
```python
# Clear variables after large operations
del large_dataset
import gc; gc.collect()

# Process in chunks for very large datasets
chunk_size = 1000
for i in range(0, len(data), chunk_size):
    chunk = data[i:i+chunk_size]
    process_chunk(chunk)
```

## Testing

### Unit Tests

```python
import unittest
from Performance import calculate_performance_metrics

class TestPerformance(unittest.TestCase):
    def test_perfect_classification(self):
        initial = {"1", "2", "3", "4"}
        selected = {"1", "2"}
        goldstandard = {"1", "2"}
        
        sens, spec, ppv, npv, tp, tn, fp, fn = calculate_performance_metrics(
            initial, selected, goldstandard
        )
        
        self.assertEqual(sens, 1.0)  # Perfect sensitivity
        self.assertEqual(spec, 1.0)  # Perfect specificity
        self.assertEqual(ppv, 1.0)   # Perfect precision
```

### Integration Tests

```python
def test_pipeline_integration():
    # Test complete pipeline with small dataset
    pmids = ["32420649", "22901346"]
    
    # Write test input
    with open("test_initial.txt", "w") as f:
        f.write("\n".join(pmids))
    
    # Run pipeline components
    records = fetch_details(pmids)
    assert len(records["PubmedArticle"]) == 2
    
    # Clean up
    os.remove("test_initial.txt")
```

## Best Practices

### API Usage
- Always handle rate limits gracefully
- Implement exponential backoff for retries
- Monitor API usage and costs
- Use batch processing to minimize requests

### Data Management
- Validate input data before processing
- Backup databases before major operations  
- Use transactions for multi-step database operations
- Implement data integrity checks

### Code Quality
- Add comprehensive error handling
- Include logging for debugging
- Use type hints for better documentation
- Follow Python naming conventions

This API documentation provides the foundation for understanding and extending the ENSURE pipeline. For specific implementation details, refer to the source code and inline documentation.
