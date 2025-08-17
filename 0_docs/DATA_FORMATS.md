# Data Format Specifications for ENSURE

This document provides specifications for all data formats used in the ENSURE systematic review automation pipeline.

## Overview

The ENSURE pipeline processes various types of data throughout the systematic review automation workflow. This document standardizes all data formats to ensure consistency, reproducibility, and interoperability across different components of the system.

## Input Data Formats

### 1. Initial PMID Lists

**Purpose**: Contains PubMed IDs for papers to be screened  
**File Pattern**: `Initial.txt`  
**Location**: Each gold standard folder (`1_goldstandards/*/`) and pipeline directory (`3_pipeline/`)

**Format Specification**:
```
File Type: Plain text (.txt)
Encoding: UTF-8
Line Endings: Unix (LF) preferred, Windows (CRLF) accepted
Structure: One PMID per line
```

**Example**:
```
32420649
22901346
34025474
25687662
35500900
27725683
29160598
```

**Validation Rules**:
- Each line must contain exactly one PMID
- PMIDs must be 7-9 digit numeric strings
- No headers, comments, or additional text
- Empty lines are ignored during processing
- Trailing/leading whitespace is trimmed

**Quality Checks**:
```python
def validate_pmid_file(filepath):
    """Validate PMID file format."""
    valid_pmids = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            if not line.isdigit() or len(line) < 7 or len(line) > 9:
                print(f"Warning: Invalid PMID on line {line_num}: {line}")
                continue
            valid_pmids.append(line)
    return valid_pmids
```

### 2. Gold Standard Selections

**Purpose**: Contains PMIDs that human reviewers determined as relevant  
**File Pattern**: `Goldstandard_Selected.txt`  
**Location**: Each gold standard folder (`1_goldstandards/*/`)

**Format Specification**: Identical to Initial PMID Lists

**Example**:
```
32420649
34025474
35500900
```

**Requirements**:
- Must be a subset of corresponding `Initial.txt`
- Represents ground truth for performance evaluation
- Created through manual expert review process
- Should include inter-rater agreement documentation

**Validation**:
```python
def validate_goldstandard(initial_file, goldstandard_file):
    """Ensure goldstandard is subset of initial PMIDs."""
    initial_pmids = set(read_pmids_from_file(initial_file))
    gold_pmids = set(read_pmids_from_file(goldstandard_file))
    
    invalid_pmids = gold_pmids - initial_pmids
    if invalid_pmids:
        raise ValueError(f"Goldstandard contains PMIDs not in initial set: {invalid_pmids}")
    
    return len(gold_pmids), len(initial_pmids)
```

### 3. Prompt Configuration Files

**Purpose**: Configuration for screening prompts and parameters  
**File Pattern**: `Prompts.xlsx`  
**Location**: Each gold standard folder and pipeline directory

**Format Specification**:
```
File Type: Excel Workbook (.xlsx)
Encoding: Default Excel encoding
Structure: Tabular data with defined columns
```

**Required Columns**:

| Column Name | Data Type | Description | Valid Values | Example |
|-------------|-----------|-------------|--------------|---------|
| `screen_titles` | Integer | Enable title screening | 0, 1 | 1 |
| `screen_abstracts` | Integer | Enable abstract screening | 0, 1 | 1 |
| `TitlePrompt` | Text | Prompt text for title screening | String | "Screen the following titles for relevance to diabetes treatment..." |
| `AbstractPrompt` | Text | Prompt text for abstract screening | String | "Evaluate these abstracts for diabetes interventions..." |
| `section_reference` | Text | Type of reference provided | "No Reference", "Topic", "Methods", "Results" | "Topic" |
| `meta_analysis_info` | Text | Meta-analysis context | "No Info", "For MA", "Title", "Criteria" | "For MA" |
| `words` | Integer | Word count of prompt | Positive integer | 75 |

**Optional Columns** (added during processing):
- `MA`: Meta-analysis identifier
- `prompt_id`: Unique prompt identifier  
- Performance metrics (automatically populated)

**Example Structure**:
```excel
screen_titles | screen_abstracts | TitlePrompt | AbstractPrompt | section_reference | meta_analysis_info | words
1            | 1                | "Screen..." | "Evaluate..."  | "Topic"          | "For MA"          | 75
0            | 1                | ""          | "Review..."    | "No Reference"   | "Criteria"        | 52
1            | 0                | "Select..." | ""             | "Methods"        | "Title"           | 68
```

**Validation Rules**:
```python
def validate_prompts_excel(filepath):
    """Validate prompt configuration file."""
    required_columns = [
        'screen_titles', 'screen_abstracts', 'TitlePrompt', 
        'AbstractPrompt', 'section_reference', 'meta_analysis_info', 'words'
    ]
    
    df = pd.read_excel(filepath)
    
    # Check required columns
    missing_cols = set(required_columns) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Validate screen flags
    for col in ['screen_titles', 'screen_abstracts']:
        if not df[col].isin([0, 1]).all():
            raise ValueError(f"Column {col} must contain only 0 or 1")
    
    # Validate prompt consistency
    title_screen_rows = df['screen_titles'] == 1
    if title_screen_rows.any() and df.loc[title_screen_rows, 'TitlePrompt'].isna().any():
        raise ValueError("TitlePrompt required when screen_titles=1")
    
    return len(df)
```

## Output Data Formats

### 1. Selected PMID Files

**Purpose**: Contains PMIDs selected by LLM during screening  
**File Patterns**: 
- `GPT_Selected_Titles.txt` (title screening results)
- `GPT_Selected_Abstracts.txt` (abstract screening results)
**Location**: Pipeline directory (`3_pipeline/`)

**Format Specification**: Identical to Initial PMID Lists

**Example**:
```
32420649
34025474
27725683
29160598
```

**Metadata**:
- Generated automatically by pipeline
- Contains subset of initial PMIDs
- Used for performance evaluation against gold standard

### 2. PubMed Data Export

**Purpose**: Retrieved paper metadata from PubMed API  
**File Pattern**: `PubMed_Data.xlsx`  
**Location**: Pipeline directory (`3_pipeline/`)

**Format Specification**:
```
File Type: Excel Workbook (.xlsx)
Structure: One row per paper
Columns: pmid, title, abstract, [additional metadata]
```

**Column Definitions**:

| Column | Data Type | Description | Example |
|--------|-----------|-------------|---------|
| `pmid` | String | PubMed ID | "32420649" |
| `title` | String | Paper title | "Systematic review of diabetes interventions" |
| `abstract` | String | Paper abstract text | "Background: This study examines..." |

**Quality Considerations**:
- Some papers may have missing abstracts
- Titles may contain special characters or formatting
- Text encoding should preserve Unicode characters

### 3. Performance Results

**Purpose**: Updated prompt configuration with calculated performance metrics  
**File Pattern**: `Prompts.xlsx` (updated in-place)  
**Location**: Pipeline directory (`3_pipeline/`)

**Additional Columns Added**:

| Column Name | Data Type | Description | Range |
|-------------|-----------|-------------|-------|
| `sensitivity_titles` | Float | Title screening sensitivity | [0.0, 1.0] |
| `specificity_titles` | Float | Title screening specificity | [0.0, 1.0] |
| `PPV_titles` | Float | Title screening positive predictive value | [0.0, 1.0] |
| `NPV_titles` | Float | Title screening negative predictive value | [0.0, 1.0] |
| `tp_titles` | Integer | True positives (titles) | ≥ 0 |
| `tn_titles` | Integer | True negatives (titles) | ≥ 0 |
| `fp_titles` | Integer | False positives (titles) | ≥ 0 |
| `fn_titles` | Integer | False negatives (titles) | ≥ 0 |
| `sensitivity_abstracts` | Float | Abstract screening sensitivity | [0.0, 1.0] |
| `specificity_abstracts` | Float | Abstract screening specificity | [0.0, 1.0] |
| `PPV_abstracts` | Float | Abstract screening positive predictive value | [0.0, 1.0] |
| `NPV_abstracts` | Float | Abstract screening negative predictive value | [0.0, 1.0] |
| `tp_abstracts` | Integer | True positives (abstracts) | ≥ 0 |
| `tn_abstracts` | Integer | True negatives (abstracts) | ≥ 0 |
| `fp_abstracts` | Integer | False positives (abstracts) | ≥ 0 |
| `fn_abstracts` | Integer | False negatives (abstracts) | ≥ 0 |
| `Output` | Text | Raw LLM responses for debugging | String |

**Validation Constraints**:
- `tp + tn + fp + fn = total_papers`
- `sensitivity = tp / (tp + fn)`
- `specificity = tn / (tn + fp)`
- `PPV = tp / (tp + fp)`
- `NPV = tn / (tn + fn)`

## Database Schema

### SQLite Database (`ensure.sqlite`)

**Purpose**: Temporary storage during pipeline execution  
**Location**: Pipeline directory (`3_pipeline/`)

**Table Schemas**:

#### Initial Table
```sql
CREATE TABLE Initial (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL UNIQUE,
    title TEXT,
    abstract TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_initial_pmid ON Initial(pmid);
```

#### Titles Table
```sql
CREATE TABLE titles (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL UNIQUE,
    title TEXT,
    abstract TEXT,
    selected_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (pmid) REFERENCES Initial(pmid)
);

CREATE INDEX idx_titles_pmid ON titles(pmid);
```

#### Abstracts Table
```sql
CREATE TABLE abstracts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL UNIQUE,
    title TEXT,
    abstract TEXT,
    selected_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (pmid) REFERENCES titles(pmid)
);

CREATE INDEX idx_abstracts_pmid ON abstracts(pmid);
```

**Data Flow**:
1. Initial papers → `Initial` table
2. Title screening results → `titles` table
3. Abstract screening results → `abstracts` table

## Prompt Engineering Data Formats

### 1. CSV Files in Prompt Engineering Pipeline

**Location**: `2_promptengineering/`

#### Initial Prompts (`initial_prompts.csv`)
```
Format: CSV with semicolon (;) delimiter
Encoding: UTF-8
Headers: Yes
```

**Column Structure**:
```csv
Prompt;section_reference;meta_analysis_info;words;additional_columns
"Screen the following titles...";Topic;For MA;75;...
"Evaluate these abstracts...";No Reference;Criteria;52;...
```

#### Generated Prompts (`generated_prompts.csv`)
- **Source**: Markov chain generation
- **Format**: Same as initial_prompts.csv
- **Purpose**: Synthetic prompt variations

#### LLM Enhanced Prompts (`generated_prompts_llm.csv`)
- **Source**: Large language model enhancement
- **Format**: Same as initial_prompts.csv  
- **Purpose**: Semantically improved prompts

#### Complete Prompts (`complete_prompts.csv`)
- **Source**: Combined dataset
- **Format**: Same as initial_prompts.csv
- **Purpose**: All prompts for analysis

#### Final Prompts (`final_prompts.csv`)
- **Source**: Optimized selection (1,440 prompts)
- **Format**: Same as initial_prompts.csv
- **Purpose**: Production-ready prompt set

### 2. Jupyter Notebook Outputs

**File**: `prompt_engineering.ipynb`  
**Outputs**: Various intermediate data files and visualizations

## Statistical Analysis Data Formats

### 1. R Analysis Input

**File**: `ENSURE_Data.xlsx`  
**Location**: `4_statisticalanalyses/`

**Purpose**: Processed data for statistical analysis in R

**Required Columns**:
- All prompt configuration columns
- All performance metric columns
- `f_score_2_1`: Calculated F-score with 2:1 sensitivity weighting

**Calculated Metrics**:
```r
# F-score calculation with 2:1 sensitivity weighting
f_score_2_1 <- ifelse(screen_titles == 1,
  (1 + 2^2) * (PPV_titles * sensitivity_titles) / ((2^2 * PPV_titles) + sensitivity_titles),
  (1 + 2^2) * (PPV_abstracts * sensitivity_abstracts) / ((2^2 * PPV_abstracts) + sensitivity_abstracts)
)
```

### 2. R Output Formats

**Generated Files**:
- `ensure_analysis.html`: Rendered analysis report
- Various plot files (`.png`, `.pdf`)
- Statistical model outputs

## Data Validation and Quality Control

### 1. Input Validation Functions

```python
def validate_pmid_format(pmid):
    """Validate individual PMID format."""
    return (
        isinstance(pmid, str) and 
        pmid.isdigit() and 
        7 <= len(pmid) <= 9
    )

def validate_file_encoding(filepath):
    """Ensure file is UTF-8 encoded."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            f.read()
        return True
    except UnicodeDecodeError:
        return False

def validate_excel_format(filepath):
    """Validate Excel file can be read."""
    try:
        pd.read_excel(filepath)
        return True
    except Exception:
        return False
```

### 2. Data Integrity Checks

```python
def check_data_consistency():
    """Comprehensive data consistency checks."""
    checks = {
        'pmid_format': validate_all_pmids(),
        'file_encoding': check_all_file_encodings(),
        'goldstandard_subset': verify_goldstandard_subsets(),
        'performance_math': validate_performance_calculations(),
        'database_integrity': check_database_constraints()
    }
    return checks
```

### 3. Error Recovery Procedures

**Missing Data Handling**:
- PMIDs not found in PubMed: Log and continue
- Missing abstracts: Use title-only screening
- Corrupted files: Attempt recovery or alert user

**Data Cleaning Pipeline**:
```python
def clean_pmid_list(pmids):
    """Clean and validate PMID list."""
    cleaned = []
    for pmid in pmids:
        pmid = str(pmid).strip()
        if validate_pmid_format(pmid):
            cleaned.append(pmid)
        else:
            logging.warning(f"Skipping invalid PMID: {pmid}")
    return cleaned
```

## File Organization Best Practices

### 1. Directory Structure Standards

```
project_root/
├── 1_goldstandards/
│   ├── 1. 38957929/
│   │   ├── Initial.txt           # Original PMID list
│   │   ├── Goldstandard_Selected.txt  # Manual selections
│   │   ├── Prompts.xlsx          # Prompt configuration
│   │   └── [supporting files]
│   └── [other meta-analyses]/
├── 2_promptengineering/
│   ├── *.csv                     # Prompt datasets
│   └── prompt_engineering.ipynb  # Development notebook
├── 3_pipeline/
│   ├── Initial.txt               # Working PMID list
│   ├── Goldstandard_Selected.txt # Working gold standard
│   ├── Prompts.xlsx              # Working configuration
│   ├── GPT_Selected_*.txt        # Results
│   ├── PubMed_Data.xlsx          # Retrieved data
│   └── ensure.sqlite             # Working database
└── 4_statisticalanalyses/
    ├── ENSURE_Data.xlsx          # Analysis input
    └── [analysis outputs]
```

### 2. Naming Conventions

**File Naming**:
- Use descriptive, consistent names
- Include processing stage indicators
- Use standardized extensions
- Avoid spaces and special characters

**Examples**:
- `Initial.txt` (not `initial pmids.txt`)
- `GPT_Selected_Titles.txt` (not `gpt_output_titles.txt`)
- `ENSURE_Data.xlsx` (not `data for analysis.xlsx`)

### 3. Version Control Considerations

**Include in Version Control**:
- Template/example files
- Documentation
- Configuration templates
- Empty directory structure

**Exclude from Version Control**:
- Large data files (`.sqlite`, large `.xlsx`)
- Generated outputs
- Temporary files
- Sensitive configuration files

## Performance and Scalability

### 1. Large Dataset Handling

**Batch Processing**:
```python
def process_in_batches(data, batch_size=1000):
    """Process large datasets in manageable batches."""
    for i in range(0, len(data), batch_size):
        batch = data[i:i + batch_size]
        yield process_batch(batch)
```

**Memory Management**:
- Use iterators for large file processing
- Clear intermediate variables
- Monitor memory usage during processing

### 2. Database Optimization

**Indexing Strategy**:
```sql
-- Performance indexes for common queries
CREATE INDEX idx_initial_pmid ON Initial(pmid);
CREATE INDEX idx_titles_pmid ON titles(pmid);
CREATE INDEX idx_abstracts_pmid ON abstracts(pmid);

-- Composite indexes for complex queries
CREATE INDEX idx_initial_title_search ON Initial(pmid, title);
```

**Query Optimization**:
- Use parameterized queries
- Implement connection pooling
- Use transactions for bulk operations

## Security and Privacy Considerations

### 1. Data Sensitivity

**PubMed Data**: Generally public, but check usage terms
**API Keys**: Never include in data files or version control
**Institutional Data**: Follow local data governance policies

### 2. Data Protection

```python
# Example: Anonymization for sharing
def anonymize_data_for_sharing(data):
    """Remove sensitive information for data sharing."""
    anonymous_data = data.copy()
    # Remove or hash any identifying information
    # Keep only necessary columns for reproduction
    return anonymous_data[['pmid', 'performance_metrics']]
```

## Data Format Migration

### 1. Version Compatibility

When updating data formats:
1. Document all changes in CHANGELOG
2. Provide migration scripts
3. Maintain backward compatibility when possible
4. Include format version identifiers

### 2. Migration Scripts

```python
def migrate_prompts_v1_to_v2(old_file, new_file):
    """Migrate prompt file from v1 to v2 format."""
    # Read old format
    old_data = pd.read_excel(old_file)
    
    # Transform to new format
    new_data = transform_prompt_structure(old_data)
    
    # Save in new format
    new_data.to_excel(new_file, index=False)
    
    # Validate migration
    validate_prompts_excel(new_file)
```

This comprehensive data format specification ensures consistent, reliable, and reproducible data handling throughout the ENSURE pipeline, supporting both current usage and future development.
