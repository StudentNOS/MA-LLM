# Data Format Specifications for ENSURE

This document provides detailed specifications for all data formats used in the ENSURE systematic review automation pipeline.

## Input Data Formats

### 1. Initial PMID Lists

**File**: `Initial.txt` (in each meta-analysis folder or pipeline directory)

**Format**: Plain text file with one PubMed ID per line

**Structure**:
```
32420649
22901346
34025474
25687662
35500900
```

**Requirements**:
- Each line contains exactly one PMID (8-digit numeric identifier)
- No headers or additional text
- UTF-8 encoding
- Unix line endings (LF) preferred
- No trailing whitespace

**Validation**:
- PMIDs should be valid and exist in PubMed
- Duplicates will be automatically handled
- Empty lines are ignored

### 2. Gold Standard Selections

**File**: `Goldstandard_Selected.txt`

**Format**: Identical to Initial PMID lists

**Purpose**: Contains PMIDs that human reviewers determined as relevant

**Structure**:
```
32420649
34025474
35500900
```

**Requirements**:
- Must be a subset of PMIDs from corresponding `Initial.txt`
- Used for calculating performance metrics (sensitivity, specificity, etc.)
- Manual curation required

### 3. Prompt Configuration

**File**: `Prompts.xlsx`

**Format**: Excel workbook (.xlsx)

**Required Columns**:

| Column Name | Data Type | Description | Example Values |
|------------|-----------|-------------|----------------|
| `screen_titles` | Integer | Whether to screen titles (1=yes, 0=no) | 1, 0 |
| `screen_abstracts` | Integer | Whether to screen abstracts (1=yes, 0=no) | 1, 0 |
| `TitlePrompt` | Text | Prompt text for title screening | "Screen the following titles..." |
| `AbstractPrompt` | Text | Prompt text for abstract screening | "Evaluate these abstracts..." |
| `section_reference` | Text | Type of reference provided | "No Reference", "Topic", "Methods" |
| `meta_analysis_info` | Text | Meta-analysis context | "No Info", "For MA", "Title", "Criteria" |
| `words` | Integer | Word count of prompt | 50, 100, 150 |

**Optional Columns** (for analysis):
- `MA`: Meta-analysis identifier
- `prompt_id`: Unique prompt identifier
- Performance metrics (automatically filled by pipeline)

**Example Structure**:
```excel
screen_titles | screen_abstracts | TitlePrompt | AbstractPrompt | section_reference | meta_analysis_info
1            | 1                | "Screen..."  | "Evaluate..."  | "Topic"          | "For MA"
0            | 1                | ""           | "Review..."    | "No Reference"   | "Criteria"
```

## Output Data Formats

### 1. Selected PMIDs

**Files**: 
- `GPT_Selected_Titles.txt` (if title screening performed)
- `GPT_Selected_Abstracts.txt` (if abstract screening performed)

**Format**: Plain text, identical to input PMID format

**Content**: PMIDs that the LLM classified as relevant

### 2. PubMed Data Export

**File**: `PubMed_Data.xlsx`

**Format**: Excel workbook with retrieved paper metadata

**Columns**:
- `pmid`: PubMed ID
- `title`: Paper title
- `abstract`: Paper abstract text
- Additional metadata as available

### 3. Performance Results

**File**: Updated `Prompts.xlsx` with additional columns

**Added Columns**:

| Column Name | Data Type | Description |
|------------|-----------|-------------|
| `sensitivity_titles` | Float | Title screening sensitivity |
| `specificity_titles` | Float | Title screening specificity |
| `PPV_titles` | Float | Title screening positive predictive value |
| `NPV_titles` | Float | Title screening negative predictive value |
| `tp_titles` | Integer | True positives (titles) |
| `tn_titles` | Integer | True negatives (titles) |
| `fp_titles` | Integer | False positives (titles) |
| `fn_titles` | Integer | False negatives (titles) |
| `sensitivity_abstracts` | Float | Abstract screening sensitivity |
| `specificity_abstracts` | Float | Abstract screening specificity |
| `PPV_abstracts` | Float | Abstract screening positive predictive value |
| `NPV_abstracts` | Float | Abstract screening negative predictive value |
| `tp_abstracts` | Integer | True positives (abstracts) |
| `tn_abstracts` | Integer | True negatives (abstracts) |
| `fp_abstracts` | Integer | False positives (abstracts) |
| `fn_abstracts` | Integer | False negatives (abstracts) |
| `Output` | Text | Raw LLM responses (for debugging) |

## Database Schema

### SQLite Database (`ensure.sqlite`)

**Table: Initial**
```sql
CREATE TABLE Initial (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL,
    title TEXT,
    abstract TEXT
);
```

**Table: titles**
```sql
CREATE TABLE titles (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL,
    title TEXT,
    abstract TEXT
);
```

**Table: abstracts**
```sql
CREATE TABLE abstracts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL,
    title TEXT,
    abstract TEXT
);
```

## Prompt Engineering Data Formats

### 1. Initial Prompts

**File**: `initial_prompts.csv`

**Format**: CSV with semicolon delimiter

**Columns**:
- `Prompt`: Base prompt text
- `section_reference`: Reference type
- `meta_analysis_info`: Meta-analysis context
- `words`: Word count
- Additional categorical variables

### 2. Generated Prompts

**Files**: 
- `generated_prompts.csv`: Markov chain generated
- `generated_prompts_llm.csv`: LLM enhanced
- `complete_prompts.csv`: Combined dataset
- `final_prompts.csv`: Final optimized set (1,440 prompts)

**Format**: CSV with semicolon delimiter

**Structure**: Same as initial prompts with additional generated variations

## Statistical Analysis Data Formats

### 1. Analysis Input

**File**: `ENSURE_Data.xlsx`

**Format**: Excel workbook for R analysis

**Required Columns**:
- All prompt configuration columns
- All performance metric columns
- `f_score_2_1`: Calculated F-score with 2:1 sensitivity weighting

### 2. R Markdown Output

**Files**: 
- `ensure_analysis.html`: Rendered analysis report
- Various plot files (`.png`, `.pdf`)

## Data Validation Guidelines

### Input Validation

1. **PMID Format**:
   ```python
   def validate_pmid(pmid):
       return pmid.isdigit() and len(pmid) >= 7 and len(pmid) <= 9
   ```

2. **File Encoding**: Always use UTF-8
3. **Line Endings**: Prefer Unix (LF) format
4. **Excel Compatibility**: Save as .xlsx format, not .xls

### Output Validation

1. **Performance Metrics**: Should be between 0 and 1
2. **Count Consistency**: tp + tn + fp + fn should equal total papers
3. **PMID Consistency**: All output PMIDs should exist in input

## Error Handling

### Common Data Issues

1. **Missing PMIDs**: Log and continue processing
2. **Invalid Characters**: Clean or flag for manual review
3. **Encoding Issues**: Convert to UTF-8 automatically
4. **Empty Abstracts**: Handle gracefully in screening logic

### Data Quality Checks

```python
# Example validation function
def validate_data_quality(data):
    checks = {
        'missing_titles': data['title'].isna().sum(),
        'missing_abstracts': data['abstract'].isna().sum(),
        'duplicate_pmids': data['pmid'].duplicated().sum(),
        'invalid_pmids': (~data['pmid'].str.isdigit()).sum()
    }
    return checks
```

## File Organization Best Practices

### Directory Structure
```
project/
├── data/
│   ├── raw/           # Original input files
│   ├── processed/     # Cleaned data
│   └── outputs/       # Generated results
├── config/            # Configuration files
├── logs/              # Processing logs
└── backup/            # Data backups
```

### Naming Conventions
- Use descriptive, consistent names
- Include timestamps for versioned data
- Use underscores rather than spaces
- Include format in extension (`.txt`, `.xlsx`, `.csv`)

### Version Control
- Track data versions with git or similar
- Document changes in metadata
- Maintain backward compatibility when possible

## Performance Considerations

### Large Datasets
- Process in batches to manage memory
- Use database queries instead of loading all data
- Implement progress tracking for long operations

### Data Storage
- Compress large text files when appropriate
- Use efficient data types (int32 vs int64)
- Consider columnar formats (Parquet) for large datasets

This specification ensures consistent data handling across all components of the ENSURE pipeline and facilitates reproducible research.
