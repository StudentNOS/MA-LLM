# ENSURE: Systematic Review Automation with Large Language Models

## Overview

ENSURE (Efficient aNd Scalable sUrvey REview) is a comprehensive research project designed to automate systematic review processes using Large Language Models (LLMs). This repository contains the complete pipeline for automated literature screening, prompt engineering, and statistical analysis used in our medical research publication.

## Project Structure

```
ensure/
├── 0_docs/                     # Project documentation and guidelines
├── 1_goldstandards/           # Gold standard datasets for validation
│   ├── 1. 38957929/          # Individual meta-analysis datasets
│   ├── 2. 39053746/          # Each containing:
│   ├── ...                   #   - PDF papers (_ft.pdf)
│   └── 10. 39779764/         #   - Supplementary materials (_supp.*)
│                             #   - Word documents (_word.docx)
│                             #   - Gold standard selections (Goldstandard_Selected.txt)
│                             #   - Initial PMID lists (Initial.txt)
│                             #   - Prompt configurations (Prompts.xlsx)
├── 2_promptengineering/       # Prompt development and optimization
│   ├── prompt_engineering.ipynb  # Main prompt engineering notebook
│   ├── initial_prompts.csv       # Base prompts for analysis
│   ├── generated_prompts.csv     # Markov chain generated prompts
│   ├── generated_prompts_llm.csv # LLM-enhanced prompts
│   ├── complete_prompts.csv      # Combined prompt dataset
│   └── final_prompts.csv         # Final optimized prompts (1,440 total)
├── 3_pipeline/                # Core processing pipeline
│   ├── ENSURE.py             # Main pipeline orchestrator
│   ├── PubMed.py             # PubMed data retrieval
│   ├── GPT.py                # LLM screening interface
│   ├── Performance.py        # Performance metrics calculation
│   ├── dbconnect.py          # Database connectivity
│   ├── PATHS.py              # Configuration paths and API keys
│   ├── Scramble.py           # Data randomization utilities
│   ├── Clear.py              # Database cleanup utilities
│   └── ensure.sqlite         # Local SQLite database
├── 4_statisticalanalyses/     # Statistical analysis and reporting
│   ├── ensure_analysis.Rmd   # R Markdown analysis file
│   ├── ENSURE_Data.xlsx      # Processed data for analysis
│   └── Power Analysis.rmd    # Statistical power analysis
├── 5_placeholders/            # Template files and utilities
│   ├── Prompts final.xlsx    # Final prompt templates
│   └── Replace.py            # Text replacement utilities
├── requirements.txt           # Python dependencies
└── README.md                 # This documentation file
```

## Research Methodology

### 1. Gold Standard Creation (`1_goldstandards/`)

Each subdirectory represents a meta-analysis with:
- **Full-text PDFs**: Original research papers
- **Supplementary materials**: Additional data and appendices
- **Word documents**: Formatted manuscripts
- **Gold standard selections**: Manually curated relevant papers (PMIDs)
- **Initial datasets**: Complete PMID lists for screening
- **Prompt configurations**: Screening prompts and parameters

### 2. Prompt Engineering (`2_promptengineering/`)

Our prompt development process includes:

- **Initial Prompts**: Base prompts derived from systematic review guidelines
- **Markov Chain Enhancement**: Syntactic variation using statistical language models
- **LLM Augmentation**: Semantic enhancement using large language models
- **Validation**: Performance testing across multiple meta-analyses

The final dataset contains 1,440 unique prompts across different screening scenarios.

### 3. Automated Pipeline (`3_pipeline/`)

The core screening pipeline consists of:

1. **Data Retrieval** (`PubMed.py`): Fetch paper metadata from PubMed
2. **LLM Screening** (`GPT.py`): Automated relevance assessment
3. **Performance Evaluation** (`Performance.py`): Calculate screening metrics
4. **Database Management** (`dbconnect.py`): Data persistence and retrieval

### 4. Statistical Analysis (`4_statisticalanalyses/`)

Comprehensive statistical evaluation including:
- Mixed-effects models for prompt performance
- Sensitivity and specificity analysis
- Power analysis for sample size determination
- Visualization and reporting

## Installation and Setup

### Prerequisites

- Python 3.8+
- R 4.0+
- SQLite 3
- OpenAI API access (or compatible LLM endpoint)

### Python Dependencies

```bash
pip install -r requirements.txt
```

Required packages:
- `pandas`: Data manipulation and analysis
- `openai`: LLM API interface
- `matplotlib`: Data visualization
- `biopython`: PubMed data retrieval
- `openpyxl`: Excel file handling

### R Dependencies

```r
install.packages(c("dplyr", "tidyr", "readxl", "lme4", "simr", 
                   "mixedpower", "ggplot2", "sjPlot"))
```

### Configuration

1. **API Configuration** (`3_pipeline/PATHS.py`):
   ```python
   api_key = 'your_openai_api_key_here'
   email = 'your_email@domain.com'  # Required for PubMed access
   ```

2. **Database Setup**:
   The SQLite database will be automatically created on first run.

## Usage

### Quick Start

1. **Configure API credentials** in `3_pipeline/PATHS.py`
2. **Prepare your PMID list** as `Initial.txt` in the pipeline directory
3. **Run the main pipeline**:
   ```bash
   cd 3_pipeline
   python ENSURE.py
   ```

### Detailed Workflow

#### 1. Data Preparation

Place your initial PMID list in `3_pipeline/Initial.txt`:
```
32420649
22901346
34025474
...
```

#### 2. Configure Screening Parameters

Update `3_pipeline/Prompts.xlsx` with your screening configuration:
- `screen_titles`: Boolean flag for title screening
- `screen_abstracts`: Boolean flag for abstract screening
- `TitlePrompt`: Prompt text for title screening
- `AbstractPrompt`: Prompt text for abstract screening

#### 3. Execute Pipeline

The main script will:
1. Randomize the initial PMID list
2. Fetch paper details from PubMed
3. Perform LLM-based screening (titles and/or abstracts)
4. Calculate performance metrics against gold standard
5. Generate Excel reports and visualizations

#### 4. Analyze Results

Use the R Markdown files in `4_statisticalanalyses/` to:
- Perform statistical analysis
- Generate performance visualizations
- Calculate power analysis for sample size planning

### Performance Metrics

The pipeline calculates standard screening performance metrics:

- **Sensitivity (Recall)**: True positives / (True positives + False negatives)
- **Specificity**: True negatives / (True negatives + False positives)
- **Positive Predictive Value (Precision)**: True positives / (True positives + False positives)
- **Negative Predictive Value**: True negatives / (True negatives + False negatives)
- **F-Score**: Harmonic mean of precision and recall

## Reproducibility Guidelines

### For Exact Reproduction

1. **Environment Setup**:
   - Use the exact Python and R versions specified
   - Install dependencies from `requirements.txt`
   - Use the same LLM model and version

2. **Data Consistency**:
   - Use the same initial PMID lists from `1_goldstandards/`
   - Apply identical randomization seeds
   - Maintain consistent API endpoints

3. **Configuration Management**:
   - Use the provided prompt configurations
   - Maintain consistent batch sizes and parameters
   - Document any deviations from default settings

### For Adaptation to New Datasets

1. **Prepare Gold Standards**:
   - Create manual relevance assessments
   - Format as PMID lists in text files
   - Validate inter-rater agreement

2. **Customize Prompts**:
   - Adapt prompts to your research domain
   - Use the prompt engineering notebook for optimization
   - Test performance across multiple configurations

3. **Validate Performance**:
   - Compare against human reviewers
   - Calculate confidence intervals
   - Perform sensitivity analysis

## File Format Specifications

### PMID Lists
- **Format**: Plain text, one PMID per line
- **Example**: `Initial.txt`, `Goldstandard_Selected.txt`
- **Encoding**: UTF-8

### Prompt Configuration
- **Format**: Excel workbook (.xlsx)
- **Required columns**: `screen_titles`, `screen_abstracts`, `TitlePrompt`, `AbstractPrompt`
- **Boolean values**: 1 (true), 0 (false)

### Performance Data
- **Format**: Excel workbook (.xlsx)
- **Generated automatically** by the pipeline
- **Includes**: Confusion matrices, performance metrics, raw outputs

## Troubleshooting

### Common Issues

1. **API Rate Limits**:
   - Implement delays between requests
   - Use batch processing for large datasets
   - Monitor API usage quotas

2. **Memory Issues**:
   - Process data in smaller batches
   - Clear intermediate results
   - Monitor database size

3. **Missing Dependencies**:
   - Verify all packages are installed
   - Check Python/R version compatibility
   - Update pip/CRAN packages

### Debug Mode

Enable detailed logging by modifying the pipeline scripts:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Contributing

### Code Standards
- Follow PEP 8 for Python code
- Use meaningful variable names
- Include docstrings for functions
- Add type hints where appropriate

### Documentation
- Update README for significant changes
- Document new features thoroughly
- Include examples for complex functions

### Testing
- Test on multiple meta-analyses
- Validate against known benchmarks
- Report performance metrics

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@article{ensure2025,
  title={ENSURE: Efficient aNd Scalable sUrvey REview - Systematic Review Automation with Large Language Models},
  author={[Author Names]},
  journal={[Journal Name]},
  year={2025},
  doi={[DOI]}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions, issues, or collaboration inquiries:
- **Email**: [contact_email@domain.com]
- **GitHub Issues**: Use the issue tracker for bug reports and feature requests
- **Research Group**: [Institution/Lab Name]

## Acknowledgments

- PubMed/NCBI for providing access to biomedical literature
- OpenAI for Large Language Model access
- R and Python communities for statistical and data analysis tools
- Systematic review methodology guidelines (PRISMA, Cochrane)

---

**Note**: This repository contains the complete implementation for our systematic review automation research. For production use, please consider additional validation, error handling, and scalability optimizations appropriate for your specific use case.
