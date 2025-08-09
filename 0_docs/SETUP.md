# ENSURE Setup Guide

This guide provides step-by-step instructions for setting up the ENSURE pipeline for systematic review automation.

## Prerequisites

### Software Requirements
- Python 3.8 or higher
- R 4.0 or higher  
- Git (for cloning the repository)
- Text editor or IDE (recommended: VS Code, PyCharm, RStudio)

### Account Requirements
- OpenAI API account with credits (for GPT models)
- Email address (for PubMed API access)

## Installation Steps

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/ensure.git
cd ensure
```

### 2. Set Up Python Environment

#### Option A: Using virtual environment (recommended)
```bash
# Create virtual environment
python -m venv venv

# Activate virtual environment
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

#### Option B: Using conda
```bash
# Create conda environment
conda create -n ensure python=3.9
conda activate ensure

# Install dependencies
pip install -r requirements.txt
```

### 3. Configure API Access

1. **Copy the configuration template**:
   ```bash
   cd 3_pipeline
   cp config_template.py PATHS.py
   ```

2. **Edit PATHS.py** with your credentials:
   ```python
   # OpenAI API key (get from https://platform.openai.com/api-keys)
   api_key = 'sk-your-actual-api-key-here'
   
   # Your email for PubMed access
   email = 'your.email@institution.edu'
   ```

### 4. Set Up R Environment

Install required R packages:
```r
# Open R or RStudio and run:
install.packages(c(
    "dplyr", "tidyr", "readxl", "lme4", 
    "simr", "mixedpower", "ggplot2", "sjPlot",
    "data.table", "knitr", "rmarkdown"
))
```

### 5. Verify Installation

Test the Python environment:
```bash
cd 3_pipeline
python -c "import pandas, openai, matplotlib; print('Python dependencies OK')"
```

Test database connectivity:
```bash
python -c "from dbconnect import ENSURE; print('Database setup OK')"
```

## Quick Start Example

### 1. Prepare Your Data

Create a file `3_pipeline/Initial.txt` with PubMed IDs:
```
32420649
22901346
34025474
25687662
```

### 2. Configure Screening

Create or edit `3_pipeline/Prompts.xlsx` with columns:
- `screen_titles`: 1 (to screen titles) or 0 (to skip)
- `screen_abstracts`: 1 (to screen abstracts) or 0 (to skip)
- `TitlePrompt`: Your title screening prompt text
- `AbstractPrompt`: Your abstract screening prompt text

Example prompt:
```
Screen the following titles for relevance to systematic reviews about diabetes treatment. 
Select titles that discuss randomized controlled trials, meta-analyses, or systematic reviews 
related to diabetes management interventions.
```

### 3. Create Gold Standard (for validation)

Create `3_pipeline/Goldstandard_Selected.txt` with PMIDs that should be selected:
```
32420649
34025474
```

### 4. Run the Pipeline

```bash
cd 3_pipeline
python ENSURE.py
```

The pipeline will:
- Fetch paper details from PubMed
- Screen titles/abstracts using your prompts
- Calculate performance metrics
- Generate Excel reports

### 5. Analyze Results

Open `4_statisticalanalyses/ensure_analysis.Rmd` in RStudio and run the analysis.

## Advanced Configuration

### Custom LLM Endpoints

To use a different LLM endpoint, modify `GPT.py`:

```python
client = OpenAI(
    base_url="https://your-custom-endpoint/v1/",
    api_key="your-api-key",
)
```

### Batch Size Optimization

Adjust batch sizes in `GPT.py` based on your API limits:

```python
# For title screening
batch_size = 50  # Reduce if hitting rate limits

# For abstract screening  
batch_size = 6   # Abstracts are longer, use smaller batches
```

### Database Configuration

The pipeline uses SQLite by default. For larger datasets, consider PostgreSQL:

1. Install PostgreSQL
2. Modify `dbconnect.py` with PostgreSQL connection string
3. Update database schema as needed

## Troubleshooting

### Common Issues

1. **API Key Error**:
   - Verify your OpenAI API key is correct
   - Check you have sufficient credits
   - Ensure no extra spaces in PATHS.py

2. **PubMed Connection Issues**:
   - Verify your email is correctly set
   - Check internet connection
   - PubMed may have rate limits

3. **Memory Issues with Large Datasets**:
   - Reduce batch sizes
   - Process in smaller chunks
   - Clear intermediate results

4. **R Package Installation Issues**:
   - Update R to latest version
   - Install packages one by one to identify issues
   - Check CRAN mirror settings

### Getting Help

1. **Check the logs**: Look for error messages in the console output
2. **Verify configuration**: Double-check all settings in PATHS.py
3. **Test with small datasets**: Start with 10-20 PMIDs to verify setup
4. **Review the documentation**: Check README.md for detailed information

## Security Considerations

### API Key Protection
- Never commit API keys to version control
- Use environment variables for production deployments
- Rotate keys regularly

### Data Privacy
- Be aware of what data you send to external APIs
- Consider local LLM deployments for sensitive data
- Review your institution's data handling policies

## Performance Optimization

### For Large Datasets
- Use database indexing for faster queries
- Implement parallel processing for batch operations
- Consider caching results to avoid re-processing

### Cost Management
- Monitor API usage and costs
- Use smaller models for initial testing
- Implement request caching where appropriate

## Next Steps

1. **Validate with your data**: Test the pipeline with a small subset of your systematic review
2. **Optimize prompts**: Use the prompt engineering notebook to improve performance  
3. **Scale up**: Once validated, run on your full dataset
4. **Analyze results**: Use the statistical analysis tools to evaluate performance
5. **Publish**: Follow reproducibility guidelines for publication

## Support

For technical support:
- Review the troubleshooting section
- Check GitHub issues for similar problems
- Contact the research team via email

Remember to cite our work if you use this pipeline in your research!
