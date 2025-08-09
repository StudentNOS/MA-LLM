# Configuration file for ENSURE pipeline
# Copy this file to PATHS.py and fill in your actual credentials

# OpenAI API Configuration
# Get your API key from: https://platform.openai.com/api-keys
api_key = 'your_openai_api_key_here'

# Alternative: For HU Berlin LLM endpoint (if applicable)
# base_url = "https://llm1-compute.cms.hu-berlin.de/v1/"
# api_key = "required-but-not-used"  # Placeholder for HU LLM

# PubMed/Entrez Configuration
# Required for accessing PubMed API
# Use your institutional or personal email
email = 'your_email@institution.edu'

# Database Configuration
database_path = 'ensure.sqlite'

# Model Configuration
# Default model for OpenAI API
default_model = 'gpt-3.5-turbo'

# For HU Berlin or other endpoints, specify the model name
# default_model = 'ModelCloud/Mistral-Large-Instruct-2407-gptq-4bit'

# Batch Processing Configuration
title_batch_size = 50      # Number of titles to process per batch
abstract_batch_size = 6    # Number of abstracts to process per batch
max_tokens = 2048         # Maximum tokens for LLM response

# Rate Limiting (to avoid API limits)
request_delay = 1.0       # Delay between requests in seconds
batch_delay = 300         # Delay between large batches in seconds (5 minutes)

# File Paths (relative to pipeline directory)
initial_pmids_file = "Initial.txt"
goldstandard_file = "Goldstandard_Selected.txt"
prompts_config_file = "Prompts.xlsx"
output_titles_file = "GPT_Selected_Titles.txt"
output_abstracts_file = "GPT_Selected_Abstracts.txt"
excel_output_file = "PubMed_Data.xlsx"
