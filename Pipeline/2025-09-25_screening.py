import os
import random
import sqlite3
import pandas as pd
import json
import re
import time
import logging
import sys
import io
import threading
from typing import List, Dict, Tuple, Optional, Any
from Bio import Entrez
from groq import Groq
from openai import OpenAI
from anthropic import Anthropic
import google.generativeai as genai
from ollama import Client as OllamaClient
from flask import Flask, render_template, request, jsonify, send_file, Response
import threading
import webbrowser
import queue
import uuid



# Remove logging configuration - using print statements with timestamps instead

app = Flask(__name__)

# Configuration and Status
DATABASE = 'ensure.sqlite'
ENTREZ_EMAIL = ''
API_KEY = ''
SELECTED_MODEL = ''
CURRENT_STATUS = 'idle'
RESULTS_DF = pd.DataFrame()
PROCESS_THREAD = None
STOP_REQUESTED = False

# Terminal output capture system
TERMINAL_OUTPUT_QUEUE = queue.Queue()
TERMINAL_CAPTURE_ACTIVE = False

# Global variables for output redirection
original_stdout = None
original_stderr = None
output_stream = None

def setup_terminal_capture():
    """
    Set up terminal output capture to redirect stdout and stderr.
    """
    global original_stdout, original_stderr, output_stream, TERMINAL_CAPTURE_ACTIVE

    if TERMINAL_CAPTURE_ACTIVE:
        return  # Already set up

    # Create a StringIO object to capture output
    output_stream = io.StringIO()

    # Save original stdout and stderr
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    # Redirect stdout and stderr to our StringIO object
    sys.stdout = output_stream
    sys.stderr = output_stream

    TERMINAL_CAPTURE_ACTIVE = True
    print("Terminal output capture started")

def restore_terminal_output():
    """
    Restore original stdout and stderr.
    """
    global original_stdout, original_stderr, output_stream, TERMINAL_CAPTURE_ACTIVE

    if not TERMINAL_CAPTURE_ACTIVE:
        return  # Already restored

    # Restore original stdout and stderr
    sys.stdout = original_stdout
    sys.stderr = original_stderr

    TERMINAL_CAPTURE_ACTIVE = False
    print("Terminal output capture stopped")

def get_captured_output():
    """
    Get the current captured terminal output.

    Returns:
        str: The captured output content
    """
    if output_stream:
        return output_stream.getvalue()
    return ""

def clear_captured_output():
    """
    Clear the captured output buffer.
    """
    if output_stream:
        output_stream.seek(0)
        output_stream.truncate(0)

def start_output_monitoring():
    """
    Start a background thread to monitor and stream terminal output.
    """
    def monitor_output():
        last_output_length = 0

        while TERMINAL_CAPTURE_ACTIVE:
            try:
                current_output = get_captured_output()

                if len(current_output) > last_output_length:
                    # New output was added
                    new_output = current_output[last_output_length:]

                    # Split into lines and send each line
                    lines = new_output.strip().split('\n')
                    for line in lines:
                        if line.strip():  # Only send non-empty lines
                            timestamp = time.strftime('%H:%M:%S')
                            formatted_message = {
                                'timestamp': timestamp,
                                'message': line.strip(),
                                'type': 'terminal'
                            }

                            try:
                                TERMINAL_OUTPUT_QUEUE.put_nowait(formatted_message)
                            except queue.Full:
                                # Remove oldest message if queue is full
                                try:
                                    TERMINAL_OUTPUT_QUEUE.get_nowait()
                                    TERMINAL_OUTPUT_QUEUE.put_nowait(formatted_message)
                                except queue.Empty:
                                    pass

                    last_output_length = len(current_output)

                time.sleep(0.1)  # Check every 100ms

            except Exception as e:
                print(f"Error in output monitoring: {str(e)}")
                break

    monitor_thread = threading.Thread(target=monitor_output, daemon=True)
    monitor_thread.start()
    print("Output monitoring started")

def get_debug_messages():
    """
    Generator function to stream debug messages to clients.
    Now includes both debug messages and terminal output.

    Yields:
        str: JSON formatted debug messages
    """
    client_id = str(uuid.uuid4())

    try:
        # Send initial connection message
        initial_message = {
            'timestamp': time.strftime('%H:%M:%S'),
            'message': f'Client {client_id[:8]} connected to debug stream',
            'type': 'info'
        }
        yield f"data: {json.dumps(initial_message)}\n\n"

        while True:
            # Check for terminal output
            try:
                terminal_message = TERMINAL_OUTPUT_QUEUE.get_nowait()
                yield f"data: {json.dumps(terminal_message)}\n\n"
                continue
            except queue.Empty:
                pass

            # No messages available, send heartbeat
            try:
                terminal_message = TERMINAL_OUTPUT_QUEUE.get(timeout=1.0)
                yield f"data: {json.dumps(terminal_message)}\n\n"
            except queue.Empty:
                # Send heartbeat
                heartbeat = {
                    'timestamp': time.strftime('%H:%M:%S'),
                    'message': 'heartbeat',
                    'type': 'heartbeat'
                }
                yield f"data: {json.dumps(heartbeat)}\n\n"

    except GeneratorExit:
        # Client disconnected
        pass
    except Exception as e:
        error_message = {
            'timestamp': time.strftime('%H:%M:%S'),
            'message': f'Stream error: {str(e)}',
            'type': 'error'
        }
        yield f"data: {json.dumps(error_message)}\n\n"

# Retry configuration
MAX_RETRIES = 3
RETRY_DELAY = 2  # seconds
BACKOFF_FACTOR = 2

class AIModelClient:
    """
    Client for interacting with various AI model providers.
    
    This class provides a unified interface for different AI providers including:
    - Groq
    - OpenAI 
    - Anthropic (Claude)
    - Google (Gemini)
    - Ollama (local models)
    
    Attributes:
        provider (str): The AI provider name (lowercase)
        model_name (str): The specific model to use
        api_key (str): API key for authentication
        client: The initialized client for the specific provider
    
    Raises:
        ValueError: If the provider is not supported
    """
    
    def __init__(self, provider: str, model_name: str, api_key: str):
        """
        Initialize the AI model client.
        
        Args:
            provider (str): AI provider name (e.g., 'groq', 'openai', 'anthropic', 'google', 'ollama')
            model_name (str): Specific model name to use
            api_key (str): API key for the provider
        """
        self.provider = provider.lower()
        self.model_name = model_name
        self.api_key = api_key
        self.supported_providers = ['groq', 'openai', 'anthropic', 'google', 'ollama']

        if self.provider not in self.supported_providers:
            supported_list = ', '.join(self.supported_providers)
            raise ValueError(f"Provider '{provider}' is not supported. Supported providers: {supported_list}. "
                           "New providers can be added by extending this class.")

        self.client = self._init_client()

    def _init_client(self):
        """
        Initialize the appropriate client based on the provider.
        
        Returns:
            The initialized client object for the specific provider
        """
        if self.provider == 'groq':
            return Groq(api_key=self.api_key)
        elif self.provider == 'openai':
            return OpenAI(api_key=self.api_key)
        elif self.provider == 'anthropic':
            return Anthropic(api_key=self.api_key)
        elif self.provider == 'google':
            genai.configure(api_key=self.api_key)
            return genai
        elif self.provider == 'ollama':
            return OllamaClient(host='http://localhost:11434')
        
        else:
            raise ValueError(f"Unknown provider: {self.provider}")

    def generate_completion(self, prompt: str, retries: int = MAX_RETRIES) -> str:
        """
        Generate a completion using the configured AI model with retry logic.
        
        This method attempts to generate a response from the AI model, with automatic
        retry logic for network errors, API rate limits, and token errors.
        
        Args:
            prompt (str): The prompt to send to the AI model
            retries (int): Number of retry attempts for failed requests
            
        Returns:
            str: The AI model's response text
            
        Raises:
            Exception: If all retry attempts fail or if there's an authentication error
        """
        last_exception = None
        
        for attempt in range(retries + 1):
            try:
                if self.provider == 'groq':
                    response = self.client.chat.completions.create(
                        model=self.model_name,
                        messages=[{"role": "user", "content": prompt}]
                    )
                    return response.choices[0].message.content

                elif self.provider == 'openai':
                    response = self.client.chat.completions.create(
                        model=self.model_name,
                        messages=[{"role": "user", "content": prompt}]
                    )
                    return response.choices[0].message.content

                elif self.provider == 'anthropic':
                    response = self.client.messages.create(
                        model=self.model_name,
                        messages=[{"role": "user", "content": prompt}]
                    )
                    return response.content[0].text

                elif self.provider == 'google':
                    model = self.client.GenerativeModel(self.model_name)
                    response = model.generate_content(prompt)
                    return response.text

                elif self.provider == 'ollama':
                    response = self.client.generate(model=self.model_name, prompt=prompt)
                    return response['response']

            except Exception as e:
                last_exception = e
                error_msg = f"Error communicating with {self.provider} {self.model_name}: {str(e)}"
                
                # Check if this is a retryable error
                is_retryable = self._is_retryable_error(e)
                
                if not is_retryable or attempt == retries:
                    if "API key" in str(e).lower() or "key" in str(e).lower():
                        error_msg += " - Please check your API key."
                    elif "model" in str(e).lower():
                        error_msg += " - The model might not be available."
                    print(f"ERROR: Final attempt failed: {error_msg}")
                    raise Exception(error_msg)

                # Wait before retrying with exponential backoff
                wait_time = RETRY_DELAY * (BACKOFF_FACTOR ** attempt)
                print(f"WARNING: Attempt {attempt + 1} failed, retrying in {wait_time}s: {str(e)}")
                time.sleep(wait_time)
        
        # This should never be reached, but just in case
        raise last_exception

    def _is_retryable_error(self, error: Exception) -> bool:
        """
        Determine if an error is retryable.
        
        Args:
            error (Exception): The exception to check
            
        Returns:
            bool: True if the error is retryable, False otherwise
        """
        error_str = str(error).lower()
        
        # Network-related errors that are retryable
        retryable_errors = [
            'connection error',
            'timeout',
            'rate limit',
            'too many requests',
            'server error',
            'internal server error',
            'service unavailable',
            'network error',
            'connection reset',
            'temporary failure'
        ]
        
        # Check for authentication errors (not retryable)
        auth_errors = [
            'api key',
            'authentication',
            'unauthorized',
            'forbidden',
            'invalid key'
        ]
        
        for auth_error in auth_errors:
            if auth_error in error_str:
                return False
                
        for retryable_error in retryable_errors:
            if retryable_error in error_str:
                return True
                
        # Default to retryable for unknown errors
        return True

def init_db():
    """
    Initialize the SQLite database for storing articles.
    
    Creates the Articles table with the following schema:
    - id: Primary key (auto-increment)
    - pmid: PubMed ID (text)
    - title: Article title (text)
    - abstract: Article abstract (text)
    - stage: Current processing stage (text, default: 'initial')
    - relevant: Relevance flag (integer, default: 0)
    """
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("DROP TABLE IF EXISTS Articles")
        cursor.execute("""
            CREATE TABLE Articles (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT,
                title TEXT,
                abstract TEXT,
                stage TEXT DEFAULT 'initial',
                relevant INTEGER DEFAULT 0
            );
        """)
        conn.commit()

def fetch_articles_by_pmids(pmids: List[str], retries: int = MAX_RETRIES) -> List[Tuple]:
    """
    Fetch article details from PubMed by PMIDs with retry logic.
    
    Args:
        pmids (List[str]): List of PubMed IDs to fetch
        retries (int): Number of retry attempts for failed requests
        
    Returns:
        List[Tuple]: List of tuples containing (pmid, title, abstract)
        
    Raises:
        Exception: If unable to fetch articles after all retries
    """
    last_exception = None
    
    for attempt in range(retries + 1):
        try:
            Entrez.email = ENTREZ_EMAIL
            handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            articles = []
            for article in records.get("PubmedArticle", []):
                pmid = article["MedlineCitation"]["PMID"]
                title = article["MedlineCitation"]["Article"].get("ArticleTitle", "")
                
                abstract = ""
                if "Abstract" in article["MedlineCitation"]["Article"]:
                    abstract = " ".join(
                        str(t) for t in article["MedlineCitation"]["Article"]["Abstract"].get("AbstractText", [])
                    )
                
                articles.append((pmid, title, abstract))
            
            return articles
            
        except Exception as e:
            last_exception = e
            error_str = str(e).lower()
            
            # Check if this is a retryable error
            is_retryable = any(term in error_str for term in ['timeout', 'connection', 'server error'])
            
            if not is_retryable or attempt == retries:
                print(f"ERROR: Failed to fetch articles after {retries + 1} attempts: {str(e)}")
                raise Exception(f"Failed to fetch articles from PubMed: {str(e)}")

            # Wait before retrying
            wait_time = RETRY_DELAY * (BACKOFF_FACTOR ** attempt)
            print(f"WARNING: Attempt {attempt + 1} failed, retrying in {wait_time}s: {str(e)}")
            time.sleep(wait_time)
    
    # This should never be reached
    raise last_exception

def search_pubmed(query: str, retmax: int = 100000, retries: int = MAX_RETRIES) -> List[str]:
    """
    Search PubMed and return list of PMIDs with retry logic.
    
    Args:
        query (str): PubMed search query
        retmax (int): Maximum number of results to return
        retries (int): Number of retry attempts for failed requests
        
    Returns:
        List[str]: List of PubMed IDs matching the query
        
    Raises:
        Exception: If unable to search PubMed after all retries
    """
    last_exception = None
    
    for attempt in range(retries + 1):
        try:
            Entrez.email = ENTREZ_EMAIL
            handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
            results = Entrez.read(handle)
            handle.close()
            return results["IdList"]
            
        except Exception as e:
            last_exception = e
            error_str = str(e).lower()
            
            # Check if this is a retryable error
            is_retryable = any(term in error_str for term in ['timeout', 'connection', 'server error'])
            
            if not is_retryable or attempt == retries:
                print(f"ERROR: Failed to search PubMed after {retries + 1} attempts: {str(e)}")
                raise Exception(f"Failed to search PubMed: {str(e)}")

            # Wait before retrying
            wait_time = RETRY_DELAY * (BACKOFF_FACTOR ** attempt)
            print(f"WARNING: Attempt {attempt + 1} failed, retrying in {wait_time}s: {str(e)}")
            time.sleep(wait_time)
    
    # This should never be reached
    raise last_exception

def save_to_db(articles: List[Tuple]):
    """
    Save articles to the database.
    
    Args:
        articles (List[Tuple]): List of tuples containing (pmid, title, abstract)
    """
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.executemany(
            "INSERT INTO Articles (pmid, title, abstract) VALUES (?, ?, ?)",
            articles
        )
        conn.commit()

def screen_articles(stage: str, prompt: str, ai_client: AIModelClient, screen_level: str = 'both',
                   batch_size: int = 5, current_prompt: int = None, total_prompts: int = None) -> List[str]:
    """
    Screen articles using AI model with retry logic.
    
    This function processes articles in batches, sending them to the AI model
    for relevance screening based on the provided criteria.
    
    Args:
        stage (str): Current processing stage ('initial' or 'titles')
        prompt (str): Screening criteria prompt for the AI model
        ai_client (AIModelClient): Configured AI client
        screen_level (str): Level of screening ('titles', 'abstracts', or 'both')
        batch_size (int): Number of articles to process in each batch
        
    Returns:
        List[str]: List of relevant PMIDs found during screening
        
    Note:
        Articles are processed in batches to manage API rate limits and memory usage.
        The function updates the database with relevance scores and advances article stages.
    """
    global CURRENT_STATUS, STOP_REQUESTED

    results = []
    print(f"INFO: Starting screen_articles with stage={stage}, screen_level={screen_level}")

    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT pmid, title, abstract FROM Articles WHERE stage = ?",
            (stage,)
        )
        articles = cursor.fetchall()
        total = len(articles)
        print(f"INFO: Found {total} articles at stage {stage}")

        for i in range(0, total, batch_size):
            if STOP_REQUESTED:
                CURRENT_STATUS = "Screening stopped by user"
                return results

            # Create base status message with article progress
            if stage == 'initial' or screen_level == 'titles':
                base_status = f"Screening titles ({min(i+batch_size, total)}/{total})"
            else:
                base_status = f"Screening abstracts ({min(i+batch_size, total)}/{total})"

            # Add prompt progress if available
            if current_prompt is not None and total_prompts is not None:
                CURRENT_STATUS = f"{base_status} von Prompt ({current_prompt + 1}/{total_prompts})"
            else:
                CURRENT_STATUS = base_status
                
            batch = articles[i:i+batch_size]
            print(f"INFO: Processing batch {i//batch_size + 1} with {len(batch)} articles")

            if stage == 'initial' or screen_level == 'titles':
                formatted_articles = [f"PMID: {pmid}\nTitle: {title}" for pmid, title, _ in batch]
                content_type = "titles"
            else:
                formatted_articles = [f"PMID: {pmid}\nTitle: {title}\nAbstract: {abstract}" for pmid, title, abstract in batch]
                content_type = "abstracts"

            formatted = "\n\n---\n\n".join(formatted_articles)

            full_prompt = (
                f"Your task is to screen articles based on a criterion.\n\n"
                f"Criterion: {prompt}\n\n"
                f"Here is a list of articles with their PMIDs and {content_type}. Identify the relevant ones.\n\n"
                f"Articles:\n{formatted}\n\n"
                f"Respond with ONLY the PMIDs of the relevant articles, separated by commas. "
                f"Do not include any other text, explanations, or formatting. "
                f"For example, your response must be exactly: 12345678,23456789,34567890"
            )

            try:
                print("INFO: Sending batch to AI model for screening")
                response_text = ai_client.generate_completion(full_prompt).strip()

                print(f"INFO: AI response received: '{response_text}'")

                # Extract all potential PMIDs from the response using regex (7-9 digits)
                extracted_pmids = re.findall(r'\b\d{7,9}\b', response_text)
                print(f"INFO: Extracted PMIDs from response: {extracted_pmids}")

                # Get the list of PMIDs in the current batch
                batch_pmids = [pmid for pmid, _, _ in batch]
                print(f"INFO: Batch PMIDs: {batch_pmids}")

                # Only consider PMIDs that are actually in the batch
                relevant_pmids = [pmid for pmid in extracted_pmids if pmid in batch_pmids]
                print(f"INFO: Relevant PMIDs (filtered): {relevant_pmids}")

                # Update database with results
                for pmid, _, _ in batch:
                    relevant = 1 if pmid in relevant_pmids else 0
                    next_stage = "titles" if stage == "initial" else "abstracts"
                    cursor.execute(
                        "UPDATE Articles SET relevant = ?, stage = ? WHERE pmid = ?",
                        (relevant, next_stage, pmid)
                    )
                    if relevant:
                        results.append(pmid)

                print(f"INFO: Batch results: {len([r for r in results if r in batch_pmids])} relevant articles")

                # Commit after successful processing
                conn.commit()
                print("INFO: Database committed successfully")

            except Exception as e:
                print(f"ERROR: Error during screening: {str(e)}")
                CURRENT_STATUS = f"Error during screening: {str(e)}"

                # Even on error, we need to advance the stage for these articles
                # so they don't get stuck in the current stage
                print("WARNING: Advancing stage for batch articles despite error...")
                for pmid, _, _ in batch:
                    next_stage = "titles" if stage == "initial" else "abstracts"
                    cursor.execute(
                        "UPDATE Articles SET relevant = 0, stage = ? WHERE pmid = ?",
                        (next_stage, pmid)
                    )
                conn.commit()
                print("INFO: Database committed after error handling")
                continue

            time.sleep(1)

        conn.commit()

    print(f"INFO: Total results for stage {stage}: {len(results)}")
    return results

def calculate_metrics(goldstandard: List[str], stage: str) -> Dict[str, float]:
    """
    Calculate screening performance metrics.
    
    Compares the AI screening results against a gold standard to calculate:
    - Sensitivity (True Positive Rate)
    - Specificity (True Negative Rate) 
    - Positive Predictive Value (PPV)
    - Confusion matrix values (TP, FP, TN, FN)
    
    Args:
        goldstandard (List[str]): List of PMIDs known to be relevant
        stage (str): Processing stage to evaluate ('titles' or 'abstracts')
        
    Returns:
        Dict[str, float]: Dictionary containing all calculated metrics
    """
    print(f"INFO: Calculating metrics for stage '{stage}'")
    print(f"INFO: Goldstandard PMIDs: {goldstandard}")

    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT pmid, relevant FROM Articles WHERE stage = ?",
            (stage,)
        )
        data = cursor.fetchall()

    print(f"INFO: Retrieved {len(data)} articles from database at stage '{stage}'")

    if not data:
        print("WARNING: No data found in database for this stage!")
        return {
            'sensitivity': 0, 'specificity': 0, 'ppv': 0,
            'tp': 0, 'fp': 0, 'tn': 0, 'fn': 0
        }

    df = pd.DataFrame(data, columns=["pmid", "relevant"])
    print(f"INFO: DataFrame shape: {df.shape}")
    print(f"INFO: Relevant counts: {df['relevant'].value_counts().to_dict()}")

    # Debug: Show which PMIDs are in database vs goldstandard
    all_db_pmids = set(df['pmid'])
    gs_set = set(goldstandard)
    print(f"INFO: PMIDs in database: {sorted(all_db_pmids)}")
    print(f"INFO: PMIDs in goldstandard: {sorted(gs_set)}")
    print(f"INFO: Overlap (should be {len(gs_set)}): {len(all_db_pmids.intersection(gs_set))}")
    print(f"INFO: Missing from database: {gs_set.difference(all_db_pmids)}")
    print(f"INFO: Extra in database: {all_db_pmids.difference(gs_set)}")

    screened_in = set(df[df['relevant'] == 1]['pmid'])
    screened_out = set(df[df['relevant'] == 0]['pmid'])
    print(f"INFO: Screened in (relevant=1): {sorted(screened_in)}")
    print(f"INFO: Screened out (relevant=0): {sorted(screened_out)[:10]}...")  # First 10 only

    tp = len(gs_set.intersection(screened_in))
    fp = len(screened_in.difference(gs_set))
    fn = len(gs_set.intersection(screened_out))
    all_pmids = set(df['pmid'])
    tn = len(all_pmids.difference(gs_set).difference(screened_in))

    print(f"INFO: TP={tp}, FP={fp}, TN={tn}, FN={fn}")
    print(f"INFO: Goldstandard articles that were missed: {sorted(gs_set.intersection(screened_out))}")

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0

    print("INFO: Final metrics - Sensitivity: {:.3f}, Specificity: {:.3f}, PPV: {:.3f}".format(sensitivity, specificity, ppv))

    return {
        'sensitivity': sensitivity,
        'specificity': specificity,
        'ppv': ppv,
        'tp': tp,
        'fp': fp,
        'tn': tn,
        'fn': fn
    }

def process_prompts(prompts_df: pd.DataFrame, initial_pmids: List[str], goldstandard_pmids: List[str], 
                   ai_client: AIModelClient):
    """
    Process multiple prompts for gold standard comparison.
    
    This function processes a set of screening prompts against a collection of articles,
    comparing results against a gold standard to evaluate performance.
    
    Args:
        prompts_df (pd.DataFrame): DataFrame containing prompts and screening configuration
        initial_pmids (List[str]): List of PMIDs to screen
        goldstandard_pmids (List[str]): List of PMIDs known to be relevant
        ai_client (AIModelClient): Configured AI client for screening
        
    The prompts_df should contain columns:
        - TitlePrompt: Prompt for title screening
        - AbstractPrompt: Prompt for abstract screening  
        - screen_titles: 1 to screen titles, 0 to skip
        - screen_abstracts: 1 to screen abstracts, 0 to skip
    """
    global CURRENT_STATUS, RESULTS_DF, STOP_REQUESTED

    RESULTS_DF = prompts_df.copy()

    # Use all initial_pmids for screening, goldstandard is only for metrics
    goldstandard_set = set(goldstandard_pmids)
    filtered_initial_pmids = [pmid for pmid in initial_pmids if pmid in goldstandard_set]

    print(f"INFO: Original initial_pmids count: {len(initial_pmids)}")
    print(f"INFO: Goldstandard_pmids count: {len(goldstandard_pmids)}")
    print(f"INFO: Filtered initial_pmids count: {len(filtered_initial_pmids)}")
    print(f"INFO: Filtered PMIDs: {filtered_initial_pmids}")

    # Fetch all articles from Initial.txt, regardless of goldstandard
    articles = fetch_articles_by_pmids(initial_pmids)

    metric_names = ['sensitivity', 'specificity', 'ppv', 'tp', 'fp', 'tn', 'fn']
    output_columns = ['Final_Relevant_PMIDs']
    for prefix in ['title', 'abstract']:
        for metric in metric_names:
            output_columns.append(f'{prefix}_{metric}')

    for col in output_columns:
        RESULTS_DF[col] = None

    empty_metrics = {metric: None for metric in metric_names}

    for idx, row in RESULTS_DF.iterrows():
        if STOP_REQUESTED:
            CURRENT_STATUS = "Processing stopped by user"
            break

        CURRENT_STATUS = f"Processing prompt set {idx+1}/{len(RESULTS_DF)}"

        init_db()
        save_to_db(articles)

        print(f"INFO: Processing prompt {idx+1}")
        print(f"INFO: screen_titles value: {row.get('screen_titles')}, type: {type(row.get('screen_titles'))}")
        print(f"INFO: screen_abstracts value: {row.get('screen_abstracts')}, type: {type(row.get('screen_abstracts'))}")
        print(f"INFO: TitlePrompt: {row.get('TitlePrompt')}")
        print(f"INFO: AbstractPrompt: {row.get('AbstractPrompt')}")

        # Convert to int for reliable comparison
        screen_titles_val = int(row['screen_titles']) if pd.notna(row.get('screen_titles')) else 0
        screen_abstracts_val = int(row['screen_abstracts']) if pd.notna(row.get('screen_abstracts')) else 0

        if screen_titles_val == 1:
            print("INFO: Screening titles...")
            screen_articles(stage='initial', prompt=row['TitlePrompt'], ai_client=ai_client,
                          current_prompt=idx, total_prompts=len(RESULTS_DF))
            title_metrics = calculate_metrics(goldstandard_pmids, 'titles')
            print(f"INFO: Title metrics: {title_metrics}")
        else:
            print("INFO: Skipping title screening")
            with sqlite3.connect(DATABASE) as conn:
                conn.cursor().execute("UPDATE Articles SET stage = 'titles' WHERE stage = 'initial'")
                conn.commit()
            title_metrics = empty_metrics.copy()

        for metric, value in title_metrics.items():
            RESULTS_DF.at[idx, f'title_{metric}'] = value

        if STOP_REQUESTED: break

        if screen_abstracts_val == 1:
            print("INFO: Screening abstracts...")
            screen_articles(stage='titles', prompt=row['AbstractPrompt'], ai_client=ai_client,
                          current_prompt=idx, total_prompts=len(RESULTS_DF))
            abstract_metrics = calculate_metrics(goldstandard_pmids, 'abstracts')
            print(f"INFO: Abstract metrics: {abstract_metrics}")
        else:
            print("INFO: Skipping abstract screening")
            abstract_metrics = empty_metrics.copy()

        for metric, value in abstract_metrics.items():
            RESULTS_DF.at[idx, f'abstract_{metric}'] = value

        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT pmid FROM Articles WHERE relevant = 1")
            final_pmids_list = [item[0] for item in cursor.fetchall()]
            RESULTS_DF.at[idx, 'Final_Relevant_PMIDs'] = ",".join(final_pmids_list)
            print(f"INFO: Final relevant PMIDs for prompt {idx+1}: {final_pmids_list}")

        if STOP_REQUESTED: break

    print(f"INFO: Final RESULTS_DF shape: {RESULTS_DF.shape}")
    print(f"INFO: Final RESULTS_DF columns: {list(RESULTS_DF.columns)}")
    print("INFO: Final RESULTS_DF content:")
    print(RESULTS_DF.to_string())

    output_path = "prompt_results.xlsx"
    RESULTS_DF.to_excel(output_path, index=False)
    print(f"INFO: Results saved to {output_path}")

    if STOP_REQUESTED:
        CURRENT_STATUS = f"Stopped! Partial results saved to {output_path}"
    else:
        CURRENT_STATUS = f"Completed! Results saved to {output_path}"

def process_freeform_search(pubmed_query: str, screening_prompt: str, screen_level: str, 
                           ai_client: AIModelClient, max_articles: int = 100000):
    """
    Process a freeform PubMed search and screening.
    
    This function performs a PubMed search, fetches articles, and screens them
    using the provided criteria without a gold standard comparison.
    
    Args:
        pubmed_query (str): PubMed search query
        screening_prompt (str): Screening criteria for the AI model
        screen_level (str): Level of screening ('titles', 'abstracts', or 'sequential')
        ai_client (AIModelClient): Configured AI client
        max_articles (int): Maximum number of articles to retrieve
    """
    global CURRENT_STATUS, RESULTS_DF, STOP_REQUESTED

    CURRENT_STATUS = "Searching PubMed..."
    pmids = search_pubmed(pubmed_query, retmax=max_articles)

    if not pmids:
        CURRENT_STATUS = "No articles found matching the PubMed query"
        return

    # Initialize RESULTS_DF early for intermediate exports
    RESULTS_DF = pd.DataFrame({
        'PubMed_Query': [pubmed_query],
        'Screening_Prompt': [screening_prompt],
        'Screen_Level': [screen_level],
        'Total_Articles': [len(pmids)],
        'Relevant_Articles': [0],
        'Relevant_PMIDs': ['']
    })

    CURRENT_STATUS = "Fetching article details..."
    articles = fetch_articles_by_pmids(pmids)

    init_db()
    save_to_db(articles)

    if screen_level == 'sequential':
        CURRENT_STATUS = "Screening titles..."
        screen_articles('initial', screening_prompt, ai_client, 'titles')

        # For sequential screening, only process abstracts of articles that were
        # marked as relevant in title screening (relevant = 1)
        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            # First, advance all articles to 'titles' stage
            cursor.execute("UPDATE Articles SET stage = 'titles' WHERE stage = 'initial'")
            # Then, only advance articles that were marked as relevant to 'abstracts' stage
            cursor.execute("UPDATE Articles SET stage = 'abstracts' WHERE stage = 'titles' AND relevant = 1")
            conn.commit()

        # Count how many articles made it to abstract screening
        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM Articles WHERE stage = 'abstracts'")
            abstracts_count = cursor.fetchone()[0]

        if abstracts_count > 0:
            CURRENT_STATUS = f"Screening abstracts of {abstracts_count} relevant titles..."
            screen_articles('abstracts', screening_prompt, ai_client, 'abstracts')

        # Get final relevant PMIDs (articles that passed both title and abstract screening)
        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT pmid FROM Articles WHERE relevant = 1")
            results = [item[0] for item in cursor.fetchall()]
    else:
        CURRENT_STATUS = "Screening articles..."
        results = screen_articles('initial', screening_prompt, ai_client, screen_level)

    # Update RESULTS_DF with final results
    RESULTS_DF.at[0, 'Relevant_Articles'] = len(results)
    RESULTS_DF.at[0, 'Relevant_PMIDs'] = ",".join(results)

    output_path = "screening_results.xlsx"
    RESULTS_DF.to_excel(output_path, index=False)

    if STOP_REQUESTED:
        CURRENT_STATUS = f"Stopped! Partial results saved to {output_path}"
    else:
        CURRENT_STATUS = f"Completed! Results saved to {output_path}"

@app.route('/')
def index():
    """
    Serve the main application interface.
    
    Returns:
        Rendered HTML template for the web interface
    """
    return render_template('2025-09-25_index.html')

def start_processing_thread(target_func, args_tuple):
    """
    Start a processing thread for long-running operations.

    Args:
        target_func: Function to run in the thread
        args_tuple: Tuple of arguments to pass to the function
    """
    global CURRENT_STATUS, PROCESS_THREAD, STOP_REQUESTED
    STOP_REQUESTED = False
    CURRENT_STATUS = "Starting screening process..."
    print("Starting background processing thread")
    PROCESS_THREAD = threading.Thread(target=target_func, args=args_tuple)
    PROCESS_THREAD.start()

@app.route('/run_comparison', methods=['POST'])
def run_comparison():
    """
    Handle gold standard comparison requests.
    
    Expects form data with:
    - entrez_email: Email for PubMed API
    - api_key: API key for AI provider
    - ai_provider: AI provider name
    - ai_model: AI model name
    - initial_file: File containing initial PMIDs
    - goldstandard_file: File containing gold standard PMIDs
    - prompts_file: Excel file with screening prompts
    
    Returns:
        JSON response indicating success or failure
    """
    global ENTREZ_EMAIL, API_KEY, SELECTED_MODEL
    ENTREZ_EMAIL = request.form['entrez_email']
    API_KEY = request.form['api_key']
    ai_provider = request.form['ai_provider']
    ai_model = request.form['ai_model']
    SELECTED_MODEL = f"{ai_provider}-{ai_model}"

    initial_file = request.files['initial_file']
    goldstandard_file = request.files['goldstandard_file']
    prompts_file = request.files['prompts_file']

    # Read and debug initial file
    initial_content = initial_file.read().decode('utf-8')
    initial_pmids = [line.strip() for line in initial_content.splitlines() if line.strip()]
    print("INFO: === DEBUG: Initial.txt Content ===")
    print(f"INFO: Raw content length: {len(initial_content)}")
    print(f"INFO: First 500 chars: {initial_content[:500]}")
    print(f"INFO: Total PMIDs found: {len(initial_pmids)}")
    print(f"INFO: First 10 PMIDs: {initial_pmids[:10]}")
    print(f"INFO: Last 10 PMIDs: {initial_pmids[-10:] if len(initial_pmids) > 10 else initial_pmids}")
    print("INFO: ================================")

    # Read and debug goldstandard file
    goldstandard_content = goldstandard_file.read().decode('utf-8')
    goldstandard_pmids = [line.strip() for line in goldstandard_content.splitlines() if line.strip()]
    print("INFO: === DEBUG: Goldstandard_Selected.txt Content ===")
    print(f"INFO: Raw content length: {len(goldstandard_content)}")
    print(f"INFO: Full content: {goldstandard_content}")
    print(f"INFO: Total PMIDs found: {len(goldstandard_pmids)}")
    print(f"INFO: All PMIDs: {goldstandard_pmids}")
    print("INFO: ================================================")

    # Read and debug prompts file
    prompts_df = pd.read_excel(prompts_file)
    print("INFO: === DEBUG: Prompts.xlsx Content ===")
    print(f"INFO: Shape: {prompts_df.shape}")
    print(f"INFO: Columns: {list(prompts_df.columns)}")
    print(f"INFO: Data types:\n{prompts_df.dtypes}")
    print(f"INFO: screen_titles unique values: {prompts_df['screen_titles'].unique()}")
    print(f"INFO: screen_abstracts unique values: {prompts_df['screen_abstracts'].unique()}")
    print("INFO: First 3 rows:")
    print(prompts_df.head(3).to_string())
    print("INFO: ================================")

    ai_client = AIModelClient(ai_provider, ai_model, API_KEY)

    start_processing_thread(process_prompts, (prompts_df, initial_pmids, goldstandard_pmids, ai_client))
    return jsonify({'status': 'started'})

@app.route('/run_freeform', methods=['POST'])
def run_freeform():
    """
    Handle freeform search requests.
    
    Expects form data with:
    - entrez_email: Email for PubMed API
    - api_key: API key for AI provider
    - ai_provider: AI provider name
    - ai_model: AI model name
    - pubmed_search: PubMed search query
    - screening_prompt: Screening criteria
    - screen_level: Level of screening
    - max_articles: Maximum articles to retrieve (optional)
    
    Returns:
        JSON response indicating success or failure
    """
    global ENTREZ_EMAIL, API_KEY, SELECTED_MODEL
    ENTREZ_EMAIL = request.form['entrez_email']
    API_KEY = request.form['api_key']
    ai_provider = request.form['ai_provider']
    ai_model = request.form['ai_model']
    SELECTED_MODEL = f"{ai_provider}-{ai_model}"

    pubmed_search = request.form['pubmed_search']
    screening_prompt = request.form['screening_prompt']
    screen_level = request.form['screen_level']
    max_articles = int(request.form.get('max_articles', 100000))

    ai_client = AIModelClient(ai_provider, ai_model, API_KEY)

    start_processing_thread(process_freeform_search, (pubmed_search, screening_prompt, screen_level, ai_client, max_articles))
    return jsonify({'status': 'started'})

@app.route('/status')
def get_status():
    """
    Get current processing status.
    
    Returns:
        JSON response with current status message
    """
    return jsonify({'status': CURRENT_STATUS})

@app.route('/stop', methods=['POST'])
def stop_processing():
    """
    Stop the current processing operation.
    
    Returns:
        JSON response confirming the stop request
    """
    global STOP_REQUESTED, PROCESS_THREAD
    STOP_REQUESTED = True
    if PROCESS_THREAD and PROCESS_THREAD.is_alive():
        PROCESS_THREAD.join(timeout=5.0)
    return jsonify({'status': 'stop_requested', 'message': CURRENT_STATUS})

@app.route('/results')
def get_results():
    """
    Download the results file.
    
    Returns:
        Excel file with screening results or error message
    """
    global RESULTS_DF
    if not RESULTS_DF.empty:
        # Determine the correct output file based on RESULTS_DF structure
        if 'Relevant_Articles' in RESULTS_DF.columns:
            # Freeform mode
            output_path = "screening_results.xlsx"
        else:
            # Goldstandard mode
            output_path = "prompt_results.xlsx"

        if os.path.exists(output_path):
            return send_file(output_path, as_attachment=True)
    return jsonify({'error': 'No results available'}), 404

@app.route('/export_intermediate')
def export_intermediate():
    """
    Export intermediate results during processing.

    Returns:
        Excel file with current progress or fallback data
    """
    global RESULTS_DF

    if not RESULTS_DF.empty:
        # Update RESULTS_DF with current database state for intermediate export
        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT pmid FROM Articles WHERE relevant = 1")
            pmids = [item[0] for item in cursor.fetchall()]

        # Update the relevant fields based on RESULTS_DF structure
        if 'Relevant_Articles' in RESULTS_DF.columns:
            # Freeform mode
            RESULTS_DF.at[0, 'Relevant_Articles'] = len(pmids)
            RESULTS_DF.at[0, 'Relevant_PMIDs'] = ",".join(pmids)
        elif 'Final_Relevant_PMIDs' in RESULTS_DF.columns:
            # Goldstandard mode - update the current row if processing
            # For now, just export as is, since goldstandard updates per prompt
            pass

        output_path = "intermediate_results.xlsx"
        RESULTS_DF.to_excel(output_path, index=False)
        return send_file(output_path, as_attachment=True)
    else:
        # Fallback: Export current database state if no RESULTS_DF available
        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT pmid FROM Articles WHERE relevant = 1")
            pmids = [item[0] for item in cursor.fetchall()]
        df = pd.DataFrame({'Relevant_PMIDs': pmids})
        output_path = "intermediate_results.xlsx"
        df.to_excel(output_path, index=False)
        return send_file(output_path, as_attachment=True)

@app.route('/debug_stream')
def debug_stream():
    """
    Stream debug messages to the frontend.

    Returns:
        Server-sent events stream with debug messages
    """
    return Response(get_debug_messages(), mimetype='text/event-stream')



@app.route('/get_terminal_output')
def get_terminal_output():
    """
    Get the current captured terminal output.

    Returns:
        JSON response with the captured terminal output
    """
    output = get_captured_output()
    return jsonify({'output': output})





if __name__ == '__main__':
    url = f"http://127.0.0.1:5000"

    if not os.path.exists('templates'):
        print("Creating 'templates' directory..., but that means you have no html file, which contains the UI, therefore you have to get it from Github.")
        os.makedirs('templates')

    # Set up terminal output capture
    setup_terminal_capture()
    start_output_monitoring()

    print("MA-LLM Screening Tool starting up")
    print(f"Server will be available at: {url}")
    print("Debug system initialized and ready")
    print("Terminal output capture active")

    # Check if we should run tests
    if len(os.sys.argv) > 1 and os.sys.argv[1] == 'test':
        print("No test function available")
    else:
        # Starte Flask-Server im Thread damit Browser-Öffnung nach Server-Start erfolgt
        import threading
        import time

        def open_browser():
            time.sleep(1.5)  # Warte bis Server gestartet ist
            webbrowser.open(url)
            print(f"'{url}' is being opened in your usual browser.")

        # Only open browser if not already running
        import threading
        browser_thread = threading.Thread(target=open_browser)
        browser_thread.daemon = True
        browser_thread.start()

        print("Starting Flask development server")
        app.run(debug=True, threaded=True, use_reloader=False)
