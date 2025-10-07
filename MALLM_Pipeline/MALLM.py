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


class TeeStream(io.TextIOBase):
    """
    A stream that writes to multiple destinations simultaneously.
    Acts like a 'tee' command - output goes to both terminal and capture.
    """

    __slots__ = ('streams', 'errors')

    def __init__(self, *streams):
        self.streams = streams
        self.errors = []  # Track any write errors

    def write(self, data):
        """Write data to all streams."""
        for stream in self.streams:
            try:
                stream.write(data)
                stream.flush()  # Ensure immediate output
            except Exception as e:
                # Track errors but don't stop processing
                self.errors.append(f"Stream error: {e}")

        return len(data)  # Return expected value for TextIOBase

    def flush(self):
        """Flush all streams."""
        for stream in self.streams:
            try:
                stream.flush()
            except Exception as e:
                self.errors.append(f"Flush error: {e}")



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
    Set up terminal output capture using tee approach.
    Maintains normal terminal output while also capturing for UI.
    """
    global original_stdout, original_stderr, output_stream, TERMINAL_CAPTURE_ACTIVE

    if TERMINAL_CAPTURE_ACTIVE:
        return  # Already set up

    # Create a StringIO object to capture output for UI
    output_stream = io.StringIO()

    # Save original stdout and stderr (for restoration later)
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    # Create tee streams that write to both terminal and capture
    sys.stdout = TeeStream(sys.stdout, output_stream)
    sys.stderr = TeeStream(sys.stderr, output_stream)

    TERMINAL_CAPTURE_ACTIVE = True
    print("Terminal output capture started (tee mode)")

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
            # Check for terminal output from the queue
            try:
                terminal_message = TERMINAL_OUTPUT_QUEUE.get(timeout=1.0)
                yield f"data: {json.dumps(terminal_message)}\n\n"
            except queue.Empty:
                # If no message, send a heartbeat to keep connection alive
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
    """
    
    def __init__(self, provider: str, model_name: str, api_key: str):
        self.provider = provider.lower()
        self.model_name = model_name
        self.api_key = api_key
        self.supported_providers = ['groq', 'openai', 'anthropic', 'google', 'ollama']
        if self.provider not in self.supported_providers:
            raise ValueError(f"Provider '{provider}' is not supported.")
        self.client = self._init_client()

    def _init_client(self):
        if self.provider == 'groq': return Groq(api_key=self.api_key)
        elif self.provider == 'openai': return OpenAI(api_key=self.api_key)
        elif self.provider == 'anthropic': return Anthropic(api_key=self.api_key)
        elif self.provider == 'google':
            genai.configure(api_key=self.api_key)
            return genai
        elif self.provider == 'ollama':
            try:
                client = OllamaClient(host='http://localhost:11434')
                # A quick test to see if the server is really responding
                client.list() 
                return client
            except Exception as e:
                print("ERROR: Failed to connect to Ollama client at http://localhost:11434.")
                print("ERROR: Please ensure Ollama is running and accessible.")
                # Create a meaningful error that stops the process.
                raise ConnectionError("Could not connect to Ollama. Please ensure it is running.")
        else: raise ValueError(f"Unknown provider: {self.provider}")

    def generate_completion(self, prompt: str, retries: int = MAX_RETRIES) -> str:
        last_exception = None
        for attempt in range(retries + 1):
            try:
                if self.provider == 'groq':
                    response = self.client.chat.completions.create(model=self.model_name, messages=[{"role": "user", "content": prompt}])
                    return response.choices[0].message.content
                elif self.provider == 'openai':
                    response = self.client.chat.completions.create(model=self.model_name, messages=[{"role": "user", "content": prompt}])
                    return response.choices[0].message.content
                elif self.provider == 'anthropic':
                    response = self.client.messages.create(model=self.model_name, max_tokens=2048, messages=[{"role": "user", "content": prompt}])
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
                if not self._is_retryable_error(e) or attempt == retries:
                    raise Exception(f"Error with {self.provider}: {e}")
                wait_time = RETRY_DELAY * (BACKOFF_FACTOR ** attempt)
                print(f"WARNING: Attempt {attempt + 1} failed, retrying in {wait_time}s: {e}")
                time.sleep(wait_time)
        raise last_exception

    def _is_retryable_error(self, error: Exception) -> bool:
        error_str = str(error).lower()
        retryable = ['connection error', 'timeout', 'rate limit', 'server error', 'service unavailable']
        non_retryable = ['api key', 'authentication', 'unauthorized', 'invalid key']
        if any(term in error_str for term in non_retryable): return False
        if any(term in error_str for term in retryable): return True
        return True

def init_db():
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("DROP TABLE IF EXISTS Articles")
        cursor.execute("""
            CREATE TABLE Articles (
                id INTEGER PRIMARY KEY AUTOINCREMENT, pmid TEXT, title TEXT,
                abstract TEXT, stage TEXT DEFAULT 'initial', relevant INTEGER DEFAULT 0);
        """)
        conn.commit()

def fetch_articles_by_pmids(pmids: List[str], retries: int = MAX_RETRIES) -> List[Tuple]:
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
                abstract = " ".join(str(t) for t in article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", []))
                articles.append((pmid, title, abstract))
            return articles
        except Exception as e:
            last_exception = e
            if 'http' not in str(e).lower() or attempt == retries:
                raise Exception(f"Failed to fetch from PubMed: {e}")
            wait_time = RETRY_DELAY * (BACKOFF_FACTOR ** attempt)
            print(f"WARNING: PubMed fetch failed, retrying in {wait_time}s: {e}")
            time.sleep(wait_time)
    raise last_exception

def search_pubmed(query: str, retmax: int = 100000, retries: int = MAX_RETRIES) -> List[str]:
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
            if 'http' not in str(e).lower() or attempt == retries:
                raise Exception(f"Failed to search PubMed: {e}")
            wait_time = RETRY_DELAY * (BACKOFF_FACTOR ** attempt)
            print(f"WARNING: PubMed search failed, retrying in {wait_time}s: {e}")
            time.sleep(wait_time)
    raise last_exception

def save_to_db(articles: List[Tuple]):
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.executemany("INSERT INTO Articles (pmid, title, abstract) VALUES (?, ?, ?)", articles)
        conn.commit()

def screen_articles(stage: str, prompt: str, ai_client: AIModelClient, screen_level: str,
                   batch_size: int = 5, current_prompt: int = None, total_prompts: int = None) -> List[str]:
    global CURRENT_STATUS, STOP_REQUESTED
    results = []
    print(f"INFO: Starting screen_articles with stage='{stage}', screen_level='{screen_level}'")
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT pmid, title, abstract FROM Articles WHERE stage = ?", (stage,))
        articles = cursor.fetchall()
        total = len(articles)
        if total == 0:
            print(f"WARNING: No articles to screen at stage {stage}.")
            return []

        for i in range(0, total, batch_size):
            if STOP_REQUESTED:
                CURRENT_STATUS = "Screening stopped by user"; return results

            batch = articles[i:i+batch_size]
            
            is_abstract_screening = (stage == 'initial' and screen_level == 'abstracts') or stage == 'abstracts'

            if is_abstract_screening:
                base_status = f"Screening abstracts ({min(i+batch_size, total)}/{total})"
                content_type = "abstracts"
                formatted_articles = [f"PMID: {pmid}\nTitle: {title}\nAbstract: {abstract}" for pmid, title, abstract in batch]
            else:  # Title screening
                base_status = f"Screening titles ({min(i+batch_size, total)}/{total})"
                content_type = "titles"
                formatted_articles = [f"PMID: {pmid}\nTitle: {title}" for pmid, title, _ in batch]

            if current_prompt is not None:
                CURRENT_STATUS = f"{base_status} for Prompt ({current_prompt + 1}/{total_prompts})"
            else:
                CURRENT_STATUS = base_status

            formatted = "\n\n---\n\n".join(formatted_articles)
            full_prompt = (f"Criterion: {prompt}\n\nArticles:\n{formatted}\n\n"
                           f"Respond ONLY with the PMIDs of relevant articles, separated by commas.")
            
            next_stage = "titles" if stage == "initial" else "abstracts"

            try:
                response_text = ai_client.generate_completion(full_prompt).strip()
                extracted_pmids = re.findall(r'\b\d{7,9}\b', response_text)
                batch_pmids = [pmid for pmid, _, _ in batch]
                relevant_pmids = [pmid for pmid in extracted_pmids if pmid in batch_pmids]

                for pmid, _, _ in batch:
                    relevant = 1 if pmid in relevant_pmids else 0
                    cursor.execute("UPDATE Articles SET relevant = ?, stage = ? WHERE pmid = ?", (relevant, next_stage, pmid))
                    if relevant: results.append(pmid)
                conn.commit()

            except Exception as e:
                print(f"ERROR: Error during screening batch: {e}")
                for pmid, _, _ in batch:
                    cursor.execute("UPDATE Articles SET relevant = 0, stage = ? WHERE pmid = ?", (next_stage, pmid))
                conn.commit()
                continue
            time.sleep(1)
    return results

def calculate_metrics(goldstandard: List[str], stage: str) -> Dict[str, float]:
    with sqlite3.connect(DATABASE) as conn:
        df = pd.read_sql_query("SELECT pmid, relevant, stage FROM Articles", conn)
    if df.empty: return {metric: "-" for metric in ['sensitivity', 'specificity', 'ppv', 'tp', 'fp', 'tn', 'fn']}
    
    gs_set = set(goldstandard)
    if stage == 'titles':
        all_pmids_in_scope = set(df['pmid'])
        screened_in = set(df[df['stage'] == 'abstracts']['pmid'])
    elif stage == 'abstracts':
        all_pmids_in_scope = set(df[df['stage'] == 'abstracts']['pmid'])
        screened_in = set(df[(df['stage'] == 'abstracts') & (df['relevant'] == 1)]['pmid'])
    else: return {metric: "-" for metric in ['sensitivity', 'specificity', 'ppv', 'tp', 'fp', 'tn', 'fn']}

    screened_out = all_pmids_in_scope - screened_in
    tp = len(gs_set.intersection(screened_in))
    fp = len(screened_in.difference(gs_set))
    fn = len(gs_set.intersection(screened_out))
    tn = len(all_pmids_in_scope.difference(gs_set).intersection(screened_out))

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0

    return {'sensitivity': sensitivity, 'specificity': specificity, 'ppv': ppv, 'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn}

def process_prompts(prompts_df: pd.DataFrame, initial_pmids: List[str], goldstandard_pmids: List[str], 
                   ai_client: AIModelClient):
    global CURRENT_STATUS, RESULTS_DF, STOP_REQUESTED
    RESULTS_DF = prompts_df.copy()
    articles = fetch_articles_by_pmids(initial_pmids)
    metric_names = ['sensitivity', 'specificity', 'ppv', 'tp', 'fp', 'tn', 'fn']
    for prefix in ['title', 'abstract']:
        for metric in metric_names: RESULTS_DF[f'{prefix}_{metric}'] = None
    RESULTS_DF['Final_Relevant_PMIDs'] = None

    for idx, row in RESULTS_DF.iterrows():
        if STOP_REQUESTED: break
        CURRENT_STATUS = f"Processing prompt set {idx+1}/{len(RESULTS_DF)}"
        init_db(); save_to_db(articles)
        
        screen_titles_val = int(row.get('screen_titles', 0))
        screen_abstracts_val = int(row.get('screen_abstracts', 0))
        
        if screen_titles_val == 1:
            screen_articles('initial', row['TitlePrompt'], ai_client, 'titles', current_prompt=idx, total_prompts=len(RESULTS_DF))
            with sqlite3.connect(DATABASE) as conn:
                conn.cursor().execute("UPDATE Articles SET stage = 'abstracts' WHERE stage = 'titles' AND relevant = 1")
            title_metrics = calculate_metrics(goldstandard_pmids, 'titles')
        else:
            with sqlite3.connect(DATABASE) as conn:
                conn.cursor().execute("UPDATE Articles SET stage = 'abstracts' WHERE stage = 'initial'")
            title_metrics = {metric: "-" for metric in metric_names}
        for metric, value in title_metrics.items(): RESULTS_DF.at[idx, f'title_{metric}'] = value
        if STOP_REQUESTED: break

        if screen_abstracts_val == 1:
            screen_articles('abstracts', row['AbstractPrompt'], ai_client, 'abstracts', current_prompt=idx, total_prompts=len(RESULTS_DF))
            abstract_metrics = calculate_metrics(goldstandard_pmids, 'abstracts')
        else:
            abstract_metrics = {metric: "-" for metric in metric_names}
            with sqlite3.connect(DATABASE) as conn:
                conn.cursor().execute("UPDATE Articles SET relevant = 1 WHERE stage = 'abstracts'")
        for metric, value in abstract_metrics.items(): RESULTS_DF.at[idx, f'abstract_{metric}'] = value

        with sqlite3.connect(DATABASE) as conn:
            pmids = pd.read_sql_query("SELECT pmid FROM Articles WHERE stage = 'abstracts' AND relevant = 1", conn)
            RESULTS_DF.at[idx, 'Final_Relevant_PMIDs'] = ",".join(pmids['pmid'].tolist()) if not pmids.empty else ""
    
    try:
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        filename = f"prompt_results_{timestamp}.xlsx"
        output_path = os.path.join(os.path.expanduser("~"), "Downloads", filename)
        
        print(f"INFO: Saving results to {output_path}")
        RESULTS_DF.to_excel(output_path, index=False)
        
        CURRENT_STATUS = "Stopped! Save partial results." if STOP_REQUESTED else f"Completed! Results saved to {output_path}"

    except (PermissionError, IOError) as e:
        print(f"FATAL ERROR: Could not save results to '{output_path}'. Reason: {e}")
        print("FATAL ERROR: Please check if you have write permissions for the Downloads folder.")
        CURRENT_STATUS = f"Error: Could not save file. Check permissions for Downloads folder."
    except Exception as e:
        print(f"FATAL ERROR: An unexpected error occurred while saving the results file: {e}")
        CURRENT_STATUS = "Error: An unexpected error occurred while saving results."

def process_freeform_search(pubmed_query: str, screening_prompt: str, screen_level: str, 
                           ai_client: AIModelClient, max_articles: int = 100000):
    global CURRENT_STATUS, RESULTS_DF, STOP_REQUESTED
    CURRENT_STATUS = "Searching PubMed..."
    pmids = search_pubmed(pubmed_query, retmax=max_articles)
    if not pmids: CURRENT_STATUS = "No articles found."; return

    RESULTS_DF = pd.DataFrame({'PubMed_Query': [pubmed_query], 'Screening_Prompt': [screening_prompt], 'Screen_Level': [screen_level],
                               'Total_Articles': [len(pmids)], 'Relevant_Articles': [0], 'Relevant_PMIDs': ['']})
    CURRENT_STATUS = "Fetching article details..."
    articles = fetch_articles_by_pmids(pmids)
    init_db(); save_to_db(articles)
    
    if screen_level == 'sequential':
        screen_articles('initial', screening_prompt, ai_client, 'titles')
        with sqlite3.connect(DATABASE) as conn:
            conn.cursor().execute("UPDATE Articles SET stage = 'abstracts' WHERE stage = 'titles' AND relevant = 1")
        screen_articles('abstracts', screening_prompt, ai_client, 'abstracts')
    else:
        screen_articles('initial', screening_prompt, ai_client, screen_level)
    
    with sqlite3.connect(DATABASE) as conn:
        final_pmids = pd.read_sql_query("SELECT pmid FROM Articles WHERE relevant = 1", conn)
    RESULTS_DF.at[0, 'Relevant_Articles'] = len(final_pmids)
    RESULTS_DF.at[0, 'Relevant_PMIDs'] = ",".join(final_pmids['pmid'].tolist())
    
    try:
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        filename = f"screening_results_{timestamp}.xlsx"
        output_path = os.path.join(os.path.expanduser("~"), "Downloads", filename)
        
        print(f"INFO: Saving results to {output_path}")
        RESULTS_DF.to_excel(output_path, index=False)
        
        CURRENT_STATUS = "Stopped! Save partial results." if STOP_REQUESTED else f"Completed! Results saved to {output_path}"

    except (PermissionError, IOError) as e:
        print(f"FATAL ERROR: Could not save results to '{output_path}'. Reason: {e}")
        print("FATAL ERROR: Please check if you have write permissions for the Downloads folder.")
        CURRENT_STATUS = f"Error: Could not save file. Check permissions for Downloads folder."
    except Exception as e:
        print(f"FATAL ERROR: An unexpected error occurred while saving the results file: {e}")
        CURRENT_STATUS = "Error: An unexpected error occurred while saving results."

@app.route('/')
def index(): return render_template('MALLM.html')

def start_processing_thread(target_func, args_tuple):
    global PROCESS_THREAD, STOP_REQUESTED
    STOP_REQUESTED = False
    PROCESS_THREAD = threading.Thread(target=target_func, args=args_tuple)
    PROCESS_THREAD.start()

@app.route('/run_comparison', methods=['POST'])
def run_comparison():
    global ENTREZ_EMAIL, API_KEY

    # --- STARTUP CHECK ---
    downloads_path = os.path.join(os.path.expanduser("~"), "Downloads")
    if not os.access(downloads_path, os.W_OK):
        CURRENT_STATUS = f"Error: No write permission for Downloads folder: {downloads_path}"
        print(f"VALIDATION ERROR: No write permission for Downloads folder: {downloads_path}")
        return jsonify({'status': 'error', 'message': CURRENT_STATUS}), 400 # Bad Request
    # --- END CHECK ---

    ENTREZ_EMAIL = request.form['entrez_email']
    ai_provider = request.form['ai_provider']

    # Handle API key for different providers
    if ai_provider.lower() == 'ollama':
        API_KEY = ''  # No API key needed for Ollama
    else:
        API_KEY = request.form['api_key']
        if not API_KEY.strip():
            CURRENT_STATUS = f"Error: API key is required for {ai_provider}"
            print(f"VALIDATION ERROR: API key is required for {ai_provider}")
            return jsonify({'status': 'error', 'message': CURRENT_STATUS}), 400

    ai_client = AIModelClient(ai_provider, request.form['ai_model'], API_KEY)
    initial_pmids = request.files['initial_file'].read().decode('utf-8').splitlines()
    goldstandard_pmids = request.files['goldstandard_file'].read().decode('utf-8').splitlines()
    prompts_df = pd.read_excel(request.files['prompts_file'])
    start_processing_thread(process_prompts, (prompts_df, initial_pmids, goldstandard_pmids, ai_client))
    return jsonify({'status': 'started'})

@app.route('/run_freeform', methods=['POST'])
def run_freeform():
    global ENTREZ_EMAIL, API_KEY

    # --- STARTUP CHECK ---
    downloads_path = os.path.join(os.path.expanduser("~"), "Downloads")
    if not os.access(downloads_path, os.W_OK):
        CURRENT_STATUS = f"Error: No write permission for Downloads folder: {downloads_path}"
        print(f"VALIDATION ERROR: No write permission for Downloads folder: {downloads_path}")
        return jsonify({'status': 'error', 'message': CURRENT_STATUS}), 400 # Bad Request
    # --- END CHECK ---

    ENTREZ_EMAIL = request.form['entrez_email']
    ai_provider = request.form['ai_provider']

    # Handle API key for different providers
    if ai_provider.lower() == 'ollama':
        API_KEY = ''  # No API key needed for Ollama
    else:
        API_KEY = request.form['api_key']
        if not API_KEY.strip():
            CURRENT_STATUS = f"Error: API key is required for {ai_provider}"
            print(f"VALIDATION ERROR: API key is required for {ai_provider}")
            return jsonify({'status': 'error', 'message': CURRENT_STATUS}), 400

    ai_client = AIModelClient(ai_provider, request.form['ai_model'], API_KEY)
    args = (request.form['pubmed_search'], request.form['screening_prompt'],
            request.form['screen_level'], ai_client, int(request.form.get('max_articles', 100000)))
    start_processing_thread(process_freeform_search, args)
    return jsonify({'status': 'started'})

@app.route('/status')
def get_status(): return jsonify({'status': CURRENT_STATUS})

@app.route('/stop', methods=['POST'])
def stop_processing():
    global STOP_REQUESTED
    STOP_REQUESTED = True
    return jsonify({'status': 'stop_requested'})

@app.route('/results')
def get_results():
    import glob

    downloads_dir = os.path.join(os.path.expanduser("~"), "Downloads")

    # Determine which type of results we're looking for
    if 'Relevant_Articles' in RESULTS_DF.columns:
        # Freeform search results
        pattern = os.path.join(downloads_dir, "screening_results_*.xlsx")
        file_list = glob.glob(pattern)
        if file_list:
            # Get the most recent file
            most_recent = max(file_list, key=os.path.getmtime)
            return send_file(most_recent, as_attachment=True)
    else:
        # Comparison screening results
        pattern = os.path.join(downloads_dir, "prompt_results_*.xlsx")
        file_list = glob.glob(pattern)
        if file_list:
            # Get the most recent file
            most_recent = max(file_list, key=os.path.getmtime)
            return send_file(most_recent, as_attachment=True)

    return (jsonify({'error': 'No results found. Please ensure processing has completed successfully.'}), 404)

@app.route('/export_intermediate')
def export_intermediate():
    global RESULTS_DF
    if RESULTS_DF.empty: return jsonify({'error': 'No results to export'}), 404
    
    if 'Relevant_Articles' in RESULTS_DF.columns and os.path.exists(DATABASE):
        with sqlite3.connect(DATABASE) as conn:
            pmids = pd.read_sql_query("SELECT pmid FROM Articles WHERE relevant = 1", conn)
            RESULTS_DF.at[0, 'Relevant_Articles'] = len(pmids)
            RESULTS_DF.at[0, 'Relevant_PMIDs'] = ",".join(pmids['pmid'].tolist())
    
    path = os.path.join(os.path.expanduser("~"), "Downloads", "intermediate_results.xlsx")
    RESULTS_DF.to_excel(path, index=False)
    return send_file(path, as_attachment=True)

@app.route('/debug_stream')
def debug_stream():
    return Response(get_debug_messages(), mimetype='text/event-stream')

@app.route('/get_terminal_output')
def get_terminal_output():
    output = get_captured_output()
    return jsonify({'output': output})

if __name__ == '__main__':
    if not getattr(sys, 'frozen', False) and not os.path.exists('templates'):
        os.makedirs('templates')
    setup_terminal_capture()
    start_output_monitoring()
    url = "http://127.0.0.1:5000"
    print(f"Server is running at {url}")
    threading.Timer(1.5, lambda: webbrowser.open(url)).start()
    app.run(debug=True, threaded=True, use_reloader=False)
