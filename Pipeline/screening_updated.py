import os
import random
import sqlite3
import pandas as pd
import json
import re
from Bio import Entrez
from groq import Groq
from openai import OpenAI
from anthropic import Anthropic
import google.generativeai as genai
from ollama import Client as OllamaClient
from flask import Flask, render_template, request, jsonify, send_file
import threading
import time

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

class AIModelClient:
    def __init__(self, provider, model_name, api_key):
        self.provider = provider.lower()
        self.model_name = model_name
        self.api_key = api_key
        self.supported_providers = ['groq', 'openai', 'anthropic', 'google', 'ollama']

        if self.provider not in self.supported_providers:
            supported_list = ', '.join(self.supported_providers)
            raise ValueError(f"Provider '{provider}' wird nicht unterstützt. Unterstützte Provider: {supported_list}. "
                           "Bei Bedarf können neue Provider durch kleine Code-Anpassungen hinzugefügt werden.")

        self.client = self._init_client()

    def _init_client(self):
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
            raise ValueError(f"Unbekannter Provider: {self.provider}")

    def generate_completion(self, prompt):
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
            error_msg = f"Fehler bei der Kommunikation mit {self.provider} {self.model_name}: {str(e)}"
            if "API key" in str(e).lower() or "key" in str(e).lower():
                error_msg += " - Überprüfen Sie Ihren API-Key."
            elif "model" in str(e).lower():
                error_msg += " - Das Modell könnte nicht verfügbar sein."
            raise Exception(error_msg)

def init_db():
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

def fetch_articles_by_pmids(pmids):
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

def search_pubmed(query, retmax=1000):
    """Search PubMed and return list of PMIDs."""
    Entrez.email = ENTREZ_EMAIL
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    results = Entrez.read(handle)
    handle.close()
    return results["IdList"]

def save_to_db(articles):
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.executemany(
            "INSERT INTO Articles (pmid, title, abstract) VALUES (?, ?, ?)",
            articles
        )
        conn.commit()

def screen_articles(stage, prompt, ai_client, screen_level='both', batch_size=5):
    global CURRENT_STATUS, STOP_REQUESTED

    results = []
    print(f"DEBUG: Starting screen_articles with stage={stage}, screen_level={screen_level}")

    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT pmid, title, abstract FROM Articles WHERE stage = ?",
            (stage,)
        )
        articles = cursor.fetchall()
        total = len(articles)
        print(f"DEBUG: Found {total} articles at stage {stage}")

        for i in range(0, total, batch_size):
            if STOP_REQUESTED:
                CURRENT_STATUS = "Screening stopped by user"
                return results

            if stage == 'initial' or screen_level == 'titles':
                CURRENT_STATUS = f"Screening titles ({min(i+batch_size, total)}/{total})"
            else:
                CURRENT_STATUS = f"Screening abstracts ({min(i+batch_size, total)}/{total})"
            batch = articles[i:i+batch_size]
            print(f"DEBUG: Processing batch {i//batch_size + 1} with {len(batch)} articles")

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
                print("=== DEBUG: Full Prompt Sent to AI ===")
                print(full_prompt)
                print("=====================================")

                response_text = ai_client.generate_completion(full_prompt).strip()

                print("=== DEBUG: Full AI Response ===")
                print(f"'{response_text}'")
                print("===============================")

                # Extract all potential PMIDs from the response using regex (7-9 digits)
                extracted_pmids = re.findall(r'\b\d{7,9}\b', response_text)
                print(f"DEBUG: Extracted PMIDs from response: {extracted_pmids}")

                # Get the list of PMIDs in the current batch
                batch_pmids = [pmid for pmid, _, _ in batch]
                print(f"DEBUG: Batch PMIDs: {batch_pmids}")

                # Only consider PMIDs that are actually in the batch
                relevant_pmids = [pmid for pmid in extracted_pmids if pmid in batch_pmids]
                print(f"DEBUG: Relevant PMIDs (filtered): {relevant_pmids}")

                # Update database even if AI call fails - mark all as not relevant and advance stage
                for pmid, _, _ in batch:
                    relevant = 1 if pmid in relevant_pmids else 0
                    next_stage = "titles" if stage == "initial" else "abstracts"
                    cursor.execute(
                        "UPDATE Articles SET relevant = ?, stage = ? WHERE pmid = ?",
                        (relevant, next_stage, pmid)
                    )
                    if relevant:
                        results.append(pmid)
                print(f"DEBUG: Batch results: {len([r for r in results if r in batch_pmids])} relevant articles")

                # Commit after successful processing
                conn.commit()
                print("DEBUG: Database committed successfully")

            except Exception as e:
                print(f"DEBUG: Error during screening: {str(e)}")
                CURRENT_STATUS = f"Error during screening: {str(e)}"

                # Even on error, we need to advance the stage for these articles
                # so they don't get stuck in the current stage
                print("DEBUG: Advancing stage for batch articles despite error...")
                for pmid, _, _ in batch:
                    next_stage = "titles" if stage == "initial" else "abstracts"
                    cursor.execute(
                        "UPDATE Articles SET relevant = 0, stage = ? WHERE pmid = ?",
                        (next_stage, pmid)
                    )
                conn.commit()
                print("DEBUG: Database committed after error handling")
                continue

            time.sleep(1)

        conn.commit()

    print(f"DEBUG: Total results for stage {stage}: {len(results)}")
    return results

def calculate_metrics(goldstandard, stage):
    print(f"DEBUG: Calculating metrics for stage '{stage}'")
    print(f"DEBUG: Goldstandard PMIDs: {goldstandard}")

    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT pmid, relevant FROM Articles WHERE stage = ?",
            (stage,)
        )
        data = cursor.fetchall()
    print(f"DEBUG: Retrieved {len(data)} articles from database at stage '{stage}'")

    if not data:
        print("DEBUG: No data found in database for this stage!")
        return {
            'sensitivity': 0, 'specificity': 0, 'ppv': 0,
            'tp': 0, 'fp': 0, 'tn': 0, 'fn': 0
        }

    df = pd.DataFrame(data, columns=["pmid", "relevant"])
    print(f"DEBUG: DataFrame shape: {df.shape}")
    print(f"DEBUG: Relevant counts: {df['relevant'].value_counts().to_dict()}")

    # Debug: Show which PMIDs are in database vs goldstandard
    all_db_pmids = set(df['pmid'])
    gs_set = set(goldstandard)
    print(f"DEBUG: PMIDs in database: {sorted(all_db_pmids)}")
    print(f"DEBUG: PMIDs in goldstandard: {sorted(gs_set)}")
    print(f"DEBUG: Overlap (should be {len(gs_set)}): {len(all_db_pmids.intersection(gs_set))}")
    print(f"DEBUG: Missing from database: {gs_set.difference(all_db_pmids)}")
    print(f"DEBUG: Extra in database: {all_db_pmids.difference(gs_set)}")

    screened_in = set(df[df['relevant'] == 1]['pmid'])
    screened_out = set(df[df['relevant'] == 0]['pmid'])
    print(f"DEBUG: Screened in (relevant=1): {sorted(screened_in)}")
    print(f"DEBUG: Screened out (relevant=0): {sorted(screened_out)[:10]}...")  # First 10 only

    tp = len(gs_set.intersection(screened_in))
    fp = len(screened_in.difference(gs_set))
    fn = len(gs_set.intersection(screened_out))
    all_pmids = set(df['pmid'])
    tn = len(all_pmids.difference(gs_set).difference(screened_in))

    print(f"DEBUG: TP={tp}, FP={fp}, TN={tn}, FN={fn}")
    print(f"DEBUG: Goldstandard articles that were missed: {sorted(gs_set.intersection(screened_out))}")

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0

    print(f"DEBUG: Final metrics - Sensitivity: {sensitivity:.3f}, Specificity: {specificity:.3f}, PPV: {ppv:.3f}")

    return {
        'sensitivity': sensitivity,
        'specificity': specificity,
        'ppv': ppv,
        'tp': tp,
        'fp': fp,
        'tn': tn,
        'fn': fn
    }

def process_prompts(prompts_df, initial_pmids, goldstandard_pmids, ai_client):
    global CURRENT_STATUS, RESULTS_DF, STOP_REQUESTED

    RESULTS_DF = prompts_df.copy()

    # Use all initial_pmids for screening, goldstandard is only for metrics
    goldstandard_set = set(goldstandard_pmids)
    filtered_initial_pmids = [pmid for pmid in initial_pmids if pmid in goldstandard_set]

    print(f"DEBUG: Original initial_pmids count: {len(initial_pmids)}")
    print(f"DEBUG: Goldstandard_pmids count: {len(goldstandard_pmids)}")
    print(f"DEBUG: Filtered initial_pmids count: {len(filtered_initial_pmids)}")
    print(f"DEBUG: Filtered PMIDs: {filtered_initial_pmids}")

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
        
        print(f"DEBUG: Processing prompt {idx+1}")
        print(f"DEBUG: screen_titles value: {row.get('screen_titles')}, type: {type(row.get('screen_titles'))}")
        print(f"DEBUG: screen_abstracts value: {row.get('screen_abstracts')}, type: {type(row.get('screen_abstracts'))}")
        print(f"DEBUG: TitlePrompt: {row.get('TitlePrompt')}")
        print(f"DEBUG: AbstractPrompt: {row.get('AbstractPrompt')}")

        # Convert to int for reliable comparison
        screen_titles_val = int(row['screen_titles']) if pd.notna(row.get('screen_titles')) else 0
        screen_abstracts_val = int(row['screen_abstracts']) if pd.notna(row.get('screen_abstracts')) else 0

        if screen_titles_val == 1:
            print("DEBUG: Screening titles...")
            screen_articles(stage='initial', prompt=row['TitlePrompt'], ai_client=ai_client)
            title_metrics = calculate_metrics(goldstandard_pmids, 'titles')
            print(f"DEBUG: Title metrics: {title_metrics}")
        else:
            print("DEBUG: Skipping title screening")
            with sqlite3.connect(DATABASE) as conn:
                conn.cursor().execute("UPDATE Articles SET stage = 'titles' WHERE stage = 'initial'")
                conn.commit()
            title_metrics = empty_metrics.copy()

        for metric, value in title_metrics.items():
            RESULTS_DF.at[idx, f'title_{metric}'] = value

        if STOP_REQUESTED: break

        if screen_abstracts_val == 1:
            print("DEBUG: Screening abstracts...")
            screen_articles(stage='titles', prompt=row['AbstractPrompt'], ai_client=ai_client)
            abstract_metrics = calculate_metrics(goldstandard_pmids, 'abstracts')
            print(f"DEBUG: Abstract metrics: {abstract_metrics}")
        else:
            print("DEBUG: Skipping abstract screening")
            abstract_metrics = empty_metrics.copy()
        
        for metric, value in abstract_metrics.items():
            RESULTS_DF.at[idx, f'abstract_{metric}'] = value
            
        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT pmid FROM Articles WHERE relevant = 1")
            final_pmids_list = [item[0] for item in cursor.fetchall()]
            RESULTS_DF.at[idx, 'Final_Relevant_PMIDs'] = ",".join(final_pmids_list)
            print(f"DEBUG: Final relevant PMIDs for prompt {idx+1}: {final_pmids_list}")

        if STOP_REQUESTED: break

    print(f"DEBUG: Final RESULTS_DF shape: {RESULTS_DF.shape}")
    print(f"DEBUG: Final RESULTS_DF columns: {list(RESULTS_DF.columns)}")
    print("DEBUG: Final RESULTS_DF content:")
    print(RESULTS_DF.to_string())

    output_path = "prompt_results.xlsx"
    RESULTS_DF.to_excel(output_path, index=False)
    print(f"DEBUG: Results saved to {output_path}")

    if STOP_REQUESTED:
        CURRENT_STATUS = f"Stopped! Partial results saved to {output_path}"
    else:
        CURRENT_STATUS = f"Completed! Results saved to {output_path}"

def process_freeform_search(pubmed_query, screening_prompt, screen_level, ai_client, max_articles=1000):
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
        CURRENT_STATUS = "Screening abstracts of relevant titles..."
        screen_articles('titles', screening_prompt, ai_client, 'abstracts')
        # Get final relevant PMIDs
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
    return render_template('index_updated.html')

def start_processing_thread(target_func, args_tuple):
    global CURRENT_STATUS, PROCESS_THREAD, STOP_REQUESTED
    STOP_REQUESTED = False
    CURRENT_STATUS = "Starting screening process..."
    PROCESS_THREAD = threading.Thread(target=target_func, args=args_tuple)
    PROCESS_THREAD.start()

@app.route('/run_comparison', methods=['POST'])
def run_comparison():
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
    print("=== DEBUG: Initial.txt Content ===")
    print(f"Raw content length: {len(initial_content)}")
    print(f"First 500 chars: {initial_content[:500]}")
    print(f"Total PMIDs found: {len(initial_pmids)}")
    print(f"First 10 PMIDs: {initial_pmids[:10]}")
    print(f"Last 10 PMIDs: {initial_pmids[-10:] if len(initial_pmids) > 10 else initial_pmids}")
    print("===================================")

    # Read and debug goldstandard file
    goldstandard_content = goldstandard_file.read().decode('utf-8')
    goldstandard_pmids = [line.strip() for line in goldstandard_content.splitlines() if line.strip()]
    print("=== DEBUG: Goldstandard_Selected.txt Content ===")
    print(f"Raw content length: {len(goldstandard_content)}")
    print(f"Full content: {goldstandard_content}")
    print(f"Total PMIDs found: {len(goldstandard_pmids)}")
    print(f"All PMIDs: {goldstandard_pmids}")
    print("================================================")

    # Read and debug prompts file
    prompts_df = pd.read_excel(prompts_file)
    print("=== DEBUG: Prompts.xlsx Content ===")
    print(f"Shape: {prompts_df.shape}")
    print(f"Columns: {list(prompts_df.columns)}")
    print(f"Data types:\n{prompts_df.dtypes}")
    print(f"screen_titles unique values: {prompts_df['screen_titles'].unique()}")
    print(f"screen_abstracts unique values: {prompts_df['screen_abstracts'].unique()}")
    print("\nFirst 3 rows:")
    print(prompts_df.head(3).to_string())
    print("===================================")

    ai_client = AIModelClient(ai_provider, ai_model, API_KEY)

    start_processing_thread(process_prompts, (prompts_df, initial_pmids, goldstandard_pmids, ai_client))
    return jsonify({'status': 'started'})

@app.route('/run_freeform', methods=['POST'])
def run_freeform():
    global ENTREZ_EMAIL, API_KEY, SELECTED_MODEL
    ENTREZ_EMAIL = request.form['entrez_email']
    API_KEY = request.form['api_key']
    ai_provider = request.form['ai_provider']
    ai_model = request.form['ai_model']
    SELECTED_MODEL = f"{ai_provider}-{ai_model}"

    pubmed_search = request.form['pubmed_search']
    screening_prompt = request.form['screening_prompt']
    screen_level = request.form['screen_level']
    max_articles = int(request.form.get('max_articles', 1000))

    ai_client = AIModelClient(ai_provider, ai_model, API_KEY)

    start_processing_thread(process_freeform_search, (pubmed_search, screening_prompt, screen_level, ai_client, max_articles))
    return jsonify({'status': 'started'})

@app.route('/status')
def get_status():
    return jsonify({'status': CURRENT_STATUS})

@app.route('/stop', methods=['POST'])
def stop_processing():
    global STOP_REQUESTED, PROCESS_THREAD
    STOP_REQUESTED = True
    if PROCESS_THREAD and PROCESS_THREAD.is_alive():
        PROCESS_THREAD.join(timeout=5.0)
    return jsonify({'status': 'stop_requested', 'message': CURRENT_STATUS})

@app.route('/results')
def get_results():
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

def test_pmid_extraction():
    """Test function to verify PMID extraction from AI responses."""
    print("=== PMID Extraction Test ===")

    # Test cases with different response formats
    test_cases = [
        {
            "response": "12345678,23456789,34567890",
            "batch_pmids": ["12345678", "23456789", "34567890", "45678901"],
            "expected": ["12345678", "23456789", "34567890"]
        },
        {
            "response": "The relevant PMIDs are: 12345678, 23456789 and 34567890",
            "batch_pmids": ["12345678", "23456789", "34567890", "45678901"],
            "expected": ["12345678", "23456789", "34567890"]
        },
        {
            "response": "Based on the criteria, I found these relevant: 12345678,23456789",
            "batch_pmids": ["12345678", "23456789", "34567890"],
            "expected": ["12345678", "23456789"]
        },
        {
            "response": "No relevant articles found.",
            "batch_pmids": ["12345678", "23456789", "34567890"],
            "expected": []
        },
        {
            "response": "12345678, 99999999, 23456789",  # Contains invalid PMID
            "batch_pmids": ["12345678", "23456789", "34567890"],
            "expected": ["12345678", "23456789"]
        }
    ]

    for i, test_case in enumerate(test_cases, 1):
        print(f"\nTest Case {i}:")
        print(f"Response: {test_case['response']}")
        print(f"Batch PMIDs: {test_case['batch_pmids']}")

        # Simulate the extraction logic
        extracted_pmids = re.findall(r'\b\d{7,9}\b', test_case['response'])
        relevant_pmids = [pmid for pmid in extracted_pmids if pmid in test_case['batch_pmids']]

        print(f"Extracted PMIDs: {extracted_pmids}")
        print(f"Relevant PMIDs: {relevant_pmids}")
        print(f"Expected: {test_case['expected']}")
        print(f"✓ PASS" if relevant_pmids == test_case['expected'] else "✗ FAIL")

    print("\n=== Test Complete ===")

@app.route('/test_parsing')
def test_parsing():
    """Web endpoint to run PMID extraction tests."""
    import io
    import sys

    # Capture print output
    old_stdout = sys.stdout
    sys.stdout = buffer = io.StringIO()

    try:
        test_pmid_extraction()
        output = buffer.getvalue()
    finally:
        sys.stdout = old_stdout

    return f"<pre>{output}</pre>"

if __name__ == '__main__':
    if not os.path.exists('templates'):
        os.makedirs('templates')

    # Check if we should run tests
    if len(os.sys.argv) > 1 and os.sys.argv[1] == 'test':
        test_pmid_extraction()
    else:
        app.run(debug=True, threaded=True)
