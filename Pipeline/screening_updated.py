import os
import random
import sqlite3
import pandas as pd
import json
from Bio import Entrez
from groq import Groq
from openai import OpenAI
from anthropic import Anthropic
import google.generativeai as genai
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
    def __init__(self, model_type, api_key):
        self.model_type = model_type
        self.api_key = api_key
        self.supported_models = {
            'groq-llama3': {'provider': 'groq', 'model_name': 'llama3-8b-8192'},
            'openai-gpt4': {'provider': 'openai', 'model_name': 'gpt-4-turbo-preview'},
            'openai-gpt35': {'provider': 'openai', 'model_name': 'gpt-3.5-turbo'},
            'anthropic-claude': {'provider': 'anthropic', 'model_name': 'claude-3-opus-20240229'},
            'google-gemini': {'provider': 'google', 'model_name': 'gemini-1.0-pro'}
        }

        if model_type not in self.supported_models:
            supported_list = ', '.join(self.supported_models.keys())
            raise ValueError(f"Modell '{model_type}' wird nicht unterstützt. Unterstützte Modelle: {supported_list}. "
                           f"Bei Bedarf können neue Modelle durch kleine Code-Anpassungen hinzugefügt werden.")

        self.model_config = self.supported_models[model_type]
        self.client = self._init_client()

    def _init_client(self):
        provider = self.model_config['provider']
        if provider == 'groq':
            return Groq(api_key=self.api_key)
        elif provider == 'openai':
            return OpenAI(api_key=self.api_key)
        elif provider == 'anthropic':
            return Anthropic(api_key=self.api_key)
        elif provider == 'google':
            genai.configure(api_key=self.api_key)
            return genai
        else:
            raise ValueError(f"Unbekannter Provider: {provider}")

    def generate_completion(self, prompt):
        try:
            provider = self.model_config['provider']
            model_name = self.model_config['model_name']

            if provider == 'groq':
                response = self.client.chat.completions.create(
                    model=model_name,
                    messages=[{"role": "user", "content": prompt}]
                )
                return response.choices[0].message.content

            elif provider == 'openai':
                response = self.client.chat.completions.create(
                    model=model_name,
                    messages=[{"role": "user", "content": prompt}]
                )
                return response.choices[0].message.content

            elif provider == 'anthropic':
                response = self.client.messages.create(
                    model=model_name,
                    messages=[{"role": "user", "content": prompt}]
                )
                return response.content[0].text

            elif provider == 'google':
                model = self.client.GenerativeModel(model_name)
                response = model.generate_content(prompt)
                return response.text

        except Exception as e:
            error_msg = f"Fehler bei der Kommunikation mit {self.model_type}: {str(e)}"
            if "API key" in str(e).lower():
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

def search_pubmed(query):
    """Search PubMed and return list of PMIDs."""
    Entrez.email = ENTREZ_EMAIL
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1000)
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
    
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT pmid, title, abstract FROM Articles WHERE stage = ?", 
            (stage,)
        )
        articles = cursor.fetchall()
        total = len(articles)
        
        for i in range(0, total, batch_size):
            if STOP_REQUESTED:
                CURRENT_STATUS = "Screening stopped by user"
                return results
                
            CURRENT_STATUS = f"Screening {stage} ({min(i+batch_size, total)}/{total})"
            batch = articles[i:i+batch_size]

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
                f"For example, your response must be exactly: 12345,67890,78912"
            )

            try:
                response_text = ai_client.generate_completion(full_prompt).strip()
                relevant_pmids = [p.strip() for p in response_text.split(",") if p.strip()]

                for pmid, _, _ in batch:
                    relevant = 1 if pmid in relevant_pmids else 0
                    next_stage = "titles" if stage == "initial" else "abstracts"
                    cursor.execute(
                        "UPDATE Articles SET relevant = ?, stage = ? WHERE pmid = ?",
                        (relevant, next_stage, pmid)
                    )
                    if relevant:
                        results.append(pmid)
            except Exception as e:
                CURRENT_STATUS = f"Error during screening: {str(e)}"
                continue
            
            time.sleep(1)
        
        conn.commit()
    
    return results

def calculate_metrics(goldstandard, stage):
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT pmid, relevant FROM Articles WHERE stage = ?",
            (stage,)
        )
        data = cursor.fetchall()
    
    df = pd.DataFrame(data, columns=["pmid", "relevant"])
    gs_set = set(goldstandard)
    
    screened_in = set(df[df['relevant'] == 1]['pmid'])
    screened_out = set(df[df['relevant'] == 0]['pmid'])

    tp = len(gs_set.intersection(screened_in))
    fp = len(screened_in.difference(gs_set))
    fn = len(gs_set.intersection(screened_out))
    all_pmids = set(df['pmid'])
    tn = len(all_pmids.difference(gs_set).difference(screened_in))

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    
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
        
        if pd.notna(row.get('screen_titles')) and row['screen_titles'] == 1:
            screen_articles(stage='initial', prompt=row['TitlePrompt'], ai_client=ai_client)
            title_metrics = calculate_metrics(goldstandard_pmids, 'titles')
        else:
            with sqlite3.connect(DATABASE) as conn:
                conn.cursor().execute("UPDATE Articles SET stage = 'titles' WHERE stage = 'initial'")
                conn.commit()
            title_metrics = empty_metrics.copy()

        for metric, value in title_metrics.items():
            RESULTS_DF.at[idx, f'title_{metric}'] = value
        
        if STOP_REQUESTED: break
            
        if pd.notna(row.get('screen_abstracts')) and row['screen_abstracts'] == 1:
            screen_articles(stage='titles', prompt=row['AbstractPrompt'], ai_client=ai_client)
            abstract_metrics = calculate_metrics(goldstandard_pmids, 'abstracts')
        else:
            abstract_metrics = empty_metrics.copy()
        
        for metric, value in abstract_metrics.items():
            RESULTS_DF.at[idx, f'abstract_{metric}'] = value
            
        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT pmid FROM Articles WHERE relevant = 1")
            final_pmids_list = [item[0] for item in cursor.fetchall()]
            RESULTS_DF.at[idx, 'Final_Relevant_PMIDs'] = ",".join(final_pmids_list)

        if STOP_REQUESTED: break
            
    output_path = "prompt_results.xlsx"
    RESULTS_DF.to_excel(output_path, index=False)
    
    if STOP_REQUESTED:
        CURRENT_STATUS = f"Stopped! Partial results saved to {output_path}"
    else:
        CURRENT_STATUS = f"Completed! Results saved to {output_path}"

def process_freeform_search(pubmed_query, screening_prompt, screen_level, ai_client):
    global CURRENT_STATUS, RESULTS_DF, STOP_REQUESTED
    
    CURRENT_STATUS = "Searching PubMed..."
    pmids = search_pubmed(pubmed_query)
    
    if not pmids:
        CURRENT_STATUS = "No articles found matching the PubMed query"
        return
    
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
    
    # Create results DataFrame
    RESULTS_DF = pd.DataFrame({
        'PubMed_Query': [pubmed_query],
        'Screening_Prompt': [screening_prompt],
        'Screen_Level': [screen_level],
        'Total_Articles': [len(pmids)],
        'Relevant_Articles': [len(results)],
        'Relevant_PMIDs': [",".join(results)]
    })
    
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
    SELECTED_MODEL = request.form['ai_model']
    
    initial_file = request.files['initial_file']
    goldstandard_file = request.files['goldstandard_file']
    prompts_file = request.files['prompts_file']
    
    initial_pmids = [line.strip() for line in initial_file.read().decode('utf-8').splitlines() if line.strip()]
    goldstandard_pmids = [line.strip() for line in goldstandard_file.read().decode('utf-8').splitlines() if line.strip()]
    prompts_df = pd.read_excel(prompts_file)
    
    ai_client = AIModelClient(SELECTED_MODEL, API_KEY)
    
    start_processing_thread(process_prompts, (prompts_df, initial_pmids, goldstandard_pmids, ai_client))
    return jsonify({'status': 'started'})

@app.route('/run_freeform', methods=['POST'])
def run_freeform():
    global ENTREZ_EMAIL, API_KEY, SELECTED_MODEL
    ENTREZ_EMAIL = request.form['entrez_email']
    API_KEY = request.form['api_key']
    SELECTED_MODEL = request.form['ai_model']
    
    pubmed_search = request.form['pubmed_search']
    screening_prompt = request.form['screening_prompt']
    screen_level = request.form['screen_level']
    
    ai_client = AIModelClient(SELECTED_MODEL, API_KEY)
    
    start_processing_thread(process_freeform_search, (pubmed_search, screening_prompt, screen_level, ai_client))
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
        output_path = "prompt_results.xlsx"
        if not os.path.exists(output_path):
            output_path = "screening_results.xlsx"
        if os.path.exists(output_path):
            return send_file(output_path, as_attachment=True)
    return jsonify({'error': 'No results available'}), 404

@app.route('/export_intermediate')
def export_intermediate():
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT pmid FROM Articles WHERE relevant = 1")
        pmids = [item[0] for item in cursor.fetchall()]
    df = pd.DataFrame({'Relevant_PMIDs': pmids})
    output_path = "intermediate_results.xlsx"
    df.to_excel(output_path, index=False)
    return send_file(output_path, as_attachment=True)

if __name__ == '__main__':
    if not os.path.exists('templates'):
        os.makedirs('templates')
    app.run(debug=True, threaded=True)
