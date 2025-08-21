import os
import random
import sqlite3
import pandas as pd
from Bio import Entrez
from groq import Groq
from flask import Flask, render_template, request, jsonify, send_file
import threading
import time

app = Flask(__name__)

# Konfiguration und Status
DATABASE = 'ensure.sqlite'
ENTREZ_EMAIL = ''
GROQ_API_KEY = ''
CURRENT_STATUS = 'idle'
RESULTS_DF = pd.DataFrame()
PROCESS_THREAD = None
STOP_REQUESTED = False

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

def fetch_articles(pmids):
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

def save_to_db(articles):
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.executemany(
            "INSERT INTO Articles (pmid, title, abstract) VALUES (?, ?, ?)",
            articles
        )
        conn.commit()

def screen_articles(stage, prompt, batch_size=5):
    global CURRENT_STATUS, STOP_REQUESTED
    
    client = Groq(api_key=GROQ_API_KEY)
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

            if stage == 'initial':
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
                response = client.chat.completions.create(
                    model="llama3-8b-8192",
                    messages=[{"role": "user", "content": full_prompt}]
                )
                response_text = response.choices[0].message.content.strip()
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

# --- MODIFIKATION BEGINNT HIER ---
def process_prompts(prompts_df, initial_pmids, goldstandard_pmids):
    global CURRENT_STATUS, RESULTS_DF, STOP_REQUESTED
    
    RESULTS_DF = prompts_df.copy()
    articles = fetch_articles(initial_pmids)
    
    # Definiere die Ergebnisspalten, die wir erstellen oder überschreiben wollen
    metric_names = ['sensitivity', 'specificity', 'ppv', 'tp', 'fp', 'tn', 'fn']
    output_columns = ['Final_Relevant_PMIDs']
    for prefix in ['title', 'abstract']:
        for metric in metric_names:
            output_columns.append(f'{prefix}_{metric}')

    # Setze diese Spalten auf None, um alte Ergebnisse aus einer vorherigen Ausführung zu löschen
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
        
        # Titel-Screening
        if pd.notna(row.get('screen_titles')) and row['screen_titles'] == 1:
            screen_articles(stage='initial', prompt=row['TitlePrompt'])
            title_metrics = calculate_metrics(goldstandard_pmids, 'titles')
        else:
            with sqlite3.connect(DATABASE) as conn:
                conn.cursor().execute("UPDATE Articles SET stage = 'titles' WHERE stage = 'initial'")
                conn.commit()
            title_metrics = empty_metrics.copy()

        for metric, value in title_metrics.items():
            RESULTS_DF.at[idx, f'title_{metric}'] = value
        
        if STOP_REQUESTED: break
            
        # Abstract-Screening
        if pd.notna(row.get('screen_abstracts')) and row['screen_abstracts'] == 1:
            screen_articles(stage='titles', prompt=row['AbstractPrompt'])
            abstract_metrics = calculate_metrics(goldstandard_pmids, 'abstracts')
        else:
            abstract_metrics = empty_metrics.copy()
        
        for metric, value in abstract_metrics.items():
            RESULTS_DF.at[idx, f'abstract_{metric}'] = value
            
        # Sammle die finalen relevanten PMIDs und schreibe sie in die "Output"-Spalte
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
# --- MODIFIKATION ENDET HIER ---

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/run', methods=['POST'])
def run_screening():
    global ENTREZ_EMAIL, GROQ_API_KEY, CURRENT_STATUS, PROCESS_THREAD, STOP_REQUESTED
    
    STOP_REQUESTED = False
    ENTREZ_EMAIL = request.form['entrez_email']
    GROQ_API_KEY = request.form['groq_api_key']
    
    initial_file = request.files['initial_file']
    goldstandard_file = request.files['goldstandard_file']
    prompts_file = request.files['prompts_file']
    
    initial_content = initial_file.read().decode('utf-8')
    goldstandard_content = goldstandard_file.read().decode('utf-8')
    
    initial_pmids = [line.strip() for line in initial_content.splitlines() if line.strip()]
    goldstandard_pmids = [line.strip() for line in goldstandard_content.splitlines() if line.strip()]
    
    prompts_df = pd.read_excel(prompts_file)
    
    CURRENT_STATUS = "Starting screening process..."
    
    PROCESS_THREAD = threading.Thread(
        target=process_prompts,
        args=(prompts_df, initial_pmids, goldstandard_pmids)
    )
    PROCESS_THREAD.start()
    
    return jsonify({'status': 'started'})

@app.route('/status', methods=['GET'])
def get_status():
    return jsonify({'status': CURRENT_STATUS})

@app.route('/stop', methods=['POST'])
def stop_processing():
    global STOP_REQUESTED, PROCESS_THREAD
    
    STOP_REQUESTED = True
    
    if PROCESS_THREAD and PROCESS_THREAD.is_alive():
        PROCESS_THREAD.join(timeout=5.0)
    
    return jsonify({'status': 'stop_requested', 'message': CURRENT_STATUS})

@app.route('/results', methods=['GET'])
def get_results():
    global RESULTS_DF
    if not RESULTS_DF.empty:
        output_path = "prompt_results.xlsx"
        RESULTS_DF.to_excel(output_path, index=False)
        return send_file(output_path, as_attachment=True)
    return jsonify({'error': 'No results available'}), 404

if __name__ == '__main__':
    app.run(debug=True, threaded=True)