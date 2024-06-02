# How to use the Pipeline

# Preliminaries:

- Python & Interface (Spyder, VSCode, PyCharm, etc.)
- SQLite
- Account bei PubMed: [https://www.ncbi.nlm.nih.gov/home/develop/api/](https://www.ncbi.nlm.nih.gov/home/develop/api/)
- Zugang zum GitHub: https://github.com/BerlinExchangeMedicine/ensure
- Zugang zu Google Drive: https://drive.google.com/drive/folders/1etCiSCeKktkCMuM-ADni9RSTgdkdY3AK

---

# Skripte und Benutzung:

— Skripte in folgender Reihenfolge ausführen

### Set-Up:

1. “**utils”-Ordner: dbconnect.py** 
    
    → Einfach laufen lassen
    
    → Output: “Database initialized with the necessary tables.” 
    
2. **“utils”-Ordner: PubMed Extraction Txt-file Based.py**
    
    → Einfach laufen lassen
    
    → Output: 
    “Total number of papers in the database: 493 
    Excel file 'PubMed_Data.xlsx' created.”
    

### Prompt Engineering direkt:

1. **02_title_screening.py**
    1. Benötigt API-Key → einfach nachfragen
        
        ```python
        client = OpenAI(api_key='upon request')
        ```
        
    2. Prompt Engineering spielt sich in folgenden Zeilen ab:
        
        ```python
        # Funktion, um ein Prompt für OpenAI zu erstellen
        def generate_prompt(search_term, titles):
            prompt = f"Search Term: {search_term}\n\nList of titles:\n"
            prompt += "\n".join(f"{i}. {title}" for i, title in enumerate(titles, 1))
            
            #PROMPT ENGINEERING
            prompt += "\n ..."
            prompt += "\n ..."
            
            # Formatierung für GPT (NICHT ÄNDERN)
            prompt += "\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes. "
            prompt += "Each entry should be the first 10 characters of each relevant title. "
            prompt += "For example: 'Example AB', 'Example CD', 'Example EF'."
            return prompt
        ```
        
        → unter “#PROMPT ENGINEERING” beliebige veränderungen vornehmen und mit “prompt +=” Prompt erweitern
        
        → Note: \n setzt einen Absatz ;)
        
    3. Skript ausführen → Output: “Enter Search Term: “
        
        ⇒ Hier die Search terms der Meta-Analyse einfügen
        
    4. Ausführen und laufen lassen
2. **03_abstract_screening.py**
    1. Benötigt ebenfalls API-Key
        
        ```python
        client = OpenAI(api_key='upon request')
        ```
        
    2. Prompt Engineering spielt sich in folgenden Zeilen ab:
        
        ```python
        # Funktion, um ein Prompt für OpenAI zu erstellen
        def generate_prompt(search_term, abstracts):
            prompt = f"Search Term: {search_term}\n\nList of abstracts:\n"
            prompt += "\n".join(abstract for _, abstract in abstracts)
            
            #PROMPT ENGINEERING
            prompt += "\n ..."
            prompt += "\n ..."
            
            # Formatierung für GPT (NICHT ÄNDERN)
            prompt += "\n\nFormat the output as a comma-separated string with each entry enclosed in single quotes. "
            prompt += "Each entry should be the first 10 characters of each relevant abstract. "
            prompt += "For example: 'Example AB', 'Example CD', 'Example EF'."
            return prompt
        ```
        
        → gleiches Spiel
        

### Evaluation & Repitition:

1. “**utils”-Ordner: Performance.py**
    
    → Output: Sensitivity & Specificity Werte 
    
    → Mit Prompt Excel Tabelle des Drive eintragen
    
    → Note: Diese Performance könnt ihr euch an jedem Schritt des Screenings geben lassen — also nach Title und Abstract jeweils
    
    ![Untitled](How%20to%20use%20the%20Pipeline%2000252f9536f24e94a4d7f2fdd7d4e868/Untitled.png)
    
2. **“utils”-Ordner: Screening Progression**
    
    → Outputs: 
    
    - Diagramm zu Rohwerten der Paper zu jedem Schritt
    - Excel Tabelle für manuelle Inspektion
3. **“utils”-Ordner: Clear_Databases.py**
    
    → Outputs: All tables have been cleared.
    
    → Nach jedem Durchlauf alle Tabellen clearen und neu beginnen
    

### ⇒ Für jeden guten Prompt 5 Mal durchgehen und Mittelwert bilden