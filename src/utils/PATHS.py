import os

#Hier den API Key eingeben
api_key = os.getenv('OPENAI_API_KEY_NOS') # oder hier openai key von till

if api_key is None:
    raise ValueError("OPENAI_API_KEY environment variable is not set")

#Hier die Enrez Mail eingeben
email='schrott.ew@t-online.de'
