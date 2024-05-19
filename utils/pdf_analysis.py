#
# USAGE EXAMPLE:
# - put PDF into pdfs subdirectory
# $ python pdf_analysis.py 'api-key' 'test.pdf'
#


import fitz  # PyMuPDF
import argparse
from chatgpt_api import chat

def extract_text_from_pdf(pdf_path):
    document = fitz.open(pdf_path)
    text = ""
    for page_num in range(len(document)):
        page = document.load_page(page_num)
        text += page.get_text()
    return text



parser = argparse.ArgumentParser(description="Extract text from a PDF and interact with ChatGPT")
parser.add_argument("api_key", type=str, help="OpenAI API key")
parser.add_argument("pdf_path", type=str, help="OpenAI API key")

args = parser.parse_args()

#print(args.api_key)

pdf_path = "pdfs/%s" % args.pdf_path
pdf_text = extract_text_from_pdf(pdf_path)

prompt = 'Erlaeutere den folgenden Text: ' + pdf_text

response = chat(prompt, args.api_key)

print(response)