# MA-LLM Pipeline — README

This repository contains the MA-LLM screening pipeline, a tool for automated screening of PubMed articles using Large Language Models (LLMs). It supports two modes:

- Screening Selection Comparison: compare LLM-based selections with a manual gold standard.
- Comparison-free Screening (Freeform): run PubMed searches and screen articles without a comparison set.

This README explains how to run the project from source, notes on the provided .exe (if any), required Python packages, common troubleshooting steps, and recommended small fixes and naming conventions.

## Table of contents
- Project Structure
- Quick start (source / ZIP)
- Running the Flask UI
- Troubleshooting (common errors including HTTP 500)
- Notes and recommendations


## Project structure (important files/folders)

- `MA-LLM Pipeline/2025-10-02_screening.py` — main Flask app and processing logic (entrypoint)
- `MA-LLM Pipeline/templates/2025-09-25_index.html` — front-end UI used by the Flask app
- `MA-LLM Pipeline/ExampleFiles/` — example PMIDs, prompts and gold-standard files (use these to test input formats)
- `requirements.txt` — Python package list
- `Readme.md`

## Quick start (from source / ZIP)

1. Clone or extract the ZIP and open a terminal in the repository root.
2. Create and activate a Python virtual environment (recommended):

```bash
python3 -m venv .venv
source .venv/bin/activate    # zsh/bash
```

3. Install dependencies:

```bash
pip install -r requirements.txt
```

4. Start the Flask UI (the web UI will open automatically):

```bash
python "MA-LLM Pipeline/2025-10-02_screening.py"
```

#### Running the provided .exe (Windows)

- If using .exe is provided in the release, it may not run out-of-the-box depending on how it was built. Ensure the file was built correctly for your platform and that the download didn't corrupt it.

#### Required Python packages

- The project expects a set of Python packages. Ensure `requirements.txt` includes at least the following (add or pin versions as needed):

- pandas
- biopython
- flask
- openai (if using OpenAI provider)
- anthropic (if using Anthropic provider)
- google-generative-ai (if using Google provider)
- ollama (if using Ollama)

## How the web UI submits work

- The front-end sends a POST to either `/run_comparison` (goldstandard mode) or `/run_freeform` (freeform mode).
- The server starts processing in a background thread and the front-end polls `/status` for progress.

## Common problems and troubleshooting

- HTTP 500 "Submission failed: HTTP error! status: 500"
	- Meaning: the server raised an unhandled exception while processing the submission. This is a server-side error, not a front-end problem.
	- What to do: open your browser DevTools → Network, find the POST to `/run_comparison` or `/run_freeform`, and inspect the Response body. The server now returns JSON with `message` and `traceback` fields to help debugging.
	- Likely causes in this codebase:
		- Required form fields or uploaded files were missing (the server expects `initial_file`, `goldstandard_file`, and `prompts_file` for comparison mode).
		- `prompts_file` is not an Excel file or has unexpected columns; `pd.read_excel` will raise on invalid input.
		- AI provider initialization failed (missing provider library, unsupported provider name, or authentication error).
		- If using `Ollama (Local)`, no API key should be required; the UI currently asks for an API key for all providers — leave the API key blank for Ollama or update the UI/server to skip the key requirement for Ollama.

- Unclear Screening Mode error message
	- If you get an error about not choosing a screening mode, the UI will show a message. Make sure `Screening Mode` is set to either `Screening Selection Comparison` (goldstandard) or `Comparison-free Screening` (freeform) before submitting.

- Problems reading `prompts.xlsx`
	- The pipeline expects columns named like `TitlePrompt`, `AbstractPrompt`, `screen_titles`, and `screen_abstracts` in the prompts Excel file. If your file uses different column names, rename them or adapt the code.

### Frontend & server behavior notes

- The UI requires `initial_file` and `goldstandard_file` (text files with PMIDs) and `prompts_file` (Excel) when you select the Comparison mode. Make sure files are uploaded in the form.
- For Freeform mode you need to provide a PubMed search query, screening prompt, and max articles.
- For the Ollama provider the server constructs a local client (hosted at `http://localhost:11434`) and should not require an API key. If the UI still shows the API key field for Ollama, you can safely leave it empty and start the server.

## License & citation