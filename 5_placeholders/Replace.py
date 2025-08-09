import pandas as pd

# Define file paths
input_file = "Prompts final.xlsx"
output_file = "Prompts.xlsx"

# Define replacement mappings
replacements = {
    r"\[ttl\]": "Preoperative radiotherapy combined with surgery versus surgery alone for primary retroperitoneal sarcoma: a meta-analysis",
    r"\[Inclusion Criteria\]": "Inclusion: The target of the search was to obtain articles that met the following inclusion criteria: direct comparison between SA and preRT + S for RPS; assessment of LR, ARFS, overall survival (OS), and treatment-related complications; and result presentation in a form that allowed quantitative synthesis.",
    r"\[Exclusion Criteria\]": "Exclusion: We exclude articles for the following reasons: inclusion of other treatments, such as chemotherapy and hyperthermia; no comparative studies regarding the two treatment arms; recurrent tumors; regression analysis of prognostic factors; and studies regarding wound problems. The references of all included papers were also examined to identify other relevant articles. Any disagreements between the reviewers were resolved through discussion."
}

# Read Excel file
df = pd.read_excel(input_file)

# Replace strings with regex enabled
df.replace(replacements, regex=True, inplace=True)

# Save the updated DataFrame to a new Excel file
df.to_excel(output_file, index=False)
