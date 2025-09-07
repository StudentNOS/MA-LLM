# ENSURE: Automated Systematic Review Screening Pipeline

## Table of Contents
1. [What is ENSURE?](#what-is-ensure)
2. [Who is this for?](#who-is-this-for)
3. [Quick Start Guide](#quick-start-guide)
4. [Detailed Setup Instructions](#detailed-setup-instructions)
5. [How to Use the Pipeline](#how-to-use-the-pipeline)
6. [Understanding Your Data](#understanding-your-data)
7. [Interpreting Results](#interpreting-results)
8. [Troubleshooting](#troubleshooting)
9. [Best Practices](#best-practices)
10. [Getting Help](#getting-help)

---

## What is ENSURE?

ENSURE is an automated screening tool designed specifically for medical researchers conducting systematic reviews and meta-analyses. Instead of manually reading through hundreds or thousands of research papers to determine their relevance, ENSURE uses artificial intelligence to automatically screen papers based on your specific criteria.

### What does ENSURE do?

- **Automatically downloads paper information** from PubMed using paper IDs
- **Screens paper titles and abstracts** using AI based on your research criteria
- **Calculates performance metrics** to show how well the AI performed compared to expert human reviewers
- **Provides detailed reports** in Excel format for further analysis

### Key Benefits

✅ **Save time**: Screen thousands of papers in minutes instead of weeks  
✅ **Reduce human error**: Consistent screening criteria applied to all papers  
✅ **Evidence-based**: Compare AI performance against expert human reviewers  
✅ **User-friendly**: Web interface - no programming knowledge required  
✅ **Flexible**: Works with any systematic review topic

---

## Who is this for?

This tool is designed for:

- **Medical researchers** conducting systematic reviews
- **Meta-analysis teams** screening large numbers of studies
- **Research groups** wanting to standardize their screening process
- **Anyone** who needs to efficiently review large sets of medical literature

### What you need to know

**You DO NOT need:**
- Programming experience
- Command line knowledge
- Database management skills

**You DO need:**
- Access to a computer
- An internet connection
- Basic Excel skills
- Your research question and inclusion/exclusion criteria

---

## Quick Start Guide

### Prerequisites (Get these first)

1. **Groq API Key**: Free AI service for paper screening
   - Go to [https://console.groq.com/](https://console.groq.com/)
   - Create a free account
   - Get your API key (it looks like: `gsk_...`)

2. **Email address**: For PubMed access (use your institutional email)

3. **Your data files** (we'll help you create these):
   - List of PubMed IDs to screen
   - Your screening prompts
   - Gold standard comparison set (optional but recommended)

### 30-Second Overview

1. **Start the tool** by running one file
2. **Upload your files** through a web interface
3. **Enter your API credentials**
4. **Click "Start Screening"**
5. **Download your results** as an Excel file

That's it! The AI does the screening work for you.

---

## Detailed Setup Instructions

### Step 1: Install Required Software

#### Option A: Simple Installation (Recommended)
1. **Download Python**:
   - Go to [https://python.org](https://python.org)
   - Download Python 3.8 or newer
   - During installation, check "Add Python to PATH"

2. **Download the ENSURE files**:
   - Download all files from the Pipeline folder
   - Save them in a new folder on your computer (e.g., `ENSURE_Pipeline`)

3. **Install dependencies**:
   - Open Command Prompt (Windows) or Terminal (Mac/Linux)
   - Navigate to your ENSURE_Pipeline folder
   - Type: `pip install pandas biopython groq flask openpyxl`
   - Press Enter and wait for installation to complete

#### Option B: Advanced Installation (For IT departments)
```bash
# Clone the repository
git clone https://github.com/your-repository/ensure.git
cd ensure/Pipeline

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Step 2: Get Your API Credentials

#### Groq API Key (Free)
1. Visit [https://console.groq.com/](https://console.groq.com/)
2. Create an account with your email
3. Go to "API Keys" section
4. Click "Create API Key"
5. Copy the key (starts with `gsk_`) - you'll need this later
6. **Important**: Keep this key private and secure

#### Email for PubMed
- Use your institutional email address
- This is required by PubMed for API access
- No special registration needed

### Step 3: Prepare Your Data Files

You need to create three files before starting:

#### File 1: Initial PMIDs (Required)
**What it is**: List of PubMed IDs you want to screen  
**File name**: `Initial.txt`  
**Format**: Plain text file, one PMID per line

**Example**:
```
32420649
22901346
34025474
25687662
```

**How to get PMIDs**:
1. Search PubMed for your topic
2. Select all relevant papers from your search
3. Go to "Send to" → "File" → "PMID List"
4. Save as text file and rename to `Initial.txt`

#### File 2: Screening Prompts (Required)
**What it is**: Your inclusion/exclusion criteria as text prompts  
**File name**: `Prompts.xlsx`  
**Format**: Excel file with specific columns

**Required columns**:
- `screen_titles`: 1 (screen titles) or 0 (skip titles)
- `screen_abstracts`: 1 (screen abstracts) or 0 (skip abstracts)
- `TitlePrompt`: Text describing what titles to include
- `AbstractPrompt`: Text describing what abstracts to include

**Example Prompts.xlsx**:
| screen_titles | screen_abstracts | TitlePrompt | AbstractPrompt |
|---------------|------------------|-------------|----------------|
| 1 | 1 | "Select titles about randomized controlled trials for diabetes treatment in adults" | "Include abstracts describing RCTs testing diabetes medications in patients over 18 years old" |

**Writing good prompts**:
- Be specific about your population (e.g., "adults over 18")
- Include study types (e.g., "randomized controlled trials")
- Mention key interventions or outcomes
- Use clear, simple language
- Test with a few examples first

#### File 3: Gold Standard (Optional but Recommended)
**What it is**: PMIDs that expert reviewers determined should be included  
**File name**: `Goldstandard_Selected.txt`  
**Format**: Plain text file, one PMID per line
**Note**: Must be a subset of your Initial.txt PMIDs

**Example**:
```
32420649
34025474
```

**Why use a gold standard**:
- Measure how well the AI performed compared to human experts
- Get statistical metrics (sensitivity, specificity, etc.)
- Validate your prompts before large-scale screening

---

## How to Use the Pipeline

### Step 1: Start the Application

1. **Open Command Prompt/Terminal**
2. **Navigate to your ENSURE_Pipeline folder**:
   ```
   cd path/to/your/ENSURE_Pipeline
   ```
3. **Start the application**:
   ```
   python screening.py
   ```
4. **Open your web browser** and go to: `http://localhost:5000`

You should see the ENSURE web interface.

### Step 2: Configure the Screening

In the web interface, you'll see a form with several sections:

#### 1. Konfiguration (Configuration)
- **Entrez Email**: Enter your institutional email address
- **Groq API Key**: Paste your API key from Step 2

#### 2. Datendateien (Data Files)
- **Initial PMIDs (.txt)**: Upload your `Initial.txt` file
- **Goldstandard PMIDs (.txt)**: Upload your `Goldstandard_Selected.txt` (optional)
- **Prompts (.xlsx)**: Upload your `Prompts.xlsx` file

#### 3. Prozess starten (Start Process)
Choose one option:
- **Start Screening with Comparison**: Use this if you have a gold standard file
- **Start Screening without Comparison**: Use this for just screening (no performance metrics)

### Step 3: Monitor Progress

After clicking start:
1. **Status updates** will appear showing progress
2. **Don't close your browser** during screening
3. **Processing time** depends on number of papers and prompts (usually 5-10 minutes per 1000 papers)
4. You can **stop the process** anytime using the "Stop Process" button

### Step 4: Download Results

When screening is complete:
1. A **"Download Results"** button will appear
2. Click it to download your Excel results file
3. The file contains all screening results and performance metrics

---

## Understanding Your Data

### Input File Requirements

#### Initial.txt Format
```
✅ Correct:
32420649
22901346
34025474

❌ Incorrect:
PMID: 32420649
32420649, 22901346
```

#### Prompts.xlsx Structure
| Column | Required | Values | Purpose |
|--------|----------|--------|---------|
| screen_titles | Yes | 1 or 0 | Enable/disable title screening |
| screen_abstracts | Yes | 1 or 0 | Enable/disable abstract screening |
| TitlePrompt | If screen_titles=1 | Text | Instructions for title screening |
| AbstractPrompt | If screen_abstracts=1 | Text | Instructions for abstract screening |

### Data Flow Process

1. **Upload PMIDs** → System downloads paper details from PubMed
2. **Title Screening** → AI reads titles and applies your criteria
3. **Abstract Screening** → AI reads abstracts of selected papers
4. **Performance Calculation** → Compare AI results to gold standard
5. **Generate Report** → Create Excel file with all results

### What the AI Sees

The AI receives your papers in this format:
```
PMID: 32420649
Title: Systematic review of diabetes interventions
Abstract: Background: This study examines...
```

And your prompt:
```
"Select papers about randomized controlled trials for diabetes treatment in adults"
```

The AI then decides: Include ✅ or Exclude ❌

---

## Interpreting Results

### Understanding the Excel Output

Your results file contains several important columns:

#### Input Columns (Your original data)
- `screen_titles`, `screen_abstracts`: Your screening settings
- `TitlePrompt`, `AbstractPrompt`: Your original prompts

#### Results Columns
- `Final_Relevant_PMIDs`: Comma-separated list of PMIDs the AI selected
- Count: Number of papers in your final inclusion set

#### Performance Metrics (if you used gold standard)

**Key Metrics Explained**:

| Metric | What it means | Good value |
|--------|---------------|------------|
| **Sensitivity** | % of relevant papers correctly identified | >90% |
| **Specificity** | % of irrelevant papers correctly excluded | >90% |
| **PPV (Precision)** | % of selected papers that are actually relevant | >80% |
| **NPV** | % of excluded papers that are actually irrelevant | >95% |

**Confusion Matrix Values**:
- `tp` (True Positives): Relevant papers correctly included
- `tn` (True Negatives): Irrelevant papers correctly excluded  
- `fp` (False Positives): Irrelevant papers incorrectly included
- `fn` (False Negatives): Relevant papers incorrectly excluded

### Reading Performance Results

#### Example Result Interpretation:
```
sensitivity_titles: 0.95 (95%)
specificity_titles: 0.88 (88%)
PPV_titles: 0.82 (82%)
tp_titles: 19
fp_titles: 4
fn_titles: 1
tn_titles: 76
```

**This means**:
- ✅ The AI found 95% of relevant papers (high sensitivity)
- ✅ The AI correctly excluded 88% of irrelevant papers
- ⚠️ Of papers the AI selected, 82% were actually relevant
- ❌ The AI missed 1 relevant paper
- ❌ The AI incorrectly included 4 irrelevant papers

### What Makes Good Performance?

#### Excellent Performance
- Sensitivity >95%
- Specificity >90%
- PPV >85%
- Very few false negatives (missed relevant papers)

#### Acceptable Performance  
- Sensitivity >90%
- Specificity >80%
- PPV >70%
- Manual review of false positives is manageable

#### Poor Performance
- Sensitivity <85%
- Many false negatives (missing important papers)
- High false positive rate requiring extensive manual review

### Improving Poor Performance

If your results are unsatisfactory:

1. **Refine your prompts**:
   - Be more specific about inclusion criteria
   - Add examples of what to include/exclude
   - Use medical terminology consistently

2. **Adjust screening strategy**:
   - Screen titles only for broad screening
   - Use abstracts for detailed screening
   - Consider two-stage screening

3. **Test with smaller datasets**:
   - Use 100-200 papers to test prompts
   - Iterate quickly with different prompt versions
   - Scale up once performance is acceptable

---

## Troubleshooting

### Common Issues and Solutions

#### 1. "API Key Error" / "Authentication Failed"
**Problem**: Invalid or missing Groq API key  
**Solutions**:
- Check your API key is copied correctly (no extra spaces)
- Verify your Groq account is active
- Generate a new API key if needed
- Check you have API credits remaining

#### 2. "PubMed Connection Error"
**Problem**: Cannot download paper information  
**Solutions**:
- Check your internet connection
- Verify your email address is correct
- Wait a few minutes (PubMed may have temporary limits)
- Try with a smaller set of PMIDs first

#### 3. "File Format Error"
**Problem**: Uploaded files are not in correct format  
**Solutions**:
- Ensure Initial.txt has one PMID per line (no headers)
- Check PMIDs are 7-9 digits only
- Verify Excel file has required columns
- Save Excel file as .xlsx format (not .xls)

#### 4. "No Results Generated"
**Problem**: Process completes but no output file  
**Solutions**:
- Check for error messages in the status
- Verify you have write permissions in the folder
- Try running with a smaller dataset
- Restart the application

#### 5. "Process Stuck/Not Responding"
**Problem**: Screening appears to freeze  
**Solutions**:
- Wait (large datasets take time)
- Check your internet connection
- Monitor status updates
- Use "Stop Process" and restart if necessary

#### 6. "Poor Performance Results"
**Problem**: Low sensitivity/specificity scores  
**Solutions**:
- Review and refine your prompts
- Check your gold standard is accurate
- Test with a subset of papers first
- Consider adjusting screening strategy

### Getting More Help

#### Log Files and Error Messages
- The application shows status messages during processing
- Note any error messages exactly as shown
- Check the terminal/command prompt for additional error details

#### System Requirements
- **Minimum**: 4GB RAM, stable internet connection
- **Recommended**: 8GB RAM for large datasets (>1000 papers)
- **Storage**: 1GB free space for results and temporary files

#### Performance Optimization
- **For large datasets**: Process in smaller batches
- **For slow internet**: Reduce batch sizes
- **For limited memory**: Close other applications during processing

---

## Best Practices

### Before You Start

#### 1. Plan Your Screening Strategy
- **Define clear inclusion/exclusion criteria** before writing prompts
- **Decide on screening stages**: titles only, abstracts only, or both
- **Prepare your gold standard** with at least 50-100 manually reviewed papers
- **Test with small datasets** before processing thousands of papers

#### 2. Write Effective Prompts
**Do**:
- Be specific about study population (age, condition, etc.)
- Include study types you want (RCT, cohort study, etc.)
- Mention key interventions or comparisons
- Use standard medical terminology
- Keep prompts focused on one concept at a time

**Don't**:
- Use vague terms like "relevant" or "important"
- Include multiple complex criteria in one prompt
- Use abbreviations without explanation
- Write overly long prompts (>200 words)

#### 3. Validate Your Approach
- **Start small**: Test with 100-200 papers
- **Compare results**: Use gold standard comparison
- **Iterate prompts**: Refine based on performance metrics
- **Document decisions**: Keep track of what works

### During Screening

#### 1. Monitor Progress
- **Watch status updates** to ensure processing continues
- **Note processing speed** to estimate completion time
- **Don't interrupt** unless absolutely necessary
- **Check for errors** in status messages

#### 2. Quality Control
- **Save intermediate results** if processing large datasets
- **Backup your input files** before starting
- **Monitor system resources** (memory, internet connection)

### After Screening

#### 1. Review Results Carefully
- **Check performance metrics** against your expectations
- **Review false positives and negatives** to understand AI decisions
- **Validate unexpected results** with manual review
- **Document any issues** for future screening rounds

#### 2. Manual Review Strategy
- **Always review false negatives** (missed relevant papers)
- **Sample false positives** to understand over-inclusion patterns
- **Consider hybrid approach**: AI pre-screening + manual final review
- **Update exclusion criteria** based on AI performance

#### 3. Reporting Your Methods
When publishing your systematic review:
- **Describe your prompts** exactly as used
- **Report performance metrics** (sensitivity, specificity)
- **Include AI screening as part of your methods**
- **Acknowledge limitations** of automated screening

### Scaling Up

#### For Large Projects (>5000 papers)
- **Process in batches** of 1000-2000 papers
- **Use high-performance computers** if available
- **Monitor API usage** and costs
- **Plan for longer processing times**

#### For Multiple Reviews
- **Standardize prompt templates** across projects
- **Maintain a prompt library** of effective prompts
- **Document lessons learned** for future use
- **Consider training team members** on the tool

---

## Getting Help

### Quick Reference

#### File Formats Checklist
- ✅ Initial.txt: One PMID per line, no headers
- ✅ Goldstandard_Selected.txt: Subset of Initial.txt PMIDs
- ✅ Prompts.xlsx: Required columns present
- ✅ All files saved in correct formats (.txt, .xlsx)

#### Common Commands
```bash
# Start the pipeline
python screening.py

# If you get "python not found"
python3 screening.py

# Check if Python is installed
python --version
```

#### API Limits and Costs
- **Groq**: Currently free with usage limits
- **PubMed**: Free but has rate limits
- **Processing time**: ~5-10 minutes per 1000 papers

### Support Resources

#### Documentation
- **This guide**: Complete usage instructions
- **API Documentation**: Technical details for developers
- **Setup Guide**: Detailed installation instructions

#### Community Support
- **GitHub Issues**: Report bugs or request features
- **Discussion Forums**: Ask questions and share experiences
- **Email Support**: Contact the development team

#### Training and Workshops
- Request training sessions for your research team
- Webinars on systematic review automation
- Best practices workshops

### Frequently Asked Questions

#### General Questions

**Q: Do I need programming experience?**  
A: No, the web interface is designed for non-programmers.

**Q: How accurate is the AI screening?**  
A: Typically 85-95% sensitivity with well-written prompts. Always validate with your gold standard.

**Q: Can I use this for non-medical systematic reviews?**  
A: Yes, but it's optimized for PubMed medical literature.

**Q: How long does screening take?**  
A: Usually 5-10 minutes per 1000 papers, depending on internet speed and AI processing time.

#### Technical Questions

**Q: What AI model does ENSURE use?**  
A: Currently uses Groq's Llama models, which are fast and effective for text screening.

**Q: Can I use my own AI service?**  
A: The code can be modified to use other AI services, but requires programming knowledge.

**Q: Is my data secure?**  
A: PubMed data is public. Your prompts are sent to Groq's API. Check your institution's data policies.

**Q: Can I run this offline?**  
A: No, it requires internet access for PubMed and AI services.

#### Usage Questions

**Q: How many papers can I screen at once?**  
A: Tested with up to 10,000 papers. For larger datasets, consider batch processing.

**Q: Can I modify my prompts during screening?**  
A: No, you need to stop and restart with new prompts.

**Q: What if some PMIDs aren't found in PubMed?**  
A: The system will skip invalid PMIDs and continue with valid ones.

**Q: Can I screen non-English papers?**  
A: Yes, but AI performance may vary. English prompts work best.

---

## Conclusion

ENSURE provides an efficient, user-friendly way to automate systematic review screening without requiring programming expertise. By following this guide, medical researchers can:

- Set up the pipeline quickly and easily
- Screen large numbers of papers automatically  
- Validate AI performance against expert reviewers
- Generate detailed reports for analysis

Remember: AI screening is a tool to assist, not replace, human expertise in systematic reviews. Always validate results and use clinical judgment in your final decisions.

For additional support, feature requests, or to report issues, please contact the ENSURE development team or visit our GitHub repository.