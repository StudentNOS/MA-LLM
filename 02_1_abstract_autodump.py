import argparse
import csv
import xml.etree.ElementTree as ET
import urllib
from urllib.request import urlopen


FIELDS = ["pubmedid", "link", "year", "drug", "journal",
          "title", "objective", "background", "methods",
          "results", "conclusions", "abstract"]

# terms: search terms
def getids(searchterms):
    url_str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="+searchterms
    url_encod = urllib.parse.quote(url_str, safe=':/?&=')
    req = urlopen(url_encod)
    print("IDS")
    print(type(req))
    root = ET.fromstring(req.read())
    ids = next(root.iter("IdList"))
    id = ids.findall("Id")
    idlist = []
    for i in id:
        idlist.append(i.text)
    return idlist

def handle_journal(root, result):
    journal = next(root.iter("Journal"))
    result["journal"] = journal.find("Title").text
    date = next(journal.iter("PubDate"))
    result["year"] = date.find("Year").text


def handle_title(root, result):
    result["title"] = next(root.iter("ArticleTitle")).text


ABS_SWITCH = {
    "OBJECTIVE": "objective",
    "METHODS": "methods",
    "RESULTS": "results",
    "CONCLUSIONS": "conclusions",
    "BACKGROUND": "background",
#    None: "results",
}


def handle_abstract(root, result):
    pieces = []
    for txt in root.iter("AbstractText"):
        cat = txt.get("NlmCategory")
        if cat in ABS_SWITCH:
            result[ABS_SWITCH[cat]] = txt.text
        else:
            pieces.append(txt.text)
    result["abstract"] = ". ".join(pieces)

def handle_drug(root, result):
    drugs = []
    for mesh in root.iter("MeshHeading"):
        for qn in mesh.findall("QualifierName"):
            if qn.text == "pharmacology":
                drugs.append(mesh.find("DescriptorName").text)
                break  # just out of MeshHeading
    result["drug"] = ", ".join(drugs)

def dump_info(pmid):
    result = {"pubmedid": pmid}
    resp = urlopen(DUMP_URL % pmid)
    result["link"] = PAPER_URL % pmid
    root = ET.fromstring(resp.read())
    handle_journal(root, result)
    handle_title(root, result)
    handle_abstract(root, result)
    handle_drug(root, result)
    for f in FIELDS:
        if f in result.keys():
            result[f] = result[f].encode("utf-8")
    
    return result


search = '"Mineralocorticoid Receptor Antagonists"[Mesh] OR "MRA"[tw] OR "aldosterone receptor blocker"[tw] OR "aldosterone blocker"[tw] OR "spironolactone"[tw] OR "eplerenone"[tw] OR "canrenone"[tw] OR "finerenone"[tw] AND ""Randomized Controlled Trial" [Publication Type]"[Mesh] OR "randomized trial"[tw] OR "randomised trial"[tw] AND "Heart Failure"[Mesh] ) OR ("Sodium-Glucose Transporter 2 Inhibitors"[Mesh] OR "SGLT-2 inhibitors"[tw] OR "SGLT2i"[tw] OR "canagliflozin"[tw] OR "empagliflozin"[tw] OR "dapagliflozin"[tw] OR "ertugliflozin"[tw] OR "sotagliflozin"[tw] OR "remogliflozin"[tw] OR "luseugliflozin"[tw] OR "ipragliflozin"[tw] AND ""Randomized Controlled Trial" [Publication Type]"[Mesh] OR "randomized trial"[tw] OR "randomised trial"[tw] AND "Heart Failure"[Mesh]'
idlist = getids(search)

print ("Length of ID list")
print (len(idlist))

DUMP_URL =  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id=%s"
PAPER_URL = "https://www.ncbi.nlm.nih.gov/pubmed/%s"

with open("ncbi_autodumped.csv", "w") as fout:
    writer = csv.DictWriter(fout, fieldnames=FIELDS)
    writer.writeheader()
    for pmid in idlist:
        print(pmid)
        writer.writerow(dump_info(pmid.strip()))

