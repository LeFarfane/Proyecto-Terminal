import requests
import csv
import time
import xml.etree.ElementTree as ET

# término de búsqueda para PubMed
search_term = "microRNA"

search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
headers = {"User-Agent": "Proyecto-Terminal (eduardo_bio12@outlook.com)"}

# obtener PMIDs a partir del término de búsqueda
pmids = []
try:
    search_params = {
        "db": "pubmed",
        "term": search_term,
        "retmode": "json",
        "retmax": 5,
    }
    search_response = requests.get(search_url, params=search_params, headers=headers)
    search_response.raise_for_status()
    pmids = search_response.json().get("esearchresult", {}).get("idlist", [])
except Exception as e:
    print(f"Error searching term '{search_term}': {e}")

# creación del documento:
with open("papers.csv", "w", newline='', encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["PMID", "Title", "Abstract"])

    for pmid in pmids:
        params = {
            "db": "pubmed",
            "id": pmid,
            "retmode": "xml",
        }
        try:
            response = requests.get(fetch_url, params=params, headers=headers)
            response.raise_for_status()
            root = ET.fromstring(response.text)
            title = root.findtext(".//ArticleTitle", default="")
            abstract = root.findtext(".//Abstract/AbstractText", default="")
            writer.writerow([pmid, title, abstract])
        except Exception as e:
            print(f"Error retrieving PMID {pmid}: {e}")
        time.sleep(0.1)
