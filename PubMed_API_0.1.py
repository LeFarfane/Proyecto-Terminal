import requests
import csv
import time
import argparse
import xml.etree.ElementTree as ET
import os

search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
headers = {"User-Agent": "Proyecto-Terminal (eduardo_bio12@outlook.com)"}

api_key = os.getenv("NCBI_API_KEY", "")

parser = argparse.ArgumentParser(description="Download PubMed papers to CSV")
parser.add_argument("--search_term", default="microRNA", help="Term to search for in PubMed")
parser.add_argument("--max_results", type=int, default=5, help="Maximum number of results to fetch")
parser.add_argument("--start_year", type=int, default=None, help="Start year for articles (optional)")
args = parser.parse_args()

# obtener PMIDs a partir del término de búsqueda
pmids = []
try:
    search_params = {
        "db": "pubmed",
        "term": args.search_term,
        "retmode": "json",
        "retmax": args.max_results,
    }
    if api_key:
        search_params["api_key"] = api_key
    if args.start_year:
        search_params["mindate"] = args.start_year
        search_params["datetype"] = "pdat"
    search_response = requests.get(search_url, params=search_params, headers=headers)
    search_response.raise_for_status()
    pmids = search_response.json().get("esearchresult", {}).get("idlist", [])
except Exception as e:
    print(f"Error searching term '{args.search_term}': {e}")

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
        if api_key:
            params["api_key"] = api_key
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
