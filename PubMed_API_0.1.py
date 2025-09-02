import requests
import csv
import time
import xml.etree.ElementTree as ET

# lista de PubMed artículos
pmids = ["16957370", "25471818", "34370220", "27807838", "21342132"]

# creación del documento:
with open("papers.csv", "w", newline='', encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["PMID", "Title", "Abstract"])

    for pmid in pmids:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "pubmed",
            "id": pmid,
            "retmode": "xml",
        }
        try:
            response = requests.get(url, params=params, headers={
                "User-Agent": "Proyecto-Terminal (your_email@example.com)"
            })
            response.raise_for_status()
            root = ET.fromstring(response.text)
            title = root.findtext(".//ArticleTitle", default="")
            abstract = root.findtext(".//Abstract/AbstractText", default="")
            writer.writerow([pmid, title, abstract])
        except Exception as e:
            print(f"Error retrieving PMID {pmid}: {e}")
        time.sleep(0.34)
