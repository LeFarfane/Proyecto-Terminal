import requests
import csv

# lista de PubMed artículos
pmids = ["16957370", "25471818", "34370220", "27807838", "21342132"]

#creación del documento:
with open("papers.csv", "w", newline='', encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["PMID", "Title", "Abstract"])

    for pmid in pmids:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "pubmed",
            "id": pmid,
            "retmode": "xml",
            "retstart": 0,
            "retmax": 100,
            "retsort": "relevance",
            "term": "PubMed",
        }
        response = requests.get(url, params=params)
        xml = response.text
        root = ET.fromstring(xml)

        writer.writerow([pmid, root.find("title").text, root.find("abstract").text])
