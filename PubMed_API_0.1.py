import requests
import csv
import time
import argparse
import xml.etree.ElementTree as ET
import os
from datetime import datetime, timezone
import email.utils

search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
headers = {"User-Agent": "Proyecto-Terminal (eduardo_bio12@outlook.com)"}

api_key = os.getenv("NCBI_API_KEY", "")


def requests_get_with_retries(url, params=None, headers=None, max_retries=5, backoff_factor=1):
    """Make GET requests with retries on network errors or 429/5xx responses.

    Respects ``Retry-After`` headers when provided by the server to comply with
    usage guidelines. Implements exponential backoff between retries.
    """
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, headers=headers)
            # Retry on 429 or server errors
            if response.status_code == 429 or response.status_code >= 500:
                retry_after = response.headers.get("Retry-After")
                if retry_after:
                    try:
                        delay = float(retry_after)
                    except ValueError:
                        try:
                            retry_date = email.utils.parsedate_to_datetime(retry_after)
                            delay = (retry_date - datetime.now(timezone.utc)).total_seconds()
                        except Exception:
                            delay = backoff_factor * (2 ** attempt)
                else:
                    delay = backoff_factor * (2 ** attempt)
                if delay > 0:
                    time.sleep(delay)
                continue
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException:
            if attempt == max_retries - 1:
                raise
            time.sleep(backoff_factor * (2 ** attempt))
    raise RuntimeError("Max retries exceeded for GET request")


def search_pmids(term, max_results, api_key, mindate="2020", maxdate="3000"):
    """Return a list of PMIDs for the search term."""
    pmids = []
    try:
        search_params = {
            "db": "pubmed",
            "term": term,
            "retmode": "json",
            "retmax": max_results,
            "datetype": "pdat",
            "mindate": mindate,
            "maxdate": maxdate,
        }
        if api_key:
            search_params["api_key"] = api_key
        search_response = requests_get_with_retries(
            search_url, params=search_params, headers=headers
        )
        search_response.raise_for_status()
        pmids = search_response.json().get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        print(f"Error searching term '{term}': {e}")
    return pmids


def download_articles(pmids, api_key, output_path):
    """Fetch article details for each PMID and write them to a CSV file."""
    batch_size = 100  # NCBI recommends batching 50-100 IDs per request
    with open(output_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PMID", "Title", "Abstract", "Authors", "Year", "Journal"])

        for i in range(0, len(pmids), batch_size):
            batch = pmids[i : i + batch_size]
            params = {
                "db": "pubmed",
                "id": ",".join(batch),
                "retmode": "xml",
            }
            if api_key:
                params["api_key"] = api_key
            try:
                response = requests_get_with_retries(
                    fetch_url, params=params, headers=headers
                )
                response.raise_for_status()
                batch_root = ET.fromstring(response.text)
                for article in batch_root.findall("PubmedArticle"):
                    root = article
                    title_elem = root.find(".//ArticleTitle")
                    title = (
                        "".join(title_elem.itertext()).strip()
                        if title_elem is not None
                        else ""
                    )

                    abstract_parts = []
                    for node in root.findall(".//Abstract/AbstractText"):
                        text = "".join(node.itertext()).strip()
                        label = node.get("Label")
                        if label:
                            abstract_parts.append(f"{label}: {text}")
                        else:
                            abstract_parts.append(text)
                    abstract = " ".join(abstract_parts)

                    authors = "; ".join(
                        [
                            f"{author.findtext('LastName', '')} {author.findtext('Initials', '')}".strip()
                            for author in root.findall(".//AuthorList/Author")
                            if author.findtext("LastName") or author.findtext("Initials")
                        ]
                    )
                    year = root.findtext(".//PubDate/Year", default="")
                    journal = root.findtext(".//Journal/Title", default="")
                    pmid = root.findtext(".//PMID", default="")
                    writer.writerow([pmid, title, abstract, authors, year, journal])
            except Exception as e:
                print(f"Error retrieving PMIDs {batch}: {e}")
            time.sleep(0.1)


def main():
    parser = argparse.ArgumentParser(description="Download PubMed papers to CSV")
    parser.add_argument(
        "--search_term", default="microRNA", help="Term to search for in PubMed"
    )
    parser.add_argument(
        "--max_results",
        type=int,
        default=5,
        help="Maximum number of results to fetch",
    )
    parser.add_argument(
        "--mindate", default="2020", help="Minimum publication year"
    )
    parser.add_argument(
        "--maxdate", default="3000", help="Maximum publication year"
    )
    args = parser.parse_args()

    pmids = search_pmids(
        args.search_term, args.max_results, api_key, args.mindate, args.maxdate
    )
    download_articles(pmids, api_key, "papers.csv")


if __name__ == "__main__":
    main()
