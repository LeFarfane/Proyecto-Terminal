import argparse
import csv
import json
import logging
import os
import re
import time
import unicodedata
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
import email.utils
import requests
from typing import List

search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
headers = {"User-Agent": "Proyecto-Terminal (eduardo_bio12@outlook.com)"}

TOOL = "Proyecto-Terminal"
EMAIL = "eduardo_bio12@outlook.com"

__version__ = "0.2"

logger = logging.getLogger(__name__)

session = requests.Session()
session.headers.update(headers)

api_key = os.getenv("NCBI_API_KEY", "")


def setup_logging(log_file: str) -> None:
    """Configure logging to console and a file."""
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)


def normalize_text(text: str) -> str:
    """Normalize unicode text and collapse internal whitespace."""
    if not text:
        return ""
    text = unicodedata.normalize("NFKC", text)
    return re.sub(r"\s+", " ", text).strip()


def requests_get_with_retries(
    url, params=None, headers=None, max_retries=5, backoff_factor=1, timeout=10
):
    """Make GET requests with retries on network errors or 429/5xx responses.

    Respects ``Retry-After`` headers when provided by the server to comply with
    usage guidelines. Implements exponential backoff between retries.
    """
    for attempt in range(max_retries):
        try:
            response = session.get(url, params=params, headers=headers, timeout=timeout)
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


def search_pmids(term, max_results, api_key, mindate="2020", maxdate="3000", sort="relevance"):
    """Return a deduplicated list of PMIDs for the search term."""
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
            "tool": TOOL,
            "email": EMAIL,
            "sort": sort,
        }
        if api_key:
            search_params["api_key"] = api_key
        search_response = requests_get_with_retries(
            search_url, params=search_params, headers=headers
        )
        search_response.raise_for_status()
        pmids = search_response.json().get("esearchresult", {}).get("idlist", [])
        pmids = list(dict.fromkeys(pmids))
    except Exception as e:
        logger.error("Error searching term '%s': %s", term, e)
    return pmids


def format_authors_apa(authors_list):
    formatted = []
    for last, initials in authors_list:
        if not last:
            continue
        initials_formatted = ". ".join(list(initials)) + "." if initials else ""
        formatted.append(f"{last}, {initials_formatted}".strip())
    if len(formatted) > 21:
        formatted = formatted[:19] + ["..."] + [formatted[-1]]
    if len(formatted) > 1:
        return ", ".join(formatted[:-1]) + ", & " + formatted[-1]
    return formatted[0] if formatted else ""


def build_apa_citation(authors, year, title, journal, volume, issue, pages, doi):
    citation = f"{authors} ({year}). {title}. {journal}"
    if volume:
        citation += f", {volume}"
    if issue:
        citation += f"({issue})"
    if pages:
        citation += f", {pages}"
    if doi:
        citation += f". https://doi.org/{doi}"
    else:
        citation += "."
    return citation


def write_readme(term: str, count: int, version: str) -> None:
    date = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
    content = f"""# Proyecto-Terminal

## Última búsqueda
- Término: {term}
- Artículos guardados: {count}

Generado el {date} con versión {version}.

## Environment Variables

Set the `NCBI_API_KEY` environment variable to use an NCBI API key for PubMed requests:

```bash
export NCBI_API_KEY=\"your_api_key_here\"
```

This increases the request rate limits when running `PubMed_API_0.1.py`.
"""
    with open("README.md", "w", encoding="utf-8") as f:
        f.write(content)


def download_articles(
    pmids,
    api_key,
    output_csv,
    output_jsonl,
    article_types,
    languages,
    batch_size=100,
):
    """Fetch article details for each PMID and write them to CSV and JSONL."""
    saved = 0
    csv_writer = None
    csvfile = None
    jsonlfile = None
    try:
        if output_csv:
            csv_exists = os.path.isfile(output_csv)
            csvfile = open(output_csv, "a", newline="", encoding="utf-8")
            csv_writer = csv.writer(csvfile)
            if not csv_exists:
                csv_writer.writerow(
                    [
                        "PMID",
                        "Title",
                        "Abstract",
                        "Authors",
                        "Year",
                        "Journal",
                        "Volume",
                        "Issue",
                        "Pages",
                        "DOI",
                        "citation_apa",
                    ]
                )
        if output_jsonl:
            jsonlfile = open(output_jsonl, "a", encoding="utf-8")

        for i in range(0, len(pmids), batch_size):
            batch = pmids[i : i + batch_size]
            params = {
                "db": "pubmed",
                "id": ",".join(batch),
                "retmode": "xml",
                "tool": TOOL,
                "email": EMAIL,
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

                    langs = [normalize_text(l.text).lower() for l in root.findall(".//Language")]
                    if languages and not any(l in languages for l in langs):
                        continue

                    pub_types = [
                        normalize_text(pt.text).lower()
                        for pt in root.findall(".//PublicationTypeList/PublicationType")
                    ]
                    if article_types and not any(pt in article_types for pt in pub_types):
                        continue

                    title_elem = root.find(".//ArticleTitle")
                    title = (
                        normalize_text("".join(title_elem.itertext()))
                        if title_elem is not None
                        else ""
                    )

                    abstract_parts = []
                    for node in root.findall(".//Abstract/AbstractText"):
                        text = normalize_text("".join(node.itertext()))
                        label = node.get("Label") or node.get("NlmCategory")
                        abstract_parts.append(f"{label}: {text}" if label else text)
                    abstract = "\n".join(abstract_parts)

                    authors_elems = root.findall(".//AuthorList/Author")
                    authors_list = []
                    for author in authors_elems:
                        last = normalize_text(author.findtext("LastName", ""))
                        initials = normalize_text(author.findtext("Initials", ""))
                        if last or initials:
                            authors_list.append((last, initials))
                    authors = "; ".join(
                        [f"{last} {initials}".strip() for last, initials in authors_list]
                    )

                    year = normalize_text(root.findtext(".//PubDate/Year", default=""))
                    journal = normalize_text(root.findtext(".//Journal/Title", default=""))
                    volume = normalize_text(
                        root.findtext(".//JournalIssue/Volume", default="")
                    )
                    issue = normalize_text(
                        root.findtext(".//JournalIssue/Issue", default="")
                    )
                    pages = normalize_text(root.findtext(".//MedlinePgn", default=""))
                    doi = normalize_text(
                        root.findtext(
                            ".//ArticleIdList/ArticleId[@IdType='doi']",
                            default="",
                        )
                    )
                    pmid = normalize_text(root.findtext(".//PMID", default=""))

                    citation = build_apa_citation(
                        format_authors_apa(authors_list),
                        year,
                        title,
                        journal,
                        volume,
                        issue,
                        pages,
                        doi,
                    )

                    row = [
                        pmid,
                        title,
                        abstract,
                        authors,
                        year,
                        journal,
                        volume,
                        issue,
                        pages,
                        doi,
                        citation,
                    ]
                    if csv_writer:
                        csv_writer.writerow(row)
                    if jsonlfile:
                        jsonlfile.write(
                            json.dumps(
                                {
                                    "PMID": pmid,
                                    "Title": title,
                                    "Abstract": abstract,
                                    "Authors": authors,
                                    "Year": year,
                                    "Journal": journal,
                                    "Volume": volume,
                                    "Issue": issue,
                                    "Pages": pages,
                                    "DOI": doi,
                                    "citation_apa": citation,
                                },
                                ensure_ascii=False,
                            )
                            + "\n"
                        )
                    saved += 1
            except Exception as e:
                logger.error("Error retrieving PMIDs %s: %s", batch, e)
            time.sleep(0.1)
    finally:
        if csvfile:
            csvfile.close()
        if jsonlfile:
            jsonlfile.close()
    return saved


def run_query(
    terms: List[str],
    operator: str,
    sort: str,
    pub_types: List[str],
    batch_size: int,
    max_results: int,
    api_key: str,
    format: str,
    append: bool,
):
    query = f" {operator} ".join(terms)
    setup_logging("PubMed_API_0.1.log")
    logger.info("Script version: %s", __version__)
    logger.info(
        "Run date: %s", datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
    )
    logger.info(
        "Parameters: %s",
        {
            "terms": terms,
            "operator": operator,
            "sort": sort,
            "pub_types": pub_types,
            "batch_size": batch_size,
            "max_results": max_results,
            "api_key_provided": bool(api_key),
            "format": format,
            "append": append,
        },
    )

    pmids = search_pmids(query, max_results, api_key, sort=sort)
    article_types = [normalize_text(t).lower() for t in pub_types]
    languages = []

    output_csv = "papers.csv" if format in ("csv", "both") else None
    output_jsonl = "papers.jsonl" if format in ("jsonl", "both") else None

    if not append:
        if output_csv and os.path.exists(output_csv):
            os.remove(output_csv)
        if output_jsonl and os.path.exists(output_jsonl):
            os.remove(output_jsonl)

    count = download_articles(
        pmids,
        api_key,
        output_csv,
        output_jsonl,
        article_types,
        languages,
        batch_size,
    )
    write_readme(query, count, __version__)
    logger.info("Saved %d articles", count)
    return count


def main(**kwargs):
    return run_query(**kwargs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download PubMed papers to CSV and JSONL"
    )
    parser.add_argument(
        "--terms",
        default="microRNA",
        help="Semicolon-separated terms to search for in PubMed",
    )
    parser.add_argument(
        "--operator",
        default="AND",
        choices=["AND", "OR", "NOT"],
        help="Operator to combine terms",
    )
    parser.add_argument(
        "--sort", default="relevance", help="Sort order for results"
    )
    parser.add_argument(
        "--pub_types",
        default="Journal Article",
        help="Comma-separated publication types to include",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=100,
        help="efetch batch size (50-200)",
    )
    parser.add_argument(
        "--max_results",
        type=int,
        default=5,
        help="Maximum number of results to fetch",
    )
    parser.add_argument(
        "--api_key", default="", help="NCBI API key"
    )
    parser.add_argument(
        "--format",
        default="both",
        choices=["csv", "jsonl", "both"],
        help="Output format",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to existing files instead of overwriting",
    )
    args = parser.parse_args()
    terms = [t.strip() for t in args.terms.split(";") if t.strip()]
    pub_types = [pt.strip() for pt in args.pub_types.split(",") if pt.strip()]
    main(
        terms=terms,
        operator=args.operator,
        sort=args.sort,
        pub_types=pub_types,
        batch_size=args.batch_size,
        max_results=args.max_results,
        api_key=args.api_key,
        format=args.format,
        append=args.append,
    )
