import argparse
import csv
import json
import logging
import os
import re
import time
import unicodedata
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from typing import List

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# URLs base de la API de PubMed y encabezado de identificación
search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
headers = {"User-Agent": "Proyecto-Terminal (eduardo_bio12@outlook.com)"}

# Identificadores requeridos por NCBI
TOOL = "Proyecto-Terminal"
EMAIL = "eduardo_bio12@outlook.com"

__version__ = "0.2"

# Logger principal del módulo
logger = logging.getLogger(__name__)

# Sesión HTTP con estrategia de reintentos
session = requests.Session()
session.headers.update(headers)
retry_strategy = Retry(
    total=5,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["HEAD", "GET", "OPTIONS"],
    backoff_factor=1,
)
adapter = HTTPAdapter(max_retries=retry_strategy)
session.mount("https://", adapter)
session.mount("http://", adapter)

# API key opcional para aumentar el límite de peticiones
api_key = os.getenv("NCBI_API_KEY", "")


def setup_logging(log_file: str) -> None:
    """Configure logging to console and a file."""
    # Prepara el logger para consola y archivo
    logger.setLevel(logging.INFO)
    if logger.handlers:
        logger.handlers.clear()
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
    # Normaliza texto Unicode y elimina espacios extra
    if not text:
        return ""
    text = unicodedata.normalize("NFKC", text)
    return re.sub(r"\s+", " ", text).strip()


def sanitize_term(term: str) -> str:
    """Wrap the term in double quotes for PubMed phrase search."""
    return f'"{term}"'

def search_pmids(term, max_results, api_key, mindate="2020", maxdate="3000", sort="relevance"):
    """Return a deduplicated list of PMIDs for the search term.

    PubMed's ESearch no longer honors ``retmax``. This function paginates
    through results using ``retstart`` until ``max_results`` PMIDs are
    collected or the record count is exhausted.
    """
    pmids = []
    retstart = 0
    count = max_results
    try:
        while len(pmids) < max_results and retstart < count:
            search_params = {
                "db": "pubmed",
                "term": term,
                "retmode": "json",
                "retstart": retstart,
                "datetype": "pdat",
                "mindate": mindate,
                "maxdate": maxdate,
                "tool": TOOL,
                "email": EMAIL,
                "sort": sort,
            }
            if api_key:
                search_params["api_key"] = api_key
            search_response = session.get(search_url, params=search_params, timeout=10)
            search_response.raise_for_status()
            result = search_response.json().get("esearchresult", {})
            ids = result.get("idlist", [])
            count = int(result.get("count", 0))
            pmids.extend(ids)
            retstart += len(ids)
            if not ids:
                break
    except Exception as e:
        logger.error("Error searching term '%s': %s", term, e)
    pmids = list(dict.fromkeys(pmids))
    return pmids[:max_results]


def format_authors_apa(authors_list):
    """Convierte lista de autores al formato APA."""
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
    """Construye la cita en formato APA."""
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
    # Actualiza el README con detalles de la última búsqueda
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


@dataclass
class Article:
    """Modelo de datos para un artículo de PubMed."""
    PMID: str
    Title: str
    Abstract: str
    Authors: str
    Year: str
    Journal: str
    Volume: str
    Issue: str
    Pages: str
    DOI: str
    citation_apa: str
    PubMedURL: str
    PublicationTypes: List[str]
    QueryTerms: List[str]


def download_articles(
    pmids,
    api_key,
    output_csv,
    output_jsonl,
    article_types,
    languages,
    query_terms,
    batch_size=100,
):
    """Fetch article details for each PMID and write them to CSV and JSONL.

    The PubMed ``efetch`` requests are issued in batches controlled by
    ``batch_size`` to avoid overly long URL queries.
    """
    # Descarga artículos en lotes y guarda los resultados
    saved = 0
    csv_writer = None
    csvfile = None
    jsonlfile = None
    try:
        if output_csv:
            csv_exists = os.path.isfile(output_csv)
            csvfile = open(output_csv, "a", newline="", encoding="utf-8")
            fieldnames = list(Article.__annotations__.keys())
            csv_writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            if not csv_exists:
                csv_writer.writeheader()
        if output_jsonl:
            jsonlfile = open(output_jsonl, "a", encoding="utf-8")

        # Procesa los PMIDs en bloques
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
                response = session.get(fetch_url, params=params, timeout=10)
                response.raise_for_status()
                batch_root = ET.fromstring(response.text)
                for article in batch_root.findall("PubmedArticle"):
                    root = article

                    # Filtra por idiomas especificados
                    langs = [normalize_text(l.text).lower() for l in root.findall(".//Language")]
                    if languages and not any(l in languages for l in langs):
                        continue

                    pub_types_raw = [
                        normalize_text(pt.text)
                        for pt in root.findall(".//PublicationTypeList/PublicationType")
                    ]
                    pub_types = [pt.lower() for pt in pub_types_raw]
                    # Filtra por tipos de publicación seleccionados
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

                    pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

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

                    article = Article(
                        PMID=pmid,
                        Title=title,
                        Abstract=abstract,
                        Authors=authors,
                        Year=year,
                        Journal=journal,
                        Volume=volume,
                        Issue=issue,
                        Pages=pages,
                        DOI=doi,
                        citation_apa=citation,
                        PubMedURL=pubmed_url,
                        PublicationTypes=pub_types_raw,
                        QueryTerms=query_terms,
                    )
                    # Guarda el artículo en los formatos solicitados
                    if csv_writer:
                        row = asdict(article)
                        row["PublicationTypes"] = "; ".join(row["PublicationTypes"])
                        row["QueryTerms"] = "; ".join(row["QueryTerms"])
                        csv_writer.writerow(row)
                    if jsonlfile:
                        json.dump(asdict(article), jsonlfile, ensure_ascii=False)
                        jsonlfile.write("\n")
                    saved += 1
            except Exception as e:
                logger.error("Error retrieving PMIDs %s: %s", batch, e)
            time.sleep(0.1)
    finally:
        # Asegura el cierre de los archivos
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
    pub_types_pubmed: List[str],
    batch_size: int,
    max_results: int,
    api_key: str,
    format: str,
    append: bool,
    output_base: str,
):
    """Ejecuta la búsqueda en PubMed y guarda los resultados."""
    sanitized_terms = [sanitize_term(t) for t in terms]
    terms_query = f" {operator} ".join([f"({t})" for t in sanitized_terms])
    if pub_types_pubmed:
        query = f"{terms_query} AND ({' OR '.join(pub_types_pubmed)})"
    else:
        query = terms_query
    # Configura el logger y registra parámetros
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
            "pub_types_pubmed": pub_types_pubmed,
            "batch_size": batch_size,
            "max_results": max_results,
            "api_key_provided": bool(api_key),
            "format": format,
            "append": append,
            "output_base": output_base,
        },
    )

    # Obtiene PMIDs y prepara filtros
    pmids = search_pmids(query, max_results, api_key, sort=sort)
    article_types = [normalize_text(t).lower() for t in pub_types]
    languages = []

    # Determina archivos de salida
    output_csv = f"{output_base}.csv" if format in ("csv", "both") else None
    output_jsonl = f"{output_base}.jsonl" if format in ("jsonl", "both") else None

    # Elimina archivos previos si no se va a anexar
    if not append:
        if output_csv and os.path.exists(output_csv):
            os.remove(output_csv)
        if output_jsonl and os.path.exists(output_jsonl):
            os.remove(output_jsonl)

    # Descarga artículos y actualiza el README
    count = download_articles(
        pmids,
        api_key,
        output_csv,
        output_jsonl,
        article_types,
        languages,
        terms,
        batch_size,
    )
    write_readme(query, count, __version__)
    logger.info("Saved %d articles", count)
    return count


def main(**kwargs):
    """Wrapper para ejecutar run_query desde otros contextos."""
    return run_query(**kwargs)


if __name__ == "__main__":
    # Manejo de argumentos de línea de comandos
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
    parser.add_argument(
        "--output_base",
        default="papers",
        help="Base name for output files",
    )
    args = parser.parse_args()
    terms = [t.strip() for t in args.terms.split(";") if t.strip()]
    pub_types = [pt.strip() for pt in args.pub_types.split(",") if pt.strip()]
    pub_types_pubmed = [f"{pt}[Publication Type]" for pt in pub_types]
    main(
        terms=terms,
        operator=args.operator,
        sort=args.sort,
        pub_types=pub_types,
        pub_types_pubmed=pub_types_pubmed,
        batch_size=args.batch_size,
        max_results=args.max_results,
        api_key=args.api_key,
        format=args.format,
        append=args.append,
        output_base=args.output_base,
    )
