#!/usr/bin/env bash
# ===========================================
# Script: 00_search/geo_to_sra.sh
# Descripción: Discover GEO datasets (IBD/UC/CD + miRNA in human), resolve links to SRA/PubMed,
#              and export artifacts for downstream pipeline (runinfo.csv, JSONL summaries, provenance).
# Fecha: $(date +%Y-%m-%d)
# ===========================================

# --- Configuración segura ---
set -euo pipefail
IFS=$'\n\t'

# --- Variables (rutas relativas al repo; no edites si ejecutas desde PT/) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"          # .../PT
OUT_SEARCH="${ROOT_DIR}/outputs/search"
OUT_RUNS="${ROOT_DIR}/outputs/run_tables"
OUT_LOGS="${ROOT_DIR}/outputs/logs"
LOGFILE="${OUT_LOGS}/00_search.log"

# Query por defecto (puedes pasar otra como 1er argumento)
DEFAULT_QUERY='("inflammatory bowel disease"[All Fields] OR "ulcerative colitis"[All Fields] OR Crohn[All Fields]) AND (miRNA OR microRNA) AND "Homo sapiens"[Organism]'
QUERY="${1:-$DEFAULT_QUERY}"

# --- Funciones ---
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE" ; }

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    log "ERROR: required command '$1' not found in PATH"
    exit 127
  fi
}

# --- Inicio ---
mkdir -p "$OUT_SEARCH" "$OUT_RUNS" "$OUT_LOGS"
: > "$LOGFILE"
log "Starting 00_search (EDirect discovery) ..."
log "Root dir: $ROOT_DIR"
log "Outputs:  $OUT_SEARCH  |  $OUT_RUNS"
log "Query:    $QUERY"

# Validar dependencias
for c in esearch elink efetch esummary xtract; do need_cmd "$c"; done
if command -v jq >/dev/null 2>&1; then HAVE_JQ=1; else HAVE_JQ=0; fi

# 1) Buscar en GEO (db=gds) con historial
log "Running esearch (db=gds, usehistory=y) ..."
esearch -db gds -query "$QUERY" -usehistory y > "${OUT_SEARCH}/gds_search.xml" || {
  log "ERROR: esearch failed"
  exit 1
}

WEBENV="$(xtract -input "${OUT_SEARCH}/gds_search.xml" -pattern ENTREZ_DIRECT -element WebEnv || true)"
QK="$(xtract -input "${OUT_SEARCH}/gds_search.xml" -pattern ENTREZ_DIRECT -element QueryKey || true)"

if [[ -z "${WEBENV}" || -z "${QK}" ]]; then
  log "No WebEnv/QueryKey found — likely zero results. Aborting gracefully."
  echo -n > "${OUT_SEARCH}/gds_uids.txt"
  exit 0
fi
log "Captured WebEnv and QueryKey."

# 2) UIDs de GDS
log "Fetching GDS UIDs ..."
efetch -db gds -format uid -webenv "$WEBENV" -query_key "$QK" > "${OUT_SEARCH}/gds_uids.txt" || {
  log "WARN: efetch UIDs failed (continuing)."
  echo -n > "${OUT_SEARCH}/gds_uids.txt"
}
UID_COUNT="$(wc -l < "${OUT_SEARCH}/gds_uids.txt" | tr -d ' ')"
log "Found ${UID_COUNT} GDS UIDs."

# 3) Resolver GDS → SRA (runinfo.csv)
log "Linking GDS → SRA and exporting runinfo.csv ..."
elink -dbfrom gds -db sra -webenv "$WEBENV" -query_key "$QK" \
| efetch -format runinfo > "${OUT_RUNS}/gds_to_sra_runinfo.csv" || {
  log "WARN: elink/efetch runinfo failed (no SRA links?)."
  echo "Run" > "${OUT_RUNS}/gds_to_sra_runinfo.csv"
}

RUN_ROWS="$(wc -l < "${OUT_RUNS}/gds_to_sra_runinfo.csv" | tr -d ' ')"
if [[ "$RUN_ROWS" -le 1 ]]; then
  log "No SRA runs resolved (runinfo has header only)."
else
  log "runinfo.csv rows: ${RUN_ROWS}"
fi

# 4) Resolver GDS → PubMed
log "Linking GDS → PubMed PMIDs ..."
elink -dbfrom gds -db pubmed -webenv "$WEBENV" -query_key "$QK" \
| efetch -format uid > "${OUT_SEARCH}/gds_to_pubmed_pmids.txt" || {
  log "WARN: elink PubMed failed (continuing)."
  echo -n > "${OUT_SEARCH}/gds_to_pubmed_pmids.txt"
}
PMID_ROWS="$(wc -l < "${OUT_SEARCH}/gds_to_pubmed_pmids.txt" | tr -d ' ')"
log "PMIDs linked: ${PMID_ROWS}"

# 5) esummary por UID → JSONL
log "Writing esummary JSONL per GDS UID ..."
: > "${OUT_SEARCH}/datasets.jsonl"
if [[ "$UID_COUNT" -gt 0 ]]; then
  while read -r uid; do
    [[ -z "$uid" ]] && continue
    esummary -db gds -id "$uid" -json >> "${OUT_SEARCH}/datasets.jsonl" || true
  done < "${OUT_SEARCH}/gds_uids.txt"
fi
log "datasets.jsonl size: $(wc -c < "${OUT_SEARCH}/datasets.jsonl" | tr -d ' ') bytes"

# 6) (Opcional) índice TSV legible si jq está disponible
if [[ "$HAVE_JQ" -eq 1 ]]; then
  log "jq detected — creating datasets_index.tsv ..."
  jq -r '
    .result[]? | select(.acc) |
    [
      .acc,
      (.title // ""),
      (.gdsType // ""),
      (.n_samples // .ssCount // 0),
      (.taxname // ""),
      (.PDAT // .pdat // ""),
      (.gpl // "")
    ] | @tsv
  ' "${OUT_SEARCH}/datasets.jsonl" > "${OUT_SEARCH}/datasets_index.tsv" || true
  log "datasets_index.tsv rows: $(wc -l < "${OUT_SEARCH}/datasets_index.tsv" | tr -d ' ')"
else
  log "jq not found — skipping datasets_index.tsv (optional)."
fi

# 7) Provenance
log "Writing provenance.txt ..."
{
  echo "query: $QUERY"
  echo "date_utc: $(date -u +%FT%TZ)"
  echo "webenv: $WEBENV"
  echo "query_key: $QK"
} > "${OUT_SEARCH}/provenance.txt"

# 8) (Opcional) Enlaces a BioProject/BioSample usando la columna Run si existe
if [[ "$RUN_ROWS" -gt 1 ]]; then
  log "Optionally linking SRA → BioProject/BioSample (UID lists) ..."
  TMP_RUNS="$(mktemp)"
  # extrae columna 'Run' del CSV (header: Run, ...)
  awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i=="Run") c=i}} NR>1 && c{print $c}' \
    "${OUT_RUNS}/gds_to_sra_runinfo.csv" | sort -u > "$TMP_RUNS" || true

  if [[ -s "$TMP_RUNS" ]]; then
    # BioProject
    elink -dbfrom sra -db bioproject -id $(paste -sd, "$TMP_RUNS") \
    | efetch -format uid > "${OUT_SEARCH}/sra_to_bioproject.txt" || true
    # BioSample
    elink -dbfrom sra -db biosample  -id $(paste -sd, "$TMP_RUNS") \
    | efetch -format uid > "${OUT_SEARCH}/sra_to_biosample.txt"  || true
  fi
  rm -f "$TMP_RUNS"
fi

# --- Fin ---
log "00_search complete."
log "Artifacts:"
log "  - ${OUT_SEARCH}/gds_search.xml"
log "  - ${OUT_SEARCH}/gds_uids.txt"
log "  - ${OUT_RUNS}/gds_to_sra_runinfo.csv"
log "  - ${OUT_SEARCH}/gds_to_pubmed_pmids.txt"
log "  - ${OUT_SEARCH}/datasets.jsonl"
[[ -f "${OUT_SEARCH}/datasets_index.tsv" ]] && log "  - ${OUT_SEARCH}/datasets_index.tsv"
log "  - ${OUT_SEARCH}/provenance.txt"
