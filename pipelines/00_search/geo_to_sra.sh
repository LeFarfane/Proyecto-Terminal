#!/usr/bin/env bash
# ===========================================
# Script: 00_search/geo_to_sra.sh
# Descripción: Discover GEO datasets (IBD/UC/CD + miRNA in human), resolver SRA/PubMed
#              y exportar runinfo.csv, JSONL y provenance.
# Nota: Modo híbrido. Usa EDirect si existe; si no, fallback con curl+Python (Replit-friendly).
# Fecha: $(date +%Y-%m-%d)
# ===========================================

set -euo pipefail
IFS=$'\n\t'

export NCBI_EDIRECT_USE_CURL=1  # force curl backend for reliability

# Prefer user-installed EDirect in $HOME/edirect
EDIRECT_HOME="${EDIRECT_HOME:-$HOME/edirect}"
if [[ -d "$EDIRECT_HOME" ]]; then
  export PATH="$EDIRECT_HOME:$PATH"
fi

# --- Rutas ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_SEARCH="${ROOT_DIR}/outputs/search"
OUT_RUNS="${ROOT_DIR}/outputs/run_tables"
OUT_LOGS="${ROOT_DIR}/outputs/logs"
LOGFILE="${OUT_LOGS}/00_search.log"

# --- Parámetros ---
# Default search term (single line, no newlines)
DEFAULT_QUERY=$(cat <<'EOF' | tr -d '\n'
("inflammatory bowel disease"[All Fields] OR "ulcerative colitis"[All Fields] OR Crohn[All Fields]) AND (miRNA OR microRNA) AND "Homo sapiens"[Organism]
EOF
)
QUERY="${1:-$DEFAULT_QUERY}"
MODE="${2:-overwrite}"   # overwrite | append
EFETCHSTOP="${EFETCHSTOP:-100000}" # efetch -stop limit

# --- Utils ---
timestamp(){ date +"%Y-%m-%d %H:%M:%S"; }
log(){ echo "[$(timestamp)] $*" | tee -a "$LOGFILE"; }

have(){ command -v "$1" >/dev/null 2>&1; }
need(){ if ! have "$1"; then log "ERROR: required command '$1' not found in PATH"; exit 127; fi; }

# Detect full EDirect suite and ensure expected options exist
have_edirect(){
  for cmd in esearch efetch elink esummary xtract; do
    command -v "$cmd" >/dev/null 2>&1 || return 1
  done
  efetch -help 2>&1 | grep -qi 'webenv' || return 1
  elink  -help 2>&1 | grep -qi 'dbfrom' || return 1
  return 0
}

FALLBACK_REASON="EDirect not found"

mkdir -p "$OUT_SEARCH" "$OUT_RUNS" "$OUT_LOGS"
[[ "$MODE" != "append" ]] && : > "$LOGFILE"

log "Starting 00_search (hybrid EDirect/curl) ..."
log "Root dir: $ROOT_DIR"
log "Outputs:  $OUT_SEARCH  |  $OUT_RUNS"
log "Query:    $QUERY"

# -------------------------------------------
# RUTA A) Con EDirect (si está disponible)
# -------------------------------------------
if have_edirect; then
  # Resolve absolute paths to avoid PATH collisions
  ESEARCH="$(command -v esearch)"
  EFETCH="$(command -v efetch)"
  ELINK="$(command -v elink)"
  ESUMMARY="$(command -v esummary)"
  XTRACT="$(command -v xtract)"

  log "EDirect detected — using native esearch/elink/efetch."
  log "esearch path: ${ESEARCH}"
  log "efetch path:  ${EFETCH}"
  log "elink path:   ${ELINK}"
  log "esearch version: $(${ESEARCH} -version 2>&1 | head -n1)"
  log "efetch version:  $(${EFETCH} -version 2>&1 | head -n1)"
  log "elink version:   $(${ELINK} -version 2>&1 | head -n1)"
  HAVE_JQ=0; have jq && HAVE_JQ=1
  EDIRECT_OK=1

  log "Fetching GDS UIDs ..."
  ${ESEARCH} -db gds -query "$QUERY" | ${EFETCH} -format uid > "${OUT_SEARCH}/gds_uids.txt" 2>> "$LOGFILE" || { log "WARN: efetch UIDs failed"; : > "${OUT_SEARCH}/gds_uids.txt"; EDIRECT_OK=0; }
  UID_COUNT="$(wc -l < "${OUT_SEARCH}/gds_uids.txt" | tr -d ' ')"
  log "Found ${UID_COUNT} GDS UIDs."

  log "Linking GDS → SRA and exporting runinfo.csv ..."
  ${ESEARCH} -db gds -query "$QUERY" | ${ELINK} -target sra | ${EFETCH} -format runinfo > "${OUT_RUNS}/gds_to_sra_runinfo.csv" 2>> "$LOGFILE" || { log "WARN: elink/efetch runinfo failed"; echo "Run" > "${OUT_RUNS}/gds_to_sra_runinfo.csv"; EDIRECT_OK=0; }
  RUN_ROWS="$(wc -l < "${OUT_RUNS}/gds_to_sra_runinfo.csv" | tr -d ' ')"
  if [[ "$RUN_ROWS" -le 1 ]]; then log "No SRA runs resolved (runinfo header only)."; else log "runinfo.csv rows: ${RUN_ROWS}"; fi

  log "Linking GDS → PubMed PMIDs ..."
  ${ESEARCH} -db gds -query "$QUERY" | ${ELINK} -target pubmed | ${EFETCH} -format uid > "${OUT_SEARCH}/gds_to_pubmed_pmids.txt" 2>> "$LOGFILE" || { log "WARN: elink PubMed failed"; : > "${OUT_SEARCH}/gds_to_pubmed_pmids.txt"; EDIRECT_OK=0; }
  PMID_ROWS="$(wc -l < "${OUT_SEARCH}/gds_to_pubmed_pmids.txt" | tr -d ' ')"
  log "PMIDs linked: ${PMID_ROWS}"

  log "Writing esummary JSONL per GDS UID ..."
  [[ "$MODE" != "append" ]] && : > "${OUT_SEARCH}/datasets.jsonl"
  if [[ "$UID_COUNT" -gt 0 ]]; then
    while read -r uid; do
      [[ -z "$uid" ]] && continue
      ${ESUMMARY} -db gds -id "$uid" -json >> "${OUT_SEARCH}/datasets.jsonl" || true
    done < "${OUT_SEARCH}/gds_uids.txt"
  fi
  log "datasets.jsonl size: $(wc -c < "${OUT_SEARCH}/datasets.jsonl" | tr -d ' ') bytes"

  if [[ "$HAVE_JQ" -eq 1 ]]; then
    log "jq detected — creating datasets_index.tsv ..."
    if [[ "$MODE" == "append" ]]; then
      jq -r '.result[]? | select(.acc) | [ .acc, (.title // ""), (.gdsType // ""), (.n_samples // .ssCount // 0), (.taxname // ""), (.PDAT // .pdat // ""), (.gpl // "") ] | @tsv' \
        "${OUT_SEARCH}/datasets.jsonl" >> "${OUT_SEARCH}/datasets_index.tsv" || true
    else
      jq -r '.result[]? | select(.acc) | [ .acc, (.title // ""), (.gdsType // ""), (.n_samples // .ssCount // 0), (.taxname // ""), (.PDAT // .pdat // ""), (.gpl // "") ] | @tsv' \
        "${OUT_SEARCH}/datasets.jsonl" > "${OUT_SEARCH}/datasets_index.tsv" || true
    fi
    log "datasets_index.tsv rows: $(wc -l < "${OUT_SEARCH}/datasets_index.tsv" | tr -d ' ')"
  else
    log "jq not found — skipping datasets_index.tsv (optional)."
  fi

  log "Writing provenance.txt ..."
  {
    echo "query: $QUERY"
    echo "date_utc: $(date -u +%FT%TZ)"
  } > "${OUT_SEARCH}/provenance.txt"

  if [[ "$RUN_ROWS" -gt 1 ]]; then
    log "Optionally linking SRA → BioProject/BioSample (UID lists) ..."
    TMP_RUNS="$(mktemp)"
    awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i=="Run") c=i}} NR>1 && c{print $c}' "${OUT_RUNS}/gds_to_sra_runinfo.csv" | sort -u > "$TMP_RUNS" || true
    if [[ -s "$TMP_RUNS" ]]; then
      ${ELINK} -dbfrom sra -db bioproject -id $(paste -sd, "$TMP_RUNS") | ${EFETCH} -format uid > "${OUT_SEARCH}/sra_to_bioproject.txt" || true
      ${ELINK} -dbfrom sra -db biosample  -id $(paste -sd, "$TMP_RUNS") | ${EFETCH} -format uid > "${OUT_SEARCH}/sra_to_biosample.txt"  || true
    fi
    rm -f "$TMP_RUNS"
  fi

  if [[ "$EDIRECT_OK" -eq 1 ]]; then
    log "00_search complete."
    log "Artifacts:"
    log "  - ${OUT_SEARCH}/gds_uids.txt"
    log "  - ${OUT_RUNS}/gds_to_sra_runinfo.csv"
    log "  - ${OUT_SEARCH}/gds_to_pubmed_pmids.txt"
    log "  - ${OUT_SEARCH}/datasets.jsonl"
    [[ -f "${OUT_SEARCH}/datasets_index.tsv" ]] && log "  - ${OUT_SEARCH}/datasets_index.tsv"
    log "  - ${OUT_SEARCH}/provenance.txt"
    exit 0
  else
    FALLBACK_REASON="EDirect pipeline failed"
  fi
fi

# -------------------------------------------
# RUTA B) Fallback sin EDirect (curl + Python)
# -------------------------------------------
log "$FALLBACK_REASON — using curl + Python fallback (compatible con Replit)."

HAVE_JQ=0; have jq && HAVE_JQ=1

# 1) esearch vía E-utilities (JSON)
ESEARCH_JSON="${OUT_SEARCH}/gds_search.json"
log "Running esearch via curl (history mode) ..."
curl -sS -G "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
  --data-urlencode "db=gds" \
  --data-urlencode "retmode=json" \
  --data-urlencode "usehistory=y" \
  --data-urlencode "term=${QUERY}" \
  -o "$ESEARCH_JSON" || { log "ERROR: esearch curl failed"; exit 1; }

# Extract WebEnv and QueryKey for subsequent efetch
read WEBENV QK COUNT < <(
  python3 - "$ESEARCH_JSON" <<'PY'
import sys, json
js = sys.argv[1]
try:
    data = json.load(open(js, "r", encoding="utf-8"))
    res = data.get("esearchresult", {})
    print(res.get("webenv", ""), res.get("querykey", ""), res.get("count", "0"))
except Exception:
    print("", "", "0")
PY
)

if [[ -z "$WEBENV" || -z "$QK" ]]; then
  log "ERROR: ESearch response invalid"; exit 1
fi

log "esearch count: ${COUNT}"
EFETCH_UIDS_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&rettype=uilist&retstart=0&retmax=${EFETCHSTOP}&query_key=${QK}&WebEnv=${WEBENV}"
curl -sS "$EFETCH_UIDS_URL" -o "${OUT_SEARCH}/gds_uids.txt" || { log "ERROR: efetch curl failed"; : > "${OUT_SEARCH}/gds_uids.txt"; }

UID_COUNT=$( [ -f "${OUT_SEARCH}/gds_uids.txt" ] && wc -l < "${OUT_SEARCH}/gds_uids.txt" | tr -d ' ' || echo 0 )
log "Found ${UID_COUNT} GDS UIDs."

  if [[ "$UID_COUNT" -eq 0 ]]; then
    echo "Run" > "${OUT_RUNS}/gds_to_sra_runinfo.csv"
    : > "${OUT_SEARCH}/gds_to_pubmed_pmids.txt"
    : > "${OUT_SEARCH}/datasets.jsonl"
    {
      echo "query: $QUERY"
      echo "date_utc: $(date -u +%FT%TZ)"
      echo "webenv: $WEBENV"
      echo "query_key: $QK"
    } > "${OUT_SEARCH}/provenance.txt"
    log "00_search complete (sin resultados)."
    exit 0
  fi

# 2) elink GDS→SRA y GDS→PubMed (JSON), en batches
mapfile -t GDS_IDS < "${OUT_SEARCH}/gds_uids.txt"
BATCH=200
join_csv(){ local IFS=,; echo "$*"; }

SRA_UIDS_TMP="$(mktemp)"; : > "$SRA_UIDS_TMP"
PMIDS_TMP="$(mktemp)";    : > "$PMIDS_TMP"

total=${#GDS_IDS[@]}; start=0
while (( start < total )); do
  end=$(( start + BATCH )); (( end > total )) && end=$total
  batch=( "${GDS_IDS[@]:start:end-start}" )
  ids_csv=$(join_csv "${batch[@]}")

  ELINK_SRA_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?retmode=json&dbfrom=gds&db=sra&id=${ids_csv}"
  curl -sS "$ELINK_SRA_URL" -o "${OUT_SEARCH}/elink_gds2sra_${start}.json" || true
  python3 - "${OUT_SEARCH}/elink_gds2sra_${start}.json" "$SRA_UIDS_TMP" << 'PY'
import sys, json
js, outp = sys.argv[1], sys.argv[2]
try:
    data = json.load(open(js, "r", encoding="utf-8"))
except Exception:
    sys.exit(0)
for ls in data.get("linksets", []):
    for dbs in ls.get("linksetdbs", []):
        if dbs.get("dbto")=="sra":
            ids=[str(i) for i in dbs.get("links", [])]
            if ids:
                with open(outp,"a",encoding="utf-8") as f: f.write("\n".join(ids)+"\n")
PY

  ELINK_PUB_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?retmode=json&dbfrom=gds&db=pubmed&id=${ids_csv}"
  curl -sS "$ELINK_PUB_URL" -o "${OUT_SEARCH}/elink_gds2pubmed_${start}.json" || true
  python3 - "${OUT_SEARCH}/elink_gds2pubmed_${start}.json" "$PMIDS_TMP" << 'PY'
import sys, json
js, outp = sys.argv[1], sys.argv[2]
try:
    data = json.load(open(js, "r", encoding="utf-8"))
except Exception:
    sys.exit(0)
for ls in data.get("linksets", []):
    for dbs in ls.get("linksetdbs", []):
        if dbs.get("dbto")=="pubmed":
            ids=[str(i) for i in dbs.get("links", [])]
            if ids:
                with open(outp,"a",encoding="utf-8") as f: f.write("\n".join(ids)+"\n")
PY

  log "Processed GDS elinks: ${end}/${total}"
  start=$end
done

sort -u "$PMIDS_TMP" > "${OUT_SEARCH}/gds_to_pubmed_pmids.txt" || :
PMID_ROWS="$(wc -l < "${OUT_SEARCH}/gds_to_pubmed_pmids.txt" | tr -d ' ')"
log "PMIDs linked: ${PMID_ROWS}"

# 3) Para SRA: de UIDs SRA → SRR (ESummary) → runinfo.csv
SRR_TMP="$(mktemp)"; : > "$SRR_TMP"
if [[ -s "$SRA_UIDS_TMP" ]]; then
  mapfile -t SRA_UIDS < <(sort -u "$SRA_UIDS_TMP")
  start=0; total=${#SRA_UIDS[@]}
  while (( start < total )); do
    end=$(( start + BATCH )); (( end > total )) && end=$total
    batch=( "${SRA_UIDS[@]:start:end-start}" )
    ids_csv=$(join_csv "${batch[@]}")

    ESUM_SRA_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&retmode=json&id=${ids_csv}"
    ESUM_SRA_JSON="${OUT_SEARCH}/esum_sra_${start}.json"
    curl -sS "$ESUM_SRA_URL" -o "$ESUM_SRA_JSON" || true

    python3 - "$ESUM_SRA_JSON" "$SRR_TMP" << 'PY'
import sys, json, xml.etree.ElementTree as ET, re
js, outp = sys.argv[1], sys.argv[2]
try:
    data = json.load(open(js,"r",encoding="utf-8"))
except Exception:
    sys.exit(0)
result = data.get("result", {})
uids = result.get("uids", [])
srrs=set()
for uid in uids:
    obj = result.get(uid, {})
    expxml = obj.get("expxml") or obj.get("expxml1") or ""
    if not expxml: continue
    try:
        root = ET.fromstring(expxml)
        for run in root.findall(".//RUN"):
            acc = run.get("accession")
            if acc and acc.startswith("SRR"): srrs.add(acc)
    except Exception:
        srrs.update(re.findall(r'\bSRR\d+\b', expxml))
with open(outp,"a",encoding="utf-8") as f:
    for s in sorted(srrs): f.write(s+"\n")
PY

    log "Parsed SRA summaries for UIDs ${start}-${end}."
    start=$end
  done
fi

RUNINFO_CSV="${OUT_RUNS}/gds_to_sra_runinfo.csv"
if [[ -s "$SRR_TMP" ]]; then
  SRR_LIST="$(paste -sd, "$SRR_TMP")"
  RUNINFO_URL="https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&acc=${SRR_LIST}"
  curl -sS "$RUNINFO_URL" -o "$RUNINFO_CSV" || echo "Run" > "$RUNINFO_CSV"
else
  echo "Run" > "$RUNINFO_CSV"
fi

RUN_ROWS="$(wc -l < "$RUNINFO_CSV" | tr -d ' ')"
if [[ "$RUN_ROWS" -le 1 ]]; then log "No SRA runs resolved (runinfo header only)."; else log "runinfo.csv rows: ${RUN_ROWS}"; fi

# 4) ESummary GDS → JSONL (por lotes)
[[ "$MODE" != "append" ]] && : > "${OUT_SEARCH}/datasets.jsonl"
mapfile -t GDS_IDS < "${OUT_SEARCH}/gds_uids.txt"
start=0; total=${#GDS_IDS[@]}
while (( start < total )); do
  end=$(( start + BATCH )); (( end > total )) && end=$total
  batch=( "${GDS_IDS[@]:start:end-start}" )
  ids_csv=$(join_csv "${batch[@]}")
  ESUM_GDS_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json&id=${ids_csv}"
  ESUM_GDS_JSON="${OUT_SEARCH}/esum_gds_${start}.json"
  curl -sS "$ESUM_GDS_URL" -o "$ESUM_GDS_JSON" || true

  python3 - "$ESUM_GDS_JSON" "${OUT_SEARCH}/datasets.jsonl" << 'PY'
import sys, json
js, outp = sys.argv[1], sys.argv[2]
try:
    data = json.load(open(js,"r",encoding="utf-8"))
except Exception:
    sys.exit(0)
res = data.get("result", {})
uids = res.get("uids", [])
with open(outp,"a",encoding="utf-8") as f:
    for uid in uids:
        obj = res.get(uid, {})
        if obj: f.write(json.dumps(obj, ensure_ascii=False)+"\n")
PY

  log "Wrote GDS summaries: ${end}/${total}"
  start=$end
done
log "datasets.jsonl size: $(wc -c < "${OUT_SEARCH}/datasets.jsonl" | tr -d ' ') bytes"

# 5) datasets_index.tsv (opcional)
if [[ "$HAVE_JQ" -eq 1 ]]; then
  log "jq detected — creating datasets_index.tsv ..."
  if [[ "$MODE" == "append" ]]; then
    jq -r 'fromjson? | select(.acc != null) | [ .acc, (.title // ""), (.gdsType // ""), (.n_samples // .ssCount // 0), (.taxname // ""), (.PDAT // .pdat // ""), (.gpl // "") ] | @tsv' \
      "${OUT_SEARCH}/datasets.jsonl" >> "${OUT_SEARCH}/datasets_index.tsv" || true
  else
    jq -r 'fromjson? | select(.acc != null) | [ .acc, (.title // ""), (.gdsType // ""), (.n_samples // .ssCount // 0), (.taxname // ""), (.PDAT // .pdat // ""), (.gpl // "") ] | @tsv' \
      "${OUT_SEARCH}/datasets.jsonl" > "${OUT_SEARCH}/datasets_index.tsv" || true
  fi
  log "datasets_index.tsv rows: $(wc -l < "${OUT_SEARCH}/datasets_index.tsv" | tr -d ' ')"
else
  log "jq not found — skipping datasets_index.tsv (optional)."
fi

# 6) provenance
{
  echo "query: $QUERY"
  echo "date_utc: $(date -u +%FT%TZ)"
  echo "webenv: $WEBENV"
  echo "query_key: $QK"
} > "${OUT_SEARCH}/provenance.txt"

# 7) opcional SRA→BioProject/BioSample
if [[ "$RUN_ROWS" -gt 1 ]]; then
  log "Optionally linking SRA → BioProject/BioSample (UID lists) ..."
  TMP_RUNS="$(mktemp)"
  awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i=="Run") c=i}} NR>1 && c{print $c}' "$RUNINFO_CSV" | sort -u > "$TMP_RUNS" || true
  if [[ -s "$TMP_RUNS" ]]; then
    SRR_CSV="$(paste -sd, "$TMP_RUNS")"
    curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=sra&db=bioproject&retmode=json&id=${SRR_CSV}" -o "${OUT_SEARCH}/sra_to_bioproject.json" || true
    curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=sra&db=biosample&retmode=json&id=${SRR_CSV}"  -o "${OUT_SEARCH}/sra_to_biosample.json"  || true
    python3 - "${OUT_SEARCH}/sra_to_bioproject.json" "${OUT_SEARCH}/sra_to_bioproject.txt" << 'PY'
import sys, json, pathlib
js, outp = sys.argv[1], sys.argv[2]
uids=set()
try:
    data=json.load(open(js,'r',encoding='utf-8'))
    for ls in data.get("linksets", []):
        for dbs in ls.get("linksetdbs", []):
            if dbs.get("dbto")=="bioproject":
                uids.update(map(str, dbs.get("links", [])))
except Exception: pass
pathlib.Path(outp).write_text("\n".join(sorted(uids))+"\n" if uids else "", encoding="utf-8")
PY
    python3 - "${OUT_SEARCH}/sra_to_biosample.json" "${OUT_SEARCH}/sra_to_biosample.txt" << 'PY'
import sys, json, pathlib
js, outp = sys.argv[1], sys.argv[2]
uids=set()
try:
    data=json.load(open(js,'r',encoding='utf-8'))
    for ls in data.get("linksets", []):
        for dbs in ls.get("linksetdbs", []):
            if dbs.get("dbto")=="biosample":
                uids.update(map(str, dbs.get("links", [])))
except Exception: pass
pathlib.Path(outp).write_text("\n".join(sorted(uids))+"\n" if uids else "", encoding="utf-8")
PY
  fi
  rm -f "$TMP_RUNS"
fi

log "00_search complete."
log "Artifacts:"
log "  - ${OUT_SEARCH}/gds_search.json (fallback) o gds_search.xml (EDirect)"
log "  - ${OUT_SEARCH}/gds_uids.txt"
log "  - ${OUT_RUNS}/gds_to_sra_runinfo.csv"
log "  - ${OUT_SEARCH}/gds_to_pubmed_pmids.txt"
log "  - ${OUT_SEARCH}/datasets.jsonl"
[[ -f "${OUT_SEARCH}/datasets_index.tsv" ]] && log "  - ${OUT_SEARCH}/datasets_index.tsv"
log "  - ${OUT_SEARCH}/provenance.txt"
