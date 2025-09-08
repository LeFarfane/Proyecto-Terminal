#!/usr/bin/env bash
# ===========================================
# Script: 01_download/download_runs.sh
# Descripción: Descarga lecturas crudas de SRA (SRR) a partir de runinfo.csv o lista SRR.
#              Usa prefetch + fasterq-dump (con reintentos) y comprime con pigz/gzip.
# Fecha: $(date +%Y-%m-%d)
# ===========================================

# --- Configuración segura ---
set -euo pipefail
IFS=$'\n\t'

# --- Variables ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"             # .../PT
OUT_LOGS="${ROOT_DIR}/outputs/logs"
LOGFILE="${OUT_LOGS}/01_download.log"

# Args:
#   $1 = ruta a runinfo.csv (con columna 'Run') o archivo con SRR por línea
#   $2 = directorio de salida para FASTQ(.gz)
#   $3 = hilos para fasterq-dump (opcional; default 4)
#   $4 = trabajos en paralelo (opcional; default 3)
INPUT_RUNS="${1:-${ROOT_DIR}/outputs/run_tables/gds_to_sra_runinfo.csv}"
OUTPUT_DIR="${2:-${ROOT_DIR}/data_raw/sra_fastq}"
FQ_THREADS="${3:-4}"
PAR_JOBS="${4:-3}"

# --- Funciones ---
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE" ; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || { log "ERROR: comando requerido no encontrado: $1"; exit 127; }; }

collect_srrs() {
  local in="$1"
  [[ -s "$in" ]] || { log "ERROR: archivo de entrada vacío o inexistente: $in"; exit 1; }
  if head -1 "$in" | grep -qi '^Run[, ]'; then
    # CSV con cabecera 'Run'
    awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i=="Run") c=i}} NR>1 && c{print $c}' "$in" \
      | sed 's/\r$//' | grep -E '^SRR[0-9]+' | sort -u
  else
    # Lista simple de SRR
    sed 's/\r$//' "$in" | grep -E '^SRR[0-9]+' | sort -u
  fi
}

download_one() {
  local SRR="$1" outdir="$2" threads="$3"

  # ya existe comprimido
  if ls "${outdir}/${SRR}"*.fastq.gz >/dev/null 2>&1; then
    log "[skip] ${SRR} ya existe (.fastq.gz)"
    return 0
  fi
  # existe sin comprimir → comprimir
  if ls "${outdir}/${SRR}"*.fastq >/dev/null 2>&1; then
    if command -v pigz >/dev/null 2>&1; then pigz -f -p "${threads}" "${outdir}/${SRR}"*.fastq
    else gzip -f "${outdir}/${SRR}"*.fastq; fi
    log "[compress] ${SRR} comprimido"
    return 0
  fi

  log "[fetch] ${SRR} → prefetch"
  if ! prefetch "${SRR}" >>"$LOGFILE" 2>&1; then
    log "[warn] prefetch falló para ${SRR} — intento directo con fasterq-dump"
  fi

  log "[dump] ${SRR} → fasterq-dump (-e ${threads})"
  if ! fasterq-dump "${SRR}" -O "${outdir}" -e "${threads}" -p --split-files >>"$LOGFILE" 2>&1; then
    log "[retry] reintento fasterq-dump para ${SRR}"
    sleep 2
    fasterq-dump "${SRR}" -O "${outdir}" -e "${threads}" -p --split-files >>"$LOGFILE" 2>&1
  fi

  # comprimir salida
  if command -v pigz >/dev/null 2>&1; then pigz -p "${threads}" "${outdir}/${SRR}"*.fastq
  else gzip "${outdir}/${SRR}"*.fastq; fi

  log "[done] ${SRR} descargado y comprimido"
}

export -f download_one
export LOGFILE

# --- Inicio ---
mkdir -p "$OUT_LOGS" "$OUTPUT_DIR"
: > "$LOGFILE"
log "Iniciando script..."
log "Entrada: $INPUT_RUNS"
log "Salida : $OUTPUT_DIR"
log "Hilos fasterq: $FQ_THREADS | Trabajos paralelos: $PAR_JOBS"

# Validar binarios requeridos
for c in prefetch fasterq-dump awk sort xargs; do need_cmd "$c"; done
command -v pigz >/dev/null 2>&1 || log "[info] pigz no encontrado; usaré gzip"

# Configuración opcional de SRA (cache local)
export NCBI_SETTINGS="${NCBI_SETTINGS:-${ROOT_DIR}/.ncbi}"
mkdir -p "$NCBI_SETTINGS"

# Recolectar SRR y guardarlos para traza
TMP_SRR="$(mktemp)"
collect_srrs "$INPUT_RUNS" > "$TMP_SRR"
TOTAL_SRR="$(wc -l < "$TMP_SRR" | tr -d ' ')"
[[ "$TOTAL_SRR" -gt 0 ]] || { log "ERROR: no se encontraron SRR en la entrada."; rm -f "$TMP_SRR"; exit 1; }
log "Total de SRR a descargar: $TOTAL_SRR"

# Descargar en paralelo
log "Descargando corridas SRA..."
xargs -I{} -P "${PAR_JOBS}" bash -c 'download_one "$@"' _ {} "$OUTPUT_DIR" "$FQ_THREADS" < "$TMP_SRR"

# Limpieza temporal
rm -f "$TMP_SRR"

# Verificación final
DL_COUNT=$(ls -1 "${OUTPUT_DIR}"/SRR*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
if [[ "$DL_COUNT" -gt 0 ]]; then
  log "Descarga completa. Archivos .fastq.gz generados: $DL_COUNT"
else
  log "WARN: No se encontraron .fastq.gz en ${OUTPUT_DIR}. Revisa el log."
fi

# --- Fin ---
log "Script finalizado correctamente."
