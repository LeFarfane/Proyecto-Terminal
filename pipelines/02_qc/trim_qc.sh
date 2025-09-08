#!/usr/bin/env bash
# ===========================================
# Script: 02_qc/trim_qc.sh
# Descripción: QC previo y posterior al trimming de lecturas small-RNA.
#              Corre FastQC (pre), recorta adaptador 3' con cutadapt y longitud mínima,
#              vuelve a correr FastQC (post) y genera un reporte MultiQC.
# Fecha: $(date +%Y-%m-%d)
# ===========================================

# --- Configuración segura ---
set -euo pipefail
IFS=$'\n\t'

# --- Variables ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"               # .../PT
LOG_DIR="${ROOT_DIR}/outputs/logs"
LOGFILE="${LOG_DIR}/02_qc.log"

# Parámetros (puedes sobreescribirlos por CLI)
# $1 = INPUT_DIR con *.fastq(.gz)
# $2 = TRIM_DIR de salida para *.trim.fastq.gz
# $3 = QC_DIR para reportes FastQC/MultiQC
# $4 = ADAPTER 3' (default Illumina small RNA)
# $5 = MINLEN (longitud mínima tras trimming)
# $6 = THREADS (para FastQC y compresión)
INPUT_DIR="${1:-${ROOT_DIR}/data_raw/sra_fastq}"
TRIM_DIR="${2:-${ROOT_DIR}/data_raw/trimmed}"
QC_DIR="${3:-${ROOT_DIR}/outputs/qc}"
ADAPTER="${4:-TGGAATTCTCGGGTGCCAAGG}"
MINLEN="${5:-15}"
THREADS="${6:-4}"

# --- Funciones ---
log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}
need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    log "ERROR: comando requerido no encontrado: $1"
    exit 127
  fi
}

# --- Inicio ---
mkdir -p "$LOG_DIR" "$TRIM_DIR" "$QC_DIR"
: > "$LOGFILE"
log "Iniciando script..."
log "INPUT_DIR = $INPUT_DIR"
log "TRIM_DIR  = $TRIM_DIR"
log "QC_DIR    = $QC_DIR"
log "ADAPTER   = $ADAPTER"
log "MINLEN    = $MINLEN"
log "THREADS   = $THREADS"

# Validar binarios
for c in fastqc cutadapt multiqc find xargs awk; do need_cmd "$c"; done
if command -v pigz >/dev/null 2>&1; then COMPRESSOR="pigz -p ${THREADS}"; else COMPRESSOR="gzip"; fi

# Validar directorio de entrada
if [ -d "$INPUT_DIR" ]; then
  log "Directorio de entrada encontrado: $INPUT_DIR"
else
  log "ERROR: Directorio de entrada no existe"
  exit 1
fi

# Listar FASTQ
FILES_COUNT="$(find "$INPUT_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) | wc -l | tr -d ' ')"
if [[ "$FILES_COUNT" -eq 0 ]]; then
  log "ERROR: No se encontraron archivos FASTQ en $INPUT_DIR"
  exit 1
fi
log "FASTQ detectados: $FILES_COUNT"

# --- QC previo ---
log "[QC pre] Ejecutando FastQC previo al trimming..."
# Ejecuta FastQC en paralelo por archivo (mejor que pasar todos a la vez)
find "$INPUT_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) -print0 \
| xargs -0 -I{} -P "$THREADS" bash -c 'fastqc -o "$0" "$1" >>"$2" 2>&1' "$QC_DIR" "{}" "$LOGFILE"
log "[QC pre] FastQC completado."

# --- Trimming ---
log "[trim] Iniciando recorte con cutadapt (adapter=$ADAPTER, minlen=$MINLEN)..."
find "$INPUT_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) -print0 \
| while IFS= read -r -d '' f; do
    bn="$(basename "$f")"
    # Normaliza nombre base (quita .fastq(.gz)/.fq(.gz))
    base="${bn%.fastq.gz}"; base="${base%.fq.gz}"; base="${base%.fastq}"; base="${base%.fq}"
    out="${TRIM_DIR}/${base}.trim.fastq.gz"
    log "[trim] ${bn} → ${out}"
    # Ejecuta cutadapt y comprime al vuelo
    # (cutadapt puede escribir gz directo si el nombre termina en .gz)
    cutadapt -a "$ADAPTER" -m "$MINLEN" -o "$out" "$f" >>"$LOGFILE" 2>&1
  done
log "[trim] Trimming completado."

# --- QC posterior ---
POST_COUNT="$(find "$TRIM_DIR" -type f -name "*.trim.fastq.gz" | wc -l | tr -d ' ')"
if [[ "$POST_COUNT" -eq 0 ]]; then
  log "ERROR: No se generaron archivos .trim.fastq.gz en $TRIM_DIR"
  exit 1
fi
log "[QC post] Ejecutando FastQC posterior al trimming sobre $POST_COUNT archivos..."
find "$TRIM_DIR" -type f -name "*.trim.fastq.gz" -print0 \
| xargs -0 -I{} -P "$THREADS" bash -c 'fastqc -o "$0" "$1" >>"$2" 2>&1' "$QC_DIR" "{}" "$LOGFILE"
log "[QC post] FastQC posterior completado."

# --- MultiQC ---
log "[multiqc] Generando reporte consolidado en $QC_DIR ..."
multiqc -o "$QC_DIR" "$QC_DIR" >>"$LOGFILE" 2>&1 || {
  log "WARN: MultiQC report falló (continúo)."
}

# --- Fin ---
log "Script finalizado correctamente."
