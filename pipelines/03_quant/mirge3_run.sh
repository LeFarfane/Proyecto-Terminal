#!/usr/bin/env bash
# ===========================================
# Script: 03_quant/mirge3_run.sh
# Descripción: Cuantifica miRNAs por muestra usando miRge3 a partir de FASTQ recortados.
#              Detecta archivos *.trim.fastq.gz en data_raw/trimmed, ejecuta miRge3 por muestra,
#              y guarda resultados en results/quant_mirge3/<sample>/ (incluye TSV de conteos).
# Fecha: $(date +%Y-%m-%d)
# ===========================================

# --- Configuración segura ---
set -euo pipefail
IFS=$'\n\t'

# --- Variables ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"                 # .../PT
LOG_DIR="${ROOT_DIR}/outputs/logs"
LOGFILE="${LOG_DIR}/03_quant_mirge3.log"

# Parámetros con defaults (puedes sobreescribir por CLI):
# $1 = TRIM_DIR (entrada)         | default: PT/data_raw/trimmed
# $2 = OUT_DIR  (salida)          | default: PT/results/quant_mirge3
# $3 = SPECIES (miRge3 --species) | default: hsa
# $4 = ADAPTER (3' adapter)       | default: TGGAATTCTCGGGTGCCAAGG
# $5 = THREADS                     | default: 4
TRIM_DIR="${1:-${ROOT_DIR}/data_raw/trimmed}"
OUT_DIR="${2:-${ROOT_DIR}/results/quant_mirge3}"
SPECIES="${3:-hsa}"
ADAPTER="${4:-TGGAATTCTCGGGTGCCAAGG}"
THREADS="${5:-4}"

# --- Funciones ---
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE" ; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || { log "ERROR: comando requerido no encontrado: $1"; exit 127; }; }

# --- Inicio ---
mkdir -p "$LOG_DIR" "$OUT_DIR"
: > "$LOGFILE"
log "Iniciando script..."
log "TRIM_DIR = $TRIM_DIR"
log "OUT_DIR  = $OUT_DIR"
log "SPECIES  = $SPECIES"
log "ADAPTER  = $ADAPTER"
log "THREADS  = $THREADS"

# Validar binarios requeridos
for c in mirge3 awk find xargs; do need_cmd "$c"; done

# Validar directorio de entrada
if [ -d "$TRIM_DIR" ]; then
  log "Directorio de entrada encontrado: $TRIM_DIR"
else
  log "ERROR: Directorio de entrada no existe"
  exit 1
fi

# Listar archivos .trim.fastq(.gz)
mapfile -t TRIM_FILES < <(find "$TRIM_DIR" -maxdepth 1 -type f \( -name "*.trim.fastq.gz" -o -name "*.trim.fq.gz" -o -name "*.trim.fastq" -o -name "*.trim.fq" \) | sort)
TOTAL="${#TRIM_FILES[@]}"
if [[ "$TOTAL" -eq 0 ]]; then
  log "ERROR: No se encontraron archivos *.trim.fastq* en $TRIM_DIR"
  exit 1
fi
log "Archivos detectados para cuantificación: $TOTAL"

# Función para procesar una muestra
process_one() {
  local in_f="$1"
  local out_root="$2"
  local species="$3"
  local adapter="$4"
  local threads="$5"
  local main_log="$6"

  local bn; bn="$(basename "$in_f")"
  # Derivar nombre de muestra quitando sufijos comunes
  local base="${bn%.trim.fastq.gz}"; base="${base%.trim.fq.gz}"; base="${base%.trim.fastq}"; base="${base%.trim.fq}"
  # Si el nombre lleva _1/_R1, normalizar (opcional)
  base="${base%_1}"; base="${base%_R1}"

  local out_dir="${out_root}/${base}"
  mkdir -p "$out_dir"

  # Skip si ya hay conteos en el directorio
  if ls "${out_dir}"/*.{tsv,txt,csv} >/dev/null 2>&1; then
    # pero solo saltar si encontramos algo que parezca conteos
    if ls "${out_dir}"/miRNA*.{tsv,txt,csv} >/dev/null 2>&1 || ls "${out_dir}"/*count*.{tsv,txt,csv} >/dev/null 2>&1; then
      echo "[$(date +'%Y-%m-%d %H:%M:%S')] [skip] ${base} ya tiene archivos de conteo en ${out_dir}" | tee -a "$main_log"
      return 0
    fi
  fi

  echo "[$(date +'%Y-%m-%d %H:%M:%S')] [run] ${base} → miRge3" | tee -a "$main_log"

  # Log específico por muestra
  local smp_log="${out_dir}/${base}.mirge3.log"

  # Ejecutar miRge3
  # Notas:
  #  -s: archivo single-end
  #  -o: directorio de salida de la muestra
  #  --species: hsa (humano)
  #  --adapter: adaptador 3' (si aplica)
  #  --threads: hilos
  #  (miRge3 gestiona referencias automáticamente por especie; si ya están descargadas, las reutiliza)
  if ! mirge3 \
      -s "$in_f" \
      -o "$out_dir" \
      --species "$species" \
      --adapter "$adapter" \
      --threads "$threads" \
      >"$smp_log" 2>&1; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [error] miRge3 falló para ${base} (ver ${smp_log})" | tee -a "$main_log"
    return 1
  fi

  # Estandarizar nombre del archivo de conteos si es posible:
  # Buscamos un TSV/TXT/CSV que contenga 'mir' y 'count' en el nombre, y lo copiamos como miRNA.counts.tsv
  local counts_found
  counts_found="$(find "$out_dir" -maxdepth 1 -type f \( -name "*mir*count*.tsv" -o -name "*miR*count*.tsv" -o -name "*count*.tsv" -o -name "*mir*.tsv" \) | head -n1 || true)"
  if [[ -n "$counts_found" ]]; then
    cp -f "$counts_found" "${out_dir}/miRNA.counts.tsv"
  else
    # intenta con *.txt o *.csv si no hubo .tsv
    counts_found="$(find "$out_dir" -maxdepth 1 -type f \( -name "*mir*count*.txt" -o -name "*miR*count*.txt" -o -name "*count*.txt" -o -name "*mir*.txt" -o -name "*.csv" \) | head -n1 || true)"
    if [[ -n "$counts_found" ]]; then
      cp -f "$counts_found" "${out_dir}/miRNA.counts.tsv"
    fi
  fi

  echo "[$(date +'%Y-%m-%d %H:%M:%S')] [done] ${base} cuantificado. Salida: ${out_dir}" | tee -a "$main_log"
}

export -f process_one

# Ejecutar en serie (robusto y con logs claros). Si prefieres paralelo, podemos añadir -P.
for f in "${TRIM_FILES[@]}"; do
  process_one "$f" "$OUT_DIR" "$SPECIES" "$ADAPTER" "$THREADS" "$LOGFILE"
done

# Resumen final
TOTAL_DIRS=$(find "$OUT_DIR" -maxdepth 1 -type d ! -path "$OUT_DIR" | wc -l | tr -d ' ')
TOTAL_COUNTS=$(find "$OUT_DIR" -type f -name "miRNA.counts.tsv" | wc -l | tr -d ' ')
log "Muestras procesadas: $TOTAL_DIRS | Archivos de conteo estandarizados: $TOTAL_COUNTS"

# --- Fin ---
log "Script finalizado correctamente."
