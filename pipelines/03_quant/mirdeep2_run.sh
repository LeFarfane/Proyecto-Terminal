#!/usr/bin/env bash
# ===========================================
# Script: 03_quant/mirdeep2_run.sh
# Descripción: Cuantifica miRNAs por muestra usando miRDeep2 (mapper.pl + miRDeep2.pl) a partir
#              de FASTQ recortados. Construye índice de Bowtie si no existe, mapea lecturas,
#              y genera salidas por muestra en results/quant_mirdeep2/<sample>/.
# Fecha: $(date +%Y-%m-%d)
# ===========================================

# --- Configuración segura ---
set -euo pipefail
IFS=$'\n\t'

# --- Variables ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"                               # .../PT
LOG_DIR="${ROOT_DIR}/outputs/logs"
LOGFILE="${LOG_DIR}/03_quant_mirdeep2.log"

# Parámetros (puedes sobreescribir por CLI):
# $1 = TRIM_DIR (dir con *.trim.fastq.gz)            [default: PT/data_raw/trimmed]
# $2 = OUT_DIR  (dir de salida por muestra)          [default: PT/results/quant_mirdeep2]
# $3 = GENOME_FA (FASTA del genoma)                  [default: PT/data_ref/human_ds/genome.fa]
# $4 = BOWTIE_INDEX_PREFIX (opcional; si no, se crea) [default: PT/data_ref/bowtie_index/human/genome]
# $5 = MAT_FA (miRBase mature.fa de hsa)             [default: PT/data_ref/mirbase/mature_hsa.fa]
# $6 = HAIRPIN_FA (miRBase hairpin.fa de hsa)        [default: PT/data_ref/mirbase/hairpin_hsa.fa]
# $7 = SPECIES code (-t)                              [default: hsa]
# $8 = MINLEN lectura colapsada (-l mapper.pl)       [default: 15]
TRIM_DIR="${1:-${ROOT_DIR}/data_raw/trimmed}"
OUT_DIR="${2:-${ROOT_DIR}/results/quant_mirdeep2}"
GENOME_FA_DEFAULT="${ROOT_DIR}/data_ref/human_ds/genome.fa"
GENOME_FA="${3:-$GENOME_FA_DEFAULT}"
BOWTIE_PREFIX_DEFAULT="${ROOT_DIR}/data_ref/bowtie_index/human/genome"
BOWTIE_INDEX_PREFIX="${4:-$BOWTIE_PREFIX_DEFAULT}"
MAT_FA="${5:-${ROOT_DIR}/data_ref/mirbase/mature_hsa.fa}"
HAIRPIN_FA="${6:-${ROOT_DIR}/data_ref/mirbase/hairpin_hsa.fa}"
SPECIES="${7:-hsa}"
MINLEN="${8:-15}"

# --- Funciones ---
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE" ; }

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    log "ERROR: comando requerido no encontrado: $1"
    exit 127
  fi
}

ensure_bowtie_index() {
  local fa="$1" prefix="$2"
  # bowtie index requiere archivos *.ebwt
  if ls "${prefix}".*.ebwt >/dev/null 2>&1; then
    log "[index] Índice Bowtie encontrado: ${prefix}*.ebwt"
    return 0
  fi
  log "[index] No existe índice Bowtie en ${prefix} → construyendo con bowtie-build ..."
  mkdir -p "$(dirname "$prefix")"
  bowtie-build "$fa" "$prefix" >>"$LOGFILE" 2>&1
  log "[index] Índice construido."
}

decompress_if_needed() {
  # echo path al FASTQ descomprimido (si entrada venía .gz, lo descomprime a una copia temporal en el dir de muestra)
  local in="$1" outdir="$2"
  if [[ "$in" == *.gz ]]; then
    local tmp="${outdir}/$(basename "${in%.gz}")"
    if [[ -s "$tmp" ]]; then
      echo "$tmp"
      return 0
    fi
    if command -v pigz >/dev/null 2>&1; then
      pigz -dc "$in" > "$tmp"
    else
      gzip -dc "$in" > "$tmp"
    fi
    echo "$tmp"
  else
    echo "$in"
  fi
}

# --- Inicio ---
mkdir -p "$LOG_DIR" "$OUT_DIR"
: > "$LOGFILE"
log "Iniciando script..."
log "TRIM_DIR            = $TRIM_DIR"
log "OUT_DIR             = $OUT_DIR"
log "GENOME_FA           = $GENOME_FA"
log "BOWTIE_INDEX_PREFIX = $BOWTIE_INDEX_PREFIX"
log "MAT_FA              = $MAT_FA"
log "HAIRPIN_FA          = $HAIRPIN_FA"
log "SPECIES (-t)        = $SPECIES"
log "MINLEN (-l)         = $MINLEN"

# Validar binarios requeridos
for c in mapper.pl miRDeep2.pl bowtie bowtie-build awk find; do need_cmd "$c"; done
command -v pigz >/dev/null 2>&1 || log "[info] pigz no encontrado; usaré gzip para descompresión"

# Validar insumos
if [ -d "$TRIM_DIR" ]; then
  log "Directorio de entrada encontrado: $TRIM_DIR"
else
  log "ERROR: Directorio de entrada no existe"
  exit 1
fi

[[ -s "$GENOME_FA" ]] || { log "ERROR: GENOME_FA no existe o está vacío: $GENOME_FA"; exit 1; }
[[ -s "$MAT_FA"    ]] || { log "ERROR: MAT_FA no existe o está vacío: $MAT_FA"; exit 1; }
[[ -s "$HAIRPIN_FA" ]] || { log "ERROR: HAIRPIN_FA no existe o está vacío: $HAIRPIN_FA"; exit 1; }

# Asegurar índice Bowtie
ensure_bowtie_index "$GENOME_FA" "$BOWTIE_INDEX_PREFIX"

# Listar archivos .trim.fastq(.gz)
mapfile -t TRIM_FILES < <(find "$TRIM_DIR" -maxdepth 1 -type f \( -name "*.trim.fastq.gz" -o -name "*.trim.fq.gz" -o -name "*.trim.fastq" -o -name "*.trim.fq" \) | sort)
TOTAL="${#TRIM_FILES[@]}"
if [[ "$TOTAL" -eq 0 ]]; then
  log "ERROR: No se encontraron archivos *.trim.fastq* en $TRIM_DIR"
  exit 1
fi
log "Archivos detectados para cuantificación (miRDeep2): $TOTAL"

# Procesar cada muestra
for in_f in "${TRIM_FILES[@]}"; do
  bn="$(basename "$in_f")"
  base="${bn%.trim.fastq.gz}"; base="${base%.trim.fq.gz}"; base="${base%.trim.fastq}"; base="${base%.trim.fq}"
  # Normaliza si trae sufijo de read
  base="${base%_1}"; base="${base%_R1}"

  smp_dir="${OUT_DIR}/${base}"
  mkdir -p "$smp_dir"
  log "[run] Muestra ${base}"

  # Descomprimir si es necesario
  fq_in="$(decompress_if_needed "$in_f" "$smp_dir")"
  log "[io] FASTQ para mapper.pl: $fq_in"

  # Paths de salida
  collapsed="${smp_dir}/reads_collapsed.fa"
  arf="${smp_dir}/reads_vs_genome.arf"
  smp_log="${smp_dir}/${base}.mirdeep2.log"

  # mapper.pl — colapsar y mapear con Bowtie (índice ya construido)
  # flags clave:
  #  -e: no mantener secuencias con Ns
  #  -h: eliminar adaptadores residuales de hairpin si los encuentra
  #  -m: colapsar secuencias idénticas
  #  -j: SALIDA ARF de alineamientos
  #  -l: longitud mínima
  #  -p: prefijo índice Bowtie
  #  -q: input es FASTQ (no FASTA)
  log "[mapper] Ejecutando mapper.pl (minlen=$MINLEN) ..."
  if ! mapper.pl "$fq_in" -e -h -m -j -l "$MINLEN" -q \
        -s "$collapsed" \
        -t "$arf" \
        -p "$BOWTIE_INDEX_PREFIX" \
        >>"$smp_log" 2>&1; then
    log "ERROR: mapper.pl falló para ${base} (ver $smp_log)"
    continue
  fi
  [[ -s "$collapsed" && -s "$arf" ]] || { log "ERROR: mapper.pl no generó archivos esperados para ${base}"; continue; }

  # miRDeep2.pl — cuantificación y detección
  log "[miRDeep2] Ejecutando miRDeep2.pl ..."
  if ! miRDeep2.pl "$collapsed" "$GENOME_FA" "$arf" "$MAT_FA" none "$HAIRPIN_FA" \
        -t "$SPECIES" -P -v >>"$smp_log" 2>&1; then
    log "ERROR: miRDeep2.pl falló para ${base} (ver $smp_log)"
    continue
  fi

  # Estandarizar archivo de conteos (buscar un CSV/TSV con expresión de miRNAs)
  # miRDeep2 genera típicamente archivos como: miRNAs_expressed_all_samples_*.csv
  counts_src="$(find "$smp_dir" -maxdepth 1 -type f \( -iname "miRNAs_expressed_*.*" -o -iname "*reads*per*miRNA*.*" -o -iname "*_result*.csv" \) | head -n1 || true)"
  if [[ -n "${counts_src}" ]]; then
    # Convertir CSV→TSV si es .csv (simple: sustituye coma por tab)
    if [[ "$counts_src" == *.csv ]]; then
      awk 'BEGIN{FS=",";OFS="\t"}{print $0}' "$counts_src" > "${smp_dir}/miRNA.counts.tsv" || true
    else
      cp -f "$counts_src" "${smp_dir}/miRNA.counts.tsv" || true
    fi
    log "[out] Conteos estandarizados: ${smp_dir}/miRNA.counts.tsv"
  else
    log "WARN: No se encontró archivo de conteos claro en ${smp_dir}; revisa el log y salidas de miRDeep2."
  fi

  # Limpieza opcional del FASTQ temporal
  if [[ "$fq_in" != "$in_f" && -f "$fq_in" ]]; then
    rm -f "$fq_in"
  fi

  log "[done] ${base} procesada."
done

# Resumen final
TOTAL_DIRS=$(find "$OUT_DIR" -maxdepth 1 -type d ! -path "$OUT_DIR" | wc -l | tr -d ' ')
TOTAL_COUNTS=$(find "$OUT_DIR" -type f -name "miRNA.counts.tsv" | wc -l | tr -d ' ')
log "Muestras procesadas: $TOTAL_DIRS | Archivos de conteo estandarizados: $TOTAL_COUNTS"

# --- Fin ---
log "Script finalizado correctamente."
