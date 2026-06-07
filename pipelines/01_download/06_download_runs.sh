#!/usr/bin/env bash
# ===========================================
# Script: 01_download/download_runs.sh
# Descripción: Descarga lecturas crudas de SRA en un solo directorio local.
#              Informa progreso en terminal y mantiene logs locales.
# Fecha: $(date +%Y-%m-%d)
# ===========================================

# --- Configuración segura ---
set -euo pipefail
IFS=$'\n\t'

# --- Variables ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"             

# Usar el directorio actual de trabajo (donde el usuario ejecuta el comando)
CURRENT_WDIR="$(pwd)"

LOGFILE="${CURRENT_WDIR}/01_download.log"
export LOGFILE  

# Si no se pasa argumento, busca srr_list.txt en la carpeta actual
INPUT_RUNS="${1:-${CURRENT_WDIR}/srr_list.txt}"
# La salida será una carpeta 'sra_fastq' en el directorio actual
OUTPUT_DIR="${2:-${CURRENT_WDIR}/sra_fastq}"

FQ_THREADS="${3:-4}"
PAR_JOBS="${4:-3}"

# --- Funciones ---
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >> "$LOGFILE" ; }
export -f log

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "❌ ERROR: comando requerido no encontrado: $1"; exit 127; }; }

collect_srrs() {
  local in="$1"
  [[ -s "$in" ]] || { echo "❌ ERROR: archivo de entrada vacío o inexistente: $in"; exit 1; }
  if head -1 "$in" | grep -qi '^Run[, ]'; then
    awk -F',' 'NR==1{for(i=1;i<=NF;i++){if($i=="Run") c=i}} NR>1 && c{print $c}' "$in" \
      | sed 's/\r$//' | grep -E '^SRR[0-9]+' | sort -u
  else
    sed 's/\r$//' "$in" | grep -E '^SRR[0-9]+' | sort -u
  fi
}

download_one() {
  local SRR="$1" outdir="$2" threads="$3"

  # Lógica de contador para la barra de progreso
  echo "" >> "$TMP_COUNT"
  local CURRENT=$(wc -l < "$TMP_COUNT" | tr -d ' ')
  local PREFIX="[$CURRENT/$TOTAL_SRR] $SRR"

  # 1. Ya existe comprimido
  if ls "${outdir}/${SRR}"*.fastq.gz >/dev/null 2>&1; then
    echo "🟢 $PREFIX -> Ya existe (.fastq.gz). Omitiendo."
    log "$PREFIX -> skip"
    return 0
  fi

  # 2. Existe sin comprimir → comprimir
  if ls "${outdir}/${SRR}"*.fastq >/dev/null 2>&1; then
    echo "🟡 $PREFIX -> FASTQ encontrado sin comprimir. Comprimiendo a .gz..."
    if command -v pigz >/dev/null 2>&1; then pigz -f -p "${threads}" "${outdir}/${SRR}"*.fastq
    else gzip -f "${outdir}/${SRR}"*.fastq; fi
    echo "✅ $PREFIX -> Comprimido con éxito."
    log "$PREFIX -> compressed"
    return 0
  fi

  # 3. Descargar nuevo
  echo "🔵 $PREFIX -> Descargando..."
  log "$PREFIX -> fetch"

  # Carpeta temporal oculta para evitar subcarpetas
  local SRA_CACHE="${outdir}/.sra_cache"
  mkdir -p "$SRA_CACHE"

  if ! prefetch "${SRR}" -O "$SRA_CACHE" >>"$LOGFILE" 2>&1; then
    log "$PREFIX -> warn: prefetch falló, intentando dump directo"
    fasterq-dump "${SRR}" -O "${outdir}" -e "${threads}" --split-files >>"$LOGFILE" 2>&1
  else
    fasterq-dump "$SRA_CACHE/${SRR}/${SRR}.sra" -O "${outdir}" -e "${threads}" --split-files >>"$LOGFILE" 2>&1
    rm -rf "$SRA_CACHE/${SRR}" 
  fi

  echo "🟣 $PREFIX -> Descargado. Comprimiendo..."
  
  if command -v pigz >/dev/null 2>&1; then pigz -p "${threads}" "${outdir}/${SRR}"*.fastq 2>/dev/null || true
  else gzip "${outdir}/${SRR}"*.fastq 2>/dev/null || true; fi

  echo "✅ $PREFIX -> Listo y guardado."
  log "$PREFIX -> done"
}

export -f download_one

# --- Inicio ---
mkdir -p "$OUTPUT_DIR"
: > "$LOGFILE"

echo "================================================"
echo " Iniciando descarga de lecturas SRA"
echo "================================================"
log "Script iniciado."

for c in prefetch fasterq-dump awk sort xargs; do need_cmd "$c"; done

# Mantenemos el caché global de NCBI en la raíz para no ensuciar las carpetas GSE
export NCBI_SETTINGS="${NCBI_SETTINGS:-${ROOT_DIR}/.ncbi}"
mkdir -p "$NCBI_SETTINGS"

TMP_SRR="$(mktemp)"
TMP_COUNT="$(mktemp)"
export TMP_COUNT

collect_srrs "$INPUT_RUNS" > "$TMP_SRR"
TOTAL_SRR="$(wc -l < "$TMP_SRR" | tr -d ' ')"
export TOTAL_SRR

[[ "$TOTAL_SRR" -gt 0 ]] || { echo "❌ ERROR: no se encontraron SRR en la entrada."; rm -f "$TMP_SRR" "$TMP_COUNT"; exit 1; }

echo "📁 Salida local: $OUTPUT_DIR"
echo "📊 Total de SRR a procesar: $TOTAL_SRR"
echo "⚙️  Paralelismo: $PAR_JOBS descargas simultáneas"
echo "------------------------------------------------"

# Descargar en paralelo
xargs -I{} -P "${PAR_JOBS}" bash -c 'download_one "$@"' _ {} "$OUTPUT_DIR" "$FQ_THREADS" < "$TMP_SRR"

# Limpieza temporal
rm -f "$TMP_SRR" "$TMP_COUNT"

# Verificación final
echo "------------------------------------------------"
DL_COUNT=$(ls -1 "${OUTPUT_DIR}"/SRR*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
if [[ "$DL_COUNT" -gt 0 ]]; then
  echo "🎉 Proceso finalizado. Archivos .fastq.gz generados: $DL_COUNT"
else
  echo "⚠️  ADVERTENCIA: No se generaron archivos .fastq.gz. Revisa el log local: $LOGFILE"
fi

echo "================================================"