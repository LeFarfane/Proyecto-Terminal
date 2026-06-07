#!/usr/bin/env bash
# run_mirge_safe_sequential.sh
# Safely iterates through child directories and processes samples one by one.

set -euo pipefail

# Auto-navigate to sra_fastq/muestras_clasificadas from the dataset root
TARGET_DIR="sra_fastq/muestras_clasificadas"
if [[ ! -d "$TARGET_DIR" ]]; then
    echo "❌ No se encontró $TARGET_DIR en: $(pwd)"
    exit 1
fi

# Warn (don't block): /mnt/ (Windows-mounted drvfs) is 10-50x slower for miRge's
# small-I/O workload. Prefer the native Linux fs (e.g. ~/bioprojects/), but allow
# the user to proceed at their own risk.
if [[ "$(pwd)" == /mnt/* ]]; then
    echo "⚠️  WARNING: Running from a Windows-mounted drive (/mnt/...) — miRge will be"
    echo "   much slower here. For best speed copy the dataset to the native Linux"
    echo "   filesystem (e.g. ~/bioprojects/) and run from there. Continuing anyway..."
fi
cd "$TARGET_DIR"
echo "📂 Trabajando en: $(pwd)"

THREADS=6

# Resolve the library-prep kit -> miRge adapter/UMI params (see ../kits.tsv).
# Set KIT=... explicitly, or DATASET=<id> to look it up. Defaults to illumina.
#   KIT=nextflex_v3 ./09_auto_subfolders_mirge3.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/kit_params.sh
source "$SCRIPT_DIR/../lib/kit_params.sh"
KIT="$(resolve_kit)"
load_kit_params "$KIT"
echo "Kit: $KIT  ->  miRge -a $MIRGE_ADAPTER ${MIRGE_UMI_ARGS[*]:-(no UMI)}"

# Function from your script to cleanly extract the sample name
get_run_id() {
  local f="$1"
  local b
  b="$(basename "$f")"
  b="${b%.fastq.gz}"
  b="${b%.fastq}"
  b="${b%.fq.gz}"
  b="${b%.fq}"
  echo "$b"
}

echo "=================================================="
echo "Scan mode selection:"
echo "  [1] Full scan   — process ALL .fastq.gz files found in each folder"
echo "  [2] List mode   — process only files listed in fasta_list.txt (if present)"
echo "=================================================="
read -rp "Enter choice [1/2]: " SCAN_CHOICE
case "$SCAN_CHOICE" in
  2) SCAN_MODE="list"  ; echo "Mode: list (fasta_list.txt)" ;;
  *) SCAN_MODE="full"  ; echo "Mode: full scan" ;;
esac

echo "Starting nested sequential processing..."

# OUTER LOOP: Iterate through all child directories
for child_dir in */; do
  # Verify it is a valid directory
  [[ -d "$child_dir" ]] || continue

  # Clean up the directory name by removing the trailing slash
  clean_dir="${child_dir%/}"

  echo "=================================================="
  echo "==> Entering folder: $clean_dir"
  echo "=================================================="

  # Use a subshell to jump into the folder safely
  (
    cd "$child_dir" || exit 1

    # Build the file list according to the chosen scan mode
    if [[ "$SCAN_MODE" == "list" ]]; then
      if [[ ! -f fasta_list.txt ]]; then
        echo "    No fasta_list.txt found in $clean_dir. Skipping..."
        exit 0
      fi
      mapfile -t fastq_files < fasta_list.txt
      # Drop blank lines and entries that don't actually exist
      valid_files=()
      for f in "${fastq_files[@]}"; do
        [[ -z "$f" ]] && continue
        if [[ -f "$f" ]]; then
          valid_files+=("$f")
        else
          echo "    WARNING: $f listed in fasta_list.txt but not found. Skipping."
        fi
      done
      fastq_files=("${valid_files[@]}")
    else
      shopt -s nullglob
      fastq_files=(*.fastq.gz)
      shopt -u nullglob
    fi

    if [[ ${#fastq_files[@]} -eq 0 ]]; then
      echo "    No .fastq.gz files to process in $clean_dir. Skipping..."
      exit 0 # Exits this folder, outer loop moves to the next one
    fi

    # INNER LOOP: Process each sample ONE by ONE
    for file in "${fastq_files[@]}"; do
      
      run_id="$(get_run_id "$file")"
      
      # Create a dedicated output folder for THIS single sample
      out_dir="mirge3_output/$run_id"
      mkdir -p "$out_dir"

      echo "    ----------------------------------------------"
      echo "    ==> Analyzing sample: $run_id"
      echo "    File: $file"
      echo "    Out : $clean_dir/$out_dir"
      echo "    ----------------------------------------------"

      # Run miRge3.0 on just this one file, using your exact logic
      stdbuf -oL -eL miRge3.0 \
        -s "$file" \
        -db miRBase \
        -lib "${MIRGE3_LIB:-${CONDA_PREFIX:-/home/genesis/miniconda3/envs/py_env}/miRge3_Lib}" \
        -on human \
        -ex 0.1 \
        -ie \
        -cpu "$THREADS" \
        -o "$out_dir" \
        -spl \
        -a "$MIRGE_ADAPTER" \
        ${MIRGE_UMI_ARGS[@]+"${MIRGE_UMI_ARGS[@]}"} \
        -nxt 20 \
        -q 20 \
        -NX \
        -m 18 \
        -gff \
        2>&1 | awk -v rid="$run_id" '{ print strftime("[%Y-%m-%d %H:%M:%S]"), "["rid"]", $0 }' | tee "$out_dir/mirge3.log"
      
      echo "    ==> Finished sample: $run_id"
    done
    
    echo "==> Finished all samples in folder: $clean_dir"
  )
done

echo "Done! All directories and samples processed successfully."