
#!/usr/bin/env bash
# run_mirge3_per_run_single.sh
# Usage:
#   ./run_mirge3_per_run_single.sh runs.txt
#
# Controls:
#   THREADS=6    -> threads used INSIDE miRge3 (-cpu)
#   PARALLEL=1   -> how many runs to process at the same time
#   FASTQ_DIR=.  -> where to look for files when you provide only SRR IDs

set -euo pipefail

LIST="${1:-}"
OUTBASE="mirge3_output"
THREADS="${THREADS:-6}"
PARALLEL="${PARALLEL:-1}"
FASTQ_DIR="${FASTQ_DIR:-.}"

# Resolve the library-prep kit -> miRge adapter/UMI params (see ../kits.tsv).
# Set KIT=... explicitly, or DATASET=<id> to look it up. Defaults to illumina.
#   KIT=nextflex_v3 ./10_single_folder_mirge3.sh runs.txt
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/kit_params.sh
source "$SCRIPT_DIR/../lib/kit_params.sh"
KIT="$(resolve_kit)"
load_kit_params "$KIT"
echo "Kit: $KIT  ->  miRge -a $MIRGE_ADAPTER ${MIRGE_UMI_ARGS[*]:-(no UMI)}" >&2

if [[ -z "$LIST" || ! -f "$LIST" ]]; then
  echo "Usage: $0 runs.txt" >&2
  exit 1
fi

mkdir -p "$OUTBASE"

trim() { echo "$1" | awk '{$1=$1;print}'; }

resolve_file() {
  local x="$1"

  if [[ -f "$x" ]]; then
    echo "$x"
    return 0
  fi

  if [[ "$x" == *.fastq.gz ]]; then
    [[ -f "$FASTQ_DIR/$x" ]] && { echo "$FASTQ_DIR/$x"; return 0; }
  else
    [[ -f "$FASTQ_DIR/$x.fastq.gz" ]] && { echo "$FASTQ_DIR/$x.fastq.gz"; return 0; }
  fi

  return 1
}

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

run_one() {
  local file="$1"
  local run_id="$2"
  local run_dir="$3"

  mkdir -p "$run_dir"

  echo "==> Running miRge3 for: $run_id" >&2
  echo "    File: $file" >&2
  echo "    Out : $run_dir" >&2
  echo "    Log : $run_dir/mirge3.log" >&2
  echo >&2

  # Stream output to terminal AND log file, line-buffered, with timestamps
  stdbuf -oL -eL miRge3.0 \
    -s "$file" \
    -db miRBase \
    -lib "${MIRGE3_LIB:-${CONDA_PREFIX:-/home/genesis/miniconda3/envs/py_env}/miRge3_Lib}" \
    -on human \
    -ex 0.1 \
    -ie \
    -cpu "$THREADS" \
    -o "$run_dir" \
    -spl \
    -a "$MIRGE_ADAPTER" \
    ${MIRGE_UMI_ARGS[@]+"${MIRGE_UMI_ARGS[@]}"} \
    -nxt 20 \
    -q 20 \
    -NX \
    -m 18 \
    -mEC \
    -gff \
    
    2>&1 | awk -v rid="$run_id" '{ print strftime("[%Y-%m-%d %H:%M:%S]"), "["rid"]", $0 }' | tee "$run_dir/mirge3.log"
}

run_count=0

while IFS= read -r line <&3; do
  line="$(trim "$line")"
  [[ -z "$line" || "$line" == \#* ]] && continue

  if ! file="$(resolve_file "$line")"; then
    echo "ERROR: Could not resolve file for line: $line" >&2
    echo "       Tried: '$line', '$FASTQ_DIR/$line', '$FASTQ_DIR/$line.fastq.gz'" >&2
    exit 1
  fi

  [[ -f "$file" ]] || { echo "ERROR: file not found: $file" >&2; exit 1; }

  run_id="$(get_run_id "$file")"
  run_dir="$OUTBASE/$run_id"

  while (( $(jobs -pr | wc -l) >= PARALLEL )); do
    wait -n
  done

  run_one "$file" "$run_id" "$run_dir" &
  run_count=$((run_count + 1))
done 3< "$LIST"

wait
echo "Done. Runs processed: $run_count" >&2
