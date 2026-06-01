#!/bin/bash
set -euo pipefail

CONFIG="${1:-mirtrace_config.csv}"
OUTDIR="mirtrace_results"
BATCH_SIZE=10
unset _JAVA_OPTIONS

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/kit_params.sh"
KIT="$(resolve_kit)"
load_kit_params "$KIT"
echo "Kit: $KIT  ->  miRTrace protocol: $MIRTRACE_PROTOCOL"

# Find the miRTrace JAR and its bundled java (bypasses the Python wrapper which
# miscalculates memory on WSL, and avoids relying on java being in PATH).
# Search from the real home dir (not ~ which may differ under wsl bash -lc).
_SEARCH_ROOT="${HOME:-$(eval echo ~)}"
MIRTRACE_JAR=$(find "$_SEARCH_ROOT/miniconda3/envs/mirtrace" -name "mirtrace.jar" 2>/dev/null | head -1)
if [ -z "$MIRTRACE_JAR" ]; then
    echo "ERROR: Could not find mirtrace.jar under $_SEARCH_ROOT/miniconda3/envs/mirtrace"; exit 1
fi
JAVA=$(find "$_SEARCH_ROOT/miniconda3/envs/mirtrace" -name "java" -type f 2>/dev/null | head -1)
JAVA="${JAVA:-java}"
echo "Using JAR: $MIRTRACE_JAR"
echo "Using java: $JAVA"

# Split config into batches and run each
total=$(wc -l < "$CONFIG")
batch=1
start=1

while [ $start -le $total ]; do
    end=$((start + BATCH_SIZE - 1))
    batch_config="mirtrace_config_batch${batch}.csv"
    batch_out="${OUTDIR}_batch${batch}"

    sed -n "${start},${end}p" "$CONFIG" > "$batch_config"

    echo "=== Running batch $batch (lines $start–$end) ==="
    "$JAVA" -Xms12g -Xmx12g -jar "$MIRTRACE_JAR" \
        --mirtrace-wrapper-name mirtrace \
        qc --species hsa -p "$MIRTRACE_PROTOCOL" -c "$batch_config" --output-dir "$batch_out" --force

    start=$((end + 1))
    batch=$((batch + 1))
done

echo ""
echo "=== All done. Results in: ${OUTDIR}_batch1 through _batch${batch} ==="
