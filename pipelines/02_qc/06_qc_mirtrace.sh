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
# Prefer the active conda env ($CONDA_PREFIX, i.e. py_env); fall back to the
# legacy py_env / mirtrace env locations so the script works during migration.
_SEARCH_ROOT="${HOME:-$(eval echo ~)}"
MIRTRACE_JAR=""
_FOUND_ROOT=""
for _root in "${CONDA_PREFIX:-}" "$_SEARCH_ROOT/miniconda3/envs/py_env" "$_SEARCH_ROOT/miniconda3/envs/mirtrace"; do
    [ -n "$_root" ] && [ -d "$_root" ] || continue
    MIRTRACE_JAR=$(find "$_root" -name "mirtrace.jar" 2>/dev/null | head -1)
    if [ -n "$MIRTRACE_JAR" ]; then _FOUND_ROOT="$_root"; break; fi
done
if [ -z "$MIRTRACE_JAR" ]; then
    echo "ERROR: Could not find mirtrace.jar (looked in CONDA_PREFIX, py_env, mirtrace envs)"; exit 1
fi
JAVA=$(find "$_FOUND_ROOT" -name "java" -type f 2>/dev/null | head -1)
JAVA="${JAVA:-java}"
echo "Using JAR: $MIRTRACE_JAR"
echo "Using java: $JAVA"

# ---------------------------------------------------------------------------
# Strategy: try one run over the whole config first (fast path → native combined
# report). If that fails (e.g. OOM on large datasets) or yields no report, fall
# back to batches of BATCH_SIZE + MultiQC merge + tidy-up.
# ---------------------------------------------------------------------------
echo ""
echo "=== Attempt 1: single run over all samples ==="
rm -rf "$OUTDIR"
single_ok=0
# Guard with 'if' so a crash doesn't trip 'set -e'; also require the HTML report.
if "$JAVA" -Xms12g -Xmx12g -jar "$MIRTRACE_JAR" \
        --mirtrace-wrapper-name mirtrace \
        qc --species hsa -p "$MIRTRACE_PROTOCOL" -c "$CONFIG" --output-dir "$OUTDIR" --force \
   && [ -f "$OUTDIR/mirtrace-report.html" ]; then
    single_ok=1
fi

if [ "$single_ok" -eq 1 ]; then
    echo ""
    echo "=== Done. Single run succeeded: $OUTDIR/mirtrace-report.html ==="
    exit 0
fi

echo ""
echo "!!! Single run failed or produced no report — falling back to batched mode !!!"
rm -rf "$OUTDIR"   # discard any partial single-run output before batching

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
echo "=== Combining batch reports with MultiQC ==="
# miRTrace can't merge runs; MultiQC's mirtrace module aggregates the per-batch
# stats into one report. Resolve multiqc the same way we resolve the JAR above
# (it lives in a conda env that may not be on PATH).
MULTIQC=$(command -v multiqc || true)
if [ -z "$MULTIQC" ]; then
    MULTIQC=$(find "$_SEARCH_ROOT/miniconda3/envs" -name "multiqc" -type f 2>/dev/null | head -1)
fi
if [ -n "$MULTIQC" ]; then
    "$MULTIQC" -f -m mirtrace -n mirtrace_combined_report.html "${OUTDIR}"_batch*
    echo "Combined report: mirtrace_combined_report.html"
else
    echo "WARNING: multiqc not found; skipping combined report. Install it or run manually:"
    echo "  multiqc -f -m mirtrace -n mirtrace_combined_report.html ${OUTDIR}_batch*"
fi

echo ""
echo "=== Tidying batch folders into ${OUTDIR}/ ==="
# Collect all per-batch result folders under one directory. Run this AFTER the
# MultiQC step above, which globs ${OUTDIR}_batch* in place.
mkdir -p "$OUTDIR"
mv "${OUTDIR}"_batch* "$OUTDIR"/ 2>/dev/null || true
mv mirtrace_config_batch*.csv "$OUTDIR"/ 2>/dev/null || true

echo ""
echo "=== All done. Per-batch results gathered in: ${OUTDIR}/ ==="
echo "=== Combined report: mirtrace_combined_report.html ==="
