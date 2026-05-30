#!/usr/bin/env bash
# download_runs.sh
# Usage:
#   ./download_runs.sh runs.txt
#   ./download_runs.sh runs.txt fastq

set -euo pipefail

LIST="${1:-}"
OUTDIR="${2:-fastq}"
LOGFILE="download_failures.log"

if [[ -z "$LIST" || ! -f "$LIST" ]]; then
  echo "Usage: $0 runs.txt [outdir=fastq]" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

# Clear previous log file if it exists so we start fresh
: > "$LOGFILE"

# Read list on FD 3 so nothing else can accidentally consume stdin
exec 3< "$LIST"
while IFS= read -r SRR <&3; do
  # trim whitespace + remove possible Windows CR
  SRR="$(echo "$SRR" | tr -d '\r' | xargs)"

  # skip blank lines / comments
  [[ -z "$SRR" || "$SRR" == \#* ]] && continue

  echo "==> Downloading $SRR"

  # We use an 'if' statement here. 
  # This prevents 'set -e' from crashing the script if fasterq-dump fails.
  if fasterq-dump "$SRR" -O "$OUTDIR" --split-files --threads 6 --progress; then
      echo "✅ Success: $SRR"
  else
      echo "❌ Error: Failed to download $SRR. Logging and skipping..." >&2
      echo "$SRR" >> "$LOGFILE"
  fi

  # polite pause (optional)
  sleep 0.25
done
exec 3<&-

# Final report
if [[ -s "$LOGFILE" ]]; then
    echo "------------------------------------------------"
    echo "⚠️  Some downloads failed. See $LOGFILE for details."
    echo "   You can retry them by running: ./download_runs.sh $LOGFILE"
else
    echo "------------------------------------------------"
    echo "🎉 All downloads completed successfully!"
fi
