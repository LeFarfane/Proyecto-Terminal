#!/usr/bin/env bash
# auto_bioproject_runinfo.sh
# Usage: ./auto_bioproject_runinfo.sh bioproject_ids.txt

set -uo pipefail

LIST_FILE="${1:-}"
if [[ -z "$LIST_FILE" || ! -f "$LIST_FILE" ]]; then
  echo "Usage: $0 bioproject_ids.txt"
  exit 1
fi

# Open the list file on FD 3 so stdin stays free
exec 3< "$LIST_FILE"

while IFS= read -r BIOPROJECT_ID <&3; do
  [[ -z "$BIOPROJECT_ID" ]] && continue
  [[ "$BIOPROJECT_ID" =~ ^# ]] && continue

  BIOPROJECT_ID="$(echo "$BIOPROJECT_ID" | xargs)"

  echo "==> Processing $BIOPROJECT_ID ..."

  mkdir -p "$BIOPROJECT_ID"

  if ! (
    cd "$BIOPROJECT_ID"
    # Also force these tools to not touch the loop’s stdin
    esearch -db bioproject -query "$BIOPROJECT_ID" < /dev/null \
      | elink -target sra \
      | efetch -format runinfo \
      > "${BIOPROJECT_ID}.csv"
  ); then
    echo "!! Failed: $BIOPROJECT_ID" | tee -a failed_ids.log
    continue
  fi

  if [[ ! -s "$BIOPROJECT_ID/$BIOPROJECT_ID.csv" ]]; then
    echo "!! Empty runinfo (no public runs?): $BIOPROJECT_ID" | tee -a empty_runinfo.log
    rm -f "$BIOPROJECT_ID/$BIOPROJECT_ID.csv"
    continue
  fi

  echo "OK: $BIOPROJECT_ID/$BIOPROJECT_ID.csv"
  sleep 0.34
done

# Close FD 3
exec 3<&-
