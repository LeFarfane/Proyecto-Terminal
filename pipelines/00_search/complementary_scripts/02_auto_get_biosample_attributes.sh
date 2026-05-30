#!/usr/bin/env bash
# biosample_ids_to_dynamic_csv.sh
# Usage:
#   ./biosample_ids_to_dynamic_csv.sh biosamples.txt > out.csv
#   ./biosample_ids_to_dynamic_csv.sh biosamples.txt out.csv
#
# Input:
#   biosamples.txt = one BioSample ID per line (SAMN...). Blank lines / comments allowed.

set -euo pipefail

LIST="${1:-}"
OUT="${2:-}"

if [[ -z "$LIST" || ! -f "$LIST" ]]; then
  echo "Usage: $0 biosamples.txt [out.csv]" >&2
  exit 1
fi

# If OUT is provided, write to it; otherwise stdout.
if [[ -n "$OUT" ]]; then
  exec > "$OUT"
fi

# CSV escaping: wrap in quotes, double any internal quotes
csvq() {
  local s="${1:-}"
  s="${s//\"/\"\"}"
  printf '"%s"' "$s"
}

# Collect BioSample IDs (skip blanks/comments)
biosamples=()
while IFS= read -r line; do
  line="$(echo "$line" | awk '{$1=$1; print}')"
  [[ -z "$line" || "$line" == \#* ]] && continue
  biosamples+=("$line")
done < "$LIST"

if (( ${#biosamples[@]} == 0 )); then
  echo "No BioSample IDs found in $LIST" >&2
  exit 1
fi

# Fetch attributes for the first sample ONCE; use it to build headers
first_id="${biosamples[0]}"

first_kv="$(
  efetch -db biosample -id "$first_id" -format xml \
    | xtract -pattern Attribute -element "@attribute_name" -element "."
)"

# Build ordered, de-duplicated attribute list from first sample
declare -A seen=()
attrs=()
while IFS=$'\t' read -r attr_name attr_val; do
  [[ -z "${attr_name:-}" ]] && continue
  if [[ -z "${seen[$attr_name]+x}" ]]; then
    seen["$attr_name"]=1
    attrs+=("$attr_name")
  fi
done <<< "$first_kv"

# Output header: biosample_id + discovered attributes
# (quote headers for safety in case any header contains commas)
csvq "biosample_id"
for a in "${attrs[@]}"; do
  echo -n ","
  csvq "$a"
done
echo

# Helper: output one row given an ID and a kv block
emit_row_from_kv() {
  local samn="$1"
  local kv="$2"

  declare -A m=()
  while IFS=$'\t' read -r k v; do
    [[ -z "${k:-}" ]] && continue
    m["$k"]="$v"
  done <<< "$kv"

  csvq "$samn"
  for a in "${attrs[@]}"; do
    echo -n ","
    csvq "${m[$a]:-}"
  done
  echo
}

# Emit first row using cached attributes (no refetch)
emit_row_from_kv "$first_id" "$first_kv"

# Emit remaining rows (fetch once per biosample)
for (( i=1; i<${#biosamples[@]}; i++ )); do
  samn="${biosamples[i]}"
  kv="$(
    efetch -db biosample -id "$samn" -format xml \
      | xtract -pattern Attribute -element "@attribute_name" -element "."
  )"
  emit_row_from_kv "$samn" "$kv"
done
