#!/usr/bin/env bash
# kit_params.sh -- shared library-prep kit -> tool-parameter mapping.
#
# Source this file, then resolve the kit and its parameters:
#
#   source /path/to/pipelines/lib/kit_params.sh
#   KIT="$(resolve_kit)"            # from $KIT, else $DATASET via kits.tsv, else illumina
#   load_kit_params "$KIT"          # sets MIRGE_ADAPTER, MIRGE_UMI_ARGS[], MIRTRACE_ADAPTER, MIRTRACE_PROTOCOL
#
# Why a map and not auto-detection: the kit is a property of the study, and
# NEXTflex v3 shares the exact 3' adapter (TGGAATTCTCGGGTGCCAAGG) with standard
# Illumina -- they differ only by 4+4 random-mer UMIs, which adapter-sniffing
# cannot see. So the kit is assigned per dataset in kits.tsv.

# Directory of this script (works regardless of caller's cwd).
_KIT_PARAMS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KITS_TSV="${KITS_TSV:-${_KIT_PARAMS_DIR}/../kits.tsv}"

# Look up the kit for a dataset id from kits.tsv. Defaults to "illumina".
kit_for_dataset() {
  local dataset="$1"
  local kit=""
  if [[ -n "$dataset" && -f "$KITS_TSV" ]]; then
    # Match first column exactly, ignore comments/blank lines.
    kit="$(awk -F'\t' -v d="$dataset" '
      /^[[:space:]]*#/ {next} NF<2 {next}
      $1==d {print $2; exit}' "$KITS_TSV")"
  fi
  echo "${kit:-illumina}"
}

# Resolve the kit from the environment:
#   $KIT (explicit) takes priority; else look up $DATASET in kits.tsv; else illumina.
resolve_kit() {
  if [[ -n "${KIT:-}" ]]; then
    echo "$KIT"
  else
    kit_for_dataset "${DATASET:-}"
  fi
}

# Given a kit name, set:
#   MIRGE_ADAPTER      value for miRge3.0 -a
#   MIRGE_UMI_ARGS     bash array, e.g. (-umi 4,4) or empty
#   MIRTRACE_ADAPTER   constant 3' adapter sequence for the miRTrace config
#   MIRTRACE_PROTOCOL  miRTrace -p value: illumina | nextflex | qiaseq | cats
load_kit_params() {
  local kit="$1"
  MIRGE_UMI_ARGS=()
  case "$kit" in
    nextflex_v3)
      # Same Illumina small-RNA 3' adapter, but trim the 4+4 random-mers.
      MIRGE_ADAPTER="illumina"
      MIRGE_UMI_ARGS=(-umi 4,4)
      MIRTRACE_ADAPTER="TGGAATTCTCGGGTGCCAAGG"
      MIRTRACE_PROTOCOL="nextflex"
      ;;
    nebnext)
      MIRGE_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
      MIRTRACE_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
      MIRTRACE_PROTOCOL="illumina"
      ;;
    illumina|"")
      MIRGE_ADAPTER="illumina"
      MIRTRACE_ADAPTER="TGGAATTCTCGGGTGCCAAGG"
      MIRTRACE_PROTOCOL="illumina"
      ;;
    *)
      echo "WARNING: unknown kit '$kit' -- falling back to 'illumina'" >&2
      MIRGE_ADAPTER="illumina"
      MIRTRACE_ADAPTER="TGGAATTCTCGGGTGCCAAGG"
      MIRTRACE_PROTOCOL="illumina"
      ;;
  esac
}
