#!/usr/bin/env bash
# =============================================================================
# 00_tool_check.sh — verify the full miRNA-seq pipeline toolchain
# =============================================================================
# This project uses TWO conda environments (see envs/*.yml and the README):
#
#   • py_env  — stages 00–04 (discovery, download, QC, quantification)
#               CLI / Java / Python tools.  Pinned to Python 3.8 by miRge3.
#   • r_env   — stages 05–08 (differential expression, enrichment, networks)
#               R 4.x + Bioconductor.
#
# Run from WSL / Linux:
#     bash pipelines/00_search/00_tool_check.sh
#
# Exits 0 if everything required is present; exits 1 and writes a log if not.
# Override env names with PY_ENV=... R_ENV=... if you renamed them.
# =============================================================================

set -uo pipefail   # deliberately NOT -e: collect ALL problems, don't bail early

PY_ENV="${PY_ENV:-py_env}"
R_ENV="${R_ENV:-r_env}"

# ── Required tooling ─────────────────────────────────────────────────────────
PY_CLI="esearch efetch elink esummary xtract prefetch fasterq-dump fastqc multiqc cutadapt bowtie miRge3.0 mirtrace seqkit java"
PY_CLI_OPT="pigz jq"                              # optional (graceful fallbacks)
PY_MODULES="pandas numpy matplotlib seaborn"
R_PKGS="tidyverse readr dplyr tidyr stringr ggplot2 DESeq2 apeglm umap pheatmap multiMiR clusterProfiler org.Hs.eg.db enrichplot pROC jsonlite httr fgsea igraph visNetwork ashr"

# ── Locate conda ─────────────────────────────────────────────────────────────
CONDA_BIN=""
for c in "$(command -v conda || true)" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
    [ -n "$c" ] && [ -x "$c" ] && { CONDA_BIN="$c"; break; }
done
if [ -z "$CONDA_BIN" ]; then
    echo "❌ conda not found. Install Miniconda, then: conda env create -f envs/py_env.yml && conda env create -f envs/r_env.yml"
    exit 1
fi

env_exists() { "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$1"; }
crun()       { "$CONDA_BIN" run -n "$1" "${@:2}"; }

MISSING=()
WARN=()

echo "============================================================"
echo " miRNA-seq pipeline dependency check"
echo "   py_env (stages 00–04)  +  r_env (stages 05–08)"
echo "============================================================"

# ── py_env ───────────────────────────────────────────────────────────────────
if ! env_exists "$PY_ENV"; then
    echo "❌ conda env '$PY_ENV' not found — create it: conda env create -f envs/py_env.yml"
    MISSING+=("env:$PY_ENV")
else
    echo ""
    echo "── $PY_ENV : CLI / Java tools ──────────────────────────────"
    for t in $PY_CLI; do
        if crun "$PY_ENV" bash -c "command -v '$t'" >/dev/null 2>&1; then
            echo "  ✅ $t"
        else
            echo "  ❌ $t"; MISSING+=("$PY_ENV:$t")
        fi
    done
    for t in $PY_CLI_OPT; do
        if crun "$PY_ENV" bash -c "command -v '$t'" >/dev/null 2>&1; then
            echo "  ✅ $t (optional)"
        else
            echo "  ⚠️  $t (optional — has a fallback)"; WARN+=("$PY_ENV:$t")
        fi
    done

    echo "── $PY_ENV : Python modules ────────────────────────────────"
    for m in $PY_MODULES; do
        if crun "$PY_ENV" python -c "import $m" >/dev/null 2>&1; then
            echo "  ✅ $m"
        else
            echo "  ❌ $m"; MISSING+=("$PY_ENV:py:$m")
        fi
    done

    if crun "$PY_ENV" bash -c '[ -d "$CONDA_PREFIX/miRge3_Lib" ]' >/dev/null 2>&1; then
        echo "  ✅ miRge3_Lib reference data (~11 GB)"
    else
        echo "  ⚠️  miRge3_Lib reference data missing under $PY_ENV"; WARN+=("$PY_ENV:miRge3_Lib")
    fi
fi

# ── r_env ────────────────────────────────────────────────────────────────────
if ! env_exists "$R_ENV"; then
    echo "❌ conda env '$R_ENV' not found — create it: conda env create -f envs/r_env.yml"
    MISSING+=("env:$R_ENV")
else
    echo ""
    echo "── $R_ENV : R / Bioconductor packages ──────────────────────"
    pkgs_csv=$(printf '"%s",' $R_PKGS); pkgs_csv="${pkgs_csv%,}"
    # Single Rscript call: print the names of any packages NOT installed.
    miss=$(crun "$R_ENV" Rscript -e "ip<-rownames(installed.packages()); want<-c($pkgs_csv); cat(want[!want %in% ip], sep='\n')" 2>/dev/null)
    if [ -z "$miss" ]; then
        echo "  ✅ all R packages present"
    else
        # echo present ones briefly, then missing
        while IFS= read -r p; do
            [ -n "$p" ] && { echo "  ❌ $p"; MISSING+=("$R_ENV:$p"); }
        done <<< "$miss"
        echo "  (remaining R packages OK)"
    fi
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
if [ "${#WARN[@]}" -gt 0 ]; then
    echo "⚠️  Optional / non-fatal: ${WARN[*]}"
fi
if [ "${#MISSING[@]}" -eq 0 ]; then
    echo "✅ All required tools present."
    echo "   Activate '$PY_ENV' for stages 00–04, '$R_ENV' for stages 05–08."
    echo "============================================================"
    exit 0
else
    LOG_DIR="tools needed log"
    LOG_FILE="$LOG_DIR/dependency_errors.log"
    mkdir -p "$LOG_DIR"
    {
        echo "--- Dependency Check Failed: $(date) ---"
        printf '  - %s\n' "${MISSING[@]}"
        echo "----------------------------------------"
    } >> "$LOG_FILE"
    echo "❌ Missing (${#MISSING[@]}): ${MISSING[*]}"
    echo "   Logged to: $LOG_FILE"
    echo "   Repair:  conda env update -n $PY_ENV -f envs/py_env.yml"
    echo "            conda env update -n $R_ENV  -f envs/r_env.yml"
    echo "============================================================"
    exit 1
fi
