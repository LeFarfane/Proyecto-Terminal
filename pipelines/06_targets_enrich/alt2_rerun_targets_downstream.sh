#!/usr/bin/env bash
# =============================================================================
# rerun_targets_downstream.sh
# Re-run the target/enrichment/IBD layers after removing the S-score gate from
# 16_multimir_targets.R (so targets now cover ALL DE-significant miRNAs:
# padj <= 0.05 & |log2FC| >= 0.58).
#
# Per dataset (outputs/<id>/) it re-runs, in order:
#   16_multimir_targets.R   -> multimir_outputs/      (targets for all DE-sig miRNAs)
#   18_pathway_enrich.R     -> pathway_outputs/       (GO/KEGG on those targets)
#   19_interactive_network.R-> network_outputs/       (miRNA-gene network)
#   22_ibd_target_overlap.R -> ibd_overlap_outputs/   (Fisher IBD overlap + comparison)
#
# NOTE: 17_multimir_targets_baseline.R is NOT re-run — the baseline has no S-score
#       gate, so it is unaffected by this change (it already queries every miRNA).
#
# Usage:
#   bash rerun_targets_downstream.sh                 # all datasets with DEA_results/
#   bash rerun_targets_downstream.sh GSE272890       # one or more specific dataset ids
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../pipelines/06_targets_enrich
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"                      # repo root
SCR="$ROOT/pipelines"
OUTROOT="$ROOT/outputs"

# Thresholds (project-wide DE definition). 4th arg to script 14 = s_score_cut;
# omitted on purpose so it uses the new default of 0 (no S-score gate).
PADJ="0.05"
LFC="0.58"
ORG="hsa"

# Which datasets: explicit ids as args, else every outputs/<id>/ with a DEA_results/
if [[ $# -gt 0 ]]; then
  DATASETS=("$@")
else
  DATASETS=()
  for d in "$OUTROOT"/*/; do
    [[ -d "${d}DEA_results" ]] && DATASETS+=("$(basename "$d")")
  done
fi

echo "Repo root : $ROOT"
echo "Datasets  : ${DATASETS[*]}"
echo "Thresholds: padj<=$PADJ  |log2FC|>=$LFC  (s_score gate OFF)"
echo "============================================================"

run() {  # run <label> <Rscript> [args...]
  echo "  -> $1"
  Rscript "$2" "${@:3}"
}

for ds in "${DATASETS[@]}"; do
  dsdir="$OUTROOT/$ds"
  if [[ ! -d "$dsdir/DEA_results" ]]; then
    echo "↷ $ds: no DEA_results/ — skipping"; continue
  fi
  echo "==================== $ds ===================="
  (
    cd "$dsdir" || exit 1
    # Derive targets for the DE-significant set FROM THE BASELINE (offline, no
    # multiMiR query). To re-query multiMiR instead, comment this out and use the
    # line below (run sequentially — do NOT parallelize, the public server throttles):
    run "derive targets from baseline (offline, no multiMiR)" \
        "$SCR/06_targets_enrich/alt1_derive_targets_from_baseline.R" DEA_results "$PADJ" "$LFC"
    # To re-query multiMiR instead (sequential — do NOT parallelize), comment out the
    # line above and use:
    # run "16 multiMiR targets (re-query)" \
    #     "$SCR/06_targets_enrich/16_multimir_targets.R" DEA_results "$PADJ" "$LFC" "" "$ORG"
    run "18 pathway enrichment (GO/KEGG)" \
        "$SCR/07_networks_clin/18_pathway_enrich.R"
    run "19 interactive network" \
        "$SCR/07_networks_clin/19_interactive_network.R"
    run "22 IBD target overlap (+ comparison)" \
        "$SCR/07_networks_clin/22_ibd_target_overlap.R"
  )
  echo "✅ $ds done"
done

echo "============================================================"
echo "🎉 Re-processing complete for: ${DATASETS[*]}"
