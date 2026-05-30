# Repository Analysis — miRNA Analysis Pipeline (corrected, full)

*Analysis date: 2026-05-28 · Scope: architecture & workflow, scientific methodology, code quality, documentation & organization.*

> **Note:** This version supersedes an earlier draft that was based on an incomplete read of the repository (a tooling issue truncated the file listing). That draft wrongly claimed quantification, differential expression, target prediction, enrichment, compression, a top-level README, and a `.gitignore` were missing. **They all exist.** This document reflects the full tree: a complete, 16-step end-to-end pipeline. The working tree is also well ahead of the last commit (`22e19d4`, 2025-09-25), with many uncommitted changes — consistent with "haven't committed in a long time."

---

## 1. Executive summary

This is a complete, end-to-end human microRNA-seq pipeline for inflammatory bowel disease (IBD: Crohn's disease, ulcerative colitis, treatment responders/non-responders, healthy controls), plus a companion literature-mining module. It runs from automated dataset discovery all the way to an interactive biomarker network dashboard:

**dataset discovery (GEO/SRA) → download → QC → quantification (miRge3) → count merge/filter → differential expression (DESeq2) → target prediction (multiMiR) → pathway enrichment (clusterProfiler) → network construction + clinical/biomarker evaluation + interactive dashboard.**

The technical quality is high for a thesis project and in places genuinely impressive: the SRA downloader has prefetch/fallback, parallelism, resume logic and compression; the DESeq2 stage does VST, PCA/UMAP QC, dual outlier-handling runs, apeglm shrinkage and a custom ranking score; the multiMiR stage handles validated/predicted/disease-association queries with miRNA-family expansion; and stage 16 emits a self-contained vis.js tripartite (miRNA→gene→pathway) dashboard. The literature module includes a fully featured local TF-IDF search engine.

The real weaknesses are now about **integration, reproducibility, and consistency**, not missing science:

- **The README has drifted from the code.** It documents idealized script names and a dedicated `cutadapt` trimming step (`trim_qc.sh`) that don't exist on disk; the actual scripts are numbered `00`–`16` with different names.
- **Inter-stage filenames/paths don't line up by default.** The DEA stage writes to `results/DEA_No_Outlier_Replacement/` etc., while `clinical_eval.R` defaults to `results/DEA/normalized_counts_deseq2.csv`; `multimir_targets.R` writes to `multimir_outputs/`, while `build_network.R` defaults to `results/targets_validated.csv`. The pipeline runs only if you pass the right arguments by hand.
- **Two execution conventions coexist:** most scripts assume "`cd` into the dataset folder and run here" (`getwd()`), while the R defaults assume a centralized `results/`/`data/` layout. There's no orchestrator to reconcile them.
- **Portability gaps:** a hard-coded miRge3 library path (`/home/genesis/miniconda3/...`), several interactive `input()`/`readLines` prompts that block automation, and a nested git repo inside `pipelines/00_search/`.
- A few concrete bugs remain (corrupted miRTrace command line; `refineria.py` substring matching and overwrite-instead-of-append).

In short: the science is essentially all here and well done. What's missing is the connective tissue (consistent paths, an orchestration layer, a refreshed README, environment pinning) that turns a set of strong stage scripts into a one-command reproducible pipeline.

---

## 2. Actual repository structure

```
Pipelines/
├── README.md                       # comprehensive pipeline instructivo (but drifted from code)
├── .gitignore                      # ignores *.log, __pycache__/, .venv/, outputs/*
├── requirements.txt                # numpy, pandas, scikit-learn, matplotlib, networkx, requests
├── pipelines/
│   ├── 00_search/                  # automated GEO/SRA dataset discovery
│   │   ├── 00_tool_check.sh                 # verify EDirect + SRA Toolkit on PATH
│   │   ├── 01_geo_ibd_curator.py            # GEO→SRA curator (argparse, blacklist, provenance)
│   │   ├── 02_bioproject_curator.py         # BioProject-centric curator
│   │   └── complementary_scripts/           # Entrez runinfo, biosample attrs, postproc, debug
│   │   └── .git/                            # ⚠ nested git repository
│   ├── 01_download/
│   │   ├── 03_download_runs.sh              # prefetch→fasterq-dump, parallel, pigz, resume
│   │   └── complementary_scripts/010,011,012  # older download / rename / compress
│   ├── 02_qc/
│   │   ├── 04_qc_fastqc.sh                  # FastQC + MultiQC
│   │   ├── 05_create_config_mirtrace.sh    # build miRTrace config CSV
│   │   └── 06_qc_mirtrace.sh               # run miRTrace   ⚠ corrupted command line
│   ├── 03_quant/
│   │   ├── 07_auto_make_list_subfolders.sh # list fastq per subfolder
│   │   ├── 08_test_auto_mirge3_subfolders.sh  # miRge3 over nested folders (trim+quant)
│   │   └── 09_single_end_runs_mirge3.sh    # miRge3 per-run from a list
│   ├── 04_merge/
│   │   ├── 11_mrina_cut_off.py             # merge miR.RPM/Counts, filter avg RPM≥100
│   │   └── 12_deseq2_metadata.py          # build aligned sample_metadata.csv
│   ├── 05_dea_r/
│   │   └── 13_run_dea.R                    # DESeq2: VST, PCA/UMAP, dual runs, apeglm, S_miRNA
│   ├── 06_targets_enrich/
│   │   └── 14_multimir_targets.R          # multiMiR validated/predicted/disease
│   ├── 07_networks_clin/
│   │   ├── 15_pathway_enrich.R            # GO + KEGG (clusterProfiler)  ⚠ documented under 06
│   │   ├── 16_interactive_network.R      # vis.js tripartite HTML dashboard
│   │   ├── build_network.R               # igraph miRNA–gene network, centrality, modules
│   │   └── clinical_eval.R               # Spearman/FDR, ROC/AUC, 5-fold CV logistic
│   └── papersrch/                        # PubMed mining + local search engine
│       ├── PubMed_API_0.1.py · pt_search.py · text_utils.py · interactive_cli.py
│       └── paper_processing/refineria.py # keyword include/exclude filter  ⚠ bugs
├── data_raw/ · outputs/ · old_data/ · "Scripts for E-direct/" · run_test_1.csv
```

Note the numbering quirks: scripts run `00`→`16` but **step 10 is skipped**, `08` is labeled "test," and the script numbers don't map cleanly onto the `00`–`07` stage folders (e.g. enrichment is script `15` but lives in the `07` folder while the README documents it under `06`).

---

## 3. Architecture & workflow

### What works well

This is a real, scientifically complete miRNA-seq workflow with thoughtful touches at almost every stage. Highlights:

- **Stage 00** is unusually sophisticated for a student pipeline: `01_geo_ibd_curator.py` (≈470 lines) automates GEO→SRA discovery with a curated MeSH query, a persistent blacklist (`ignore_list.txt`), rich per-dataset metadata extraction, Markdown dataset docs, and provenance/run logging. `00_tool_check.sh` pre-flights the EDirect/SRA toolchain.
- **Stage 01** (`03_download_runs.sh`) is robust: `set -euo pipefail`, command checks, `prefetch` with a direct-`fasterq-dump` fallback, parallel downloads via `xargs -P`, skip/resume if files already exist, and `pigz`/`gzip` compression — exactly the operational details that matter at scale.
- **Stage 03** uses miRge3.0, which performs adapter trimming, collapsing, and miRBase quantification in one step (so the absence of a standalone `cutadapt` script is fine in practice — trimming happens inside miRge3 with `-a illumina -q 20 -m 18`).
- **Stage 05** (`13_run_dea.R`) is the analytical centerpiece and is well done (see §4).
- **Stage 07** closes the loop with an igraph network, a clinical/biomarker evaluation, and a polished standalone interactive dashboard.

### The main structural problems

**1. The README no longer matches the code.** The instructivo documents clean script names — `02_qc/trim_qc.sh`, `03_quant/mirge3_run.sh`/`mirdeep2_run.sh`, `04_merge/merge_counts.py`, `05_dea_r/run_dea.R`, `06_targets_enrich/multimir_targets.R`/`enrichment.R` — but the disk has `04_qc_fastqc.sh`, `08/09_*mirge3*`, `11_mrina_cut_off.py`/`12_deseq2_metadata.py`, `13_run_dea.R`, `14_multimir_targets.R`, `15_pathway_enrich.R`. It also describes a dedicated cutadapt trimming script and a miRDeep2 path that aren't present. A reader following the README literally will not find the files.

**2. Inter-stage I/O contracts don't align by default.** The stages were clearly developed somewhat independently, and their default input/output paths don't match:

| Producer | Writes (default) | Consumer | Expects (default) |
|---|---|---|---|
| `13_run_dea.R` | `results/DEA_No_Outlier_Replacement/DEA_Strict_*.csv`, `…/vst_normalized_matrix.csv` | `clinical_eval.R` | `results/DEA/normalized_counts_deseq2.csv`, `DEA_DESeq2_*_vs_*.csv` |
| `14_multimir_targets.R` | `multimir_outputs/targets_validated.csv` | `build_network.R` | `results/targets_validated.csv` |
| `11_mrina_cut_off.py` | `Master_Filtered_Counts_DESeq2.csv` (cwd) | `13_run_dea.R` | `Master_Filtered_Counts_DESeq2.csv` (cwd) ✓ |

Some links do work (e.g. `multimir_targets.R` reads `DEA_*.csv` and the `S_miRNA` column that `13_run_dea.R` produces; `16_interactive_network.R` uses recursive file-finding and is robust). But others require the user to pass explicit arguments, and the normalized-count filename mismatch between `13` and `clinical_eval.R` means the clinical step won't find its input out of the box.

**3. No orchestration, two conventions.** Every stage is a hand-run script. Most expect you to `cd` into a dataset directory and operate on `getwd()` (creating `mirge3_output/`, `multimir_outputs/`, `pathway_outputs/`, `network_outputs/` locally), while the older R scripts (`build_network.R`, `clinical_eval.R`) compute a `ROOT_DIR` and read/write a centralized `results/`+`data/` layout. There's no Snakemake/Nextflow/Make to declare dependencies, reconcile these conventions, or re-run only what changed.

**4. Repo hygiene.** There is a **nested git repository** at `pipelines/00_search/.git` (not registered as a submodule), which will make the parent repo treat that subtree inconsistently. `outputs/` is correctly gitignored now, but historical outputs may already be in history. `"Scripts for E-direct/"` (spaces in the name) and `run_test_1.csv` at the root look like scratch material.

---

## 4. Scientific methodology

### Quantification — appropriate

miRge3.0 against miRBase is a standard, well-regarded choice for small-RNA quantification, and running it per-sample with built-in trimming/collapsing is reasonable. The parameters (`-q 20`, `-m 18`, `-a illumina`, `-ex 0.1`) are sensible defaults for human small-RNA libraries. The only portability issue is the hard-coded `-lib /home/genesis/miniconda3/envs/mirge3/miRge3_Lib` (see §5).

### Count filtering — reasonable but aggressive, and global

`11_mrina_cut_off.py` keeps miRNAs whose **mean RPM across all samples ≥ 100**. That is a fairly strict abundance gate, and because it's a *global* average it can drop miRNAs that are highly expressed in one group but near-zero elsewhere — which is exactly the kind of group-specific signal a case/control study may care about. Consider a less aggressive threshold and/or a "expressed in ≥ N samples of any group" rule (e.g. CPM/RPM ≥ k in ≥ smallest-group-size samples), which is the more standard edgeR/DESeq2 pre-filtering idiom.

### Differential expression — strong

`13_run_dea.R` is the best-engineered analysis in the repo: variance-stabilizing transformation, PCA and UMAP QC plots, a heatmap, and then **two DESeq2 runs** — a "strict" model (`minReplicatesForReplace = Inf`) and the default Cook's-distance outlier-replacement model — which is a thoughtful way to show robustness. It uses apeglm LFC shrinkage for effect sizes while retaining the original Wald `stat`, builds a design that automatically includes available covariates (batch/age/sex/tissue/platform), and computes a custom `S_miRNA = log2FoldChange × −log10(padj)` ranking score. Good statistical hygiene overall. Minor notes: `fitType="local"` is forced (fine, but worth confirming the dispersion fit visually — the script does save dispersion plots, good); and the interactive `readLines` prompt for the condition column blocks non-interactive runs.

### Target prediction & enrichment — solid

`14_multimir_targets.R` queries multiMiR for validated and predicted interactions plus disease/drug associations (IBD/UC/CD), chunks queries politely, and includes a neat "slash fix" that expands grouped miRNA family labels (e.g. `hsa-miR-1/2/3`) into individual IDs before querying. `15_pathway_enrich.R` does GO (BP/CC/MF) and KEGG enrichment via clusterProfiler with BH correction and symbol→Entrez conversion. Both are standard, defensible choices.

### Network & clinical evaluation — good, with two cautions

`build_network.R` (igraph bipartite network + centrality + Louvain modules + Cytoscape SIF) and `16_interactive_network.R` (tripartite vis.js dashboard) are good exploratory tools. `clinical_eval.R` runs Spearman correlations with BH FDR, per-miRNA ROC/AUC, and a 5-fold CV multi-miRNA logistic model. Two methodological cautions carry over:

- **Feature-selection leakage** in the multi-miRNA model: the "top 5 by AUC" are chosen using AUCs computed on all samples and then evaluated by CV on the same samples, which inflates `AUC_CV`. Move feature selection inside each fold (nested CV) or use a held-out test set.
- **Per-miRNA single-feature AUCs are in-sample/optimistic** and `direction="<"` is fixed (assumes higher expression = case); consider `direction="auto"` and clear caveats.

Also note Louvain modularity on a bipartite graph is an approximation — fine for hypothesis generation, less so if modules become a headline result.

### Literature module — capable, with filtering bugs

The PubMed fetcher and the local TF-IDF `SearchEngine` (separate title/abstract indices, boolean candidate selection, recency/review/domain/proximity boosts, faceting, co-author network) are well beyond a typical helper script. `refineria.py`, however, has real bugs that will distort the curated corpus (see §5).

---

## 5. Code quality

### Strengths

Defensive shell scripting is the norm in the ingestion stages (`set -euo pipefail`, dependency checks, failure logs, resume logic, rate-limit pauses). The Python is generally clean and the heavier R scripts are modular, logged, and parameterized with sensible defaults. The DEA and multiMiR scripts in particular show care (shrinkage, dual runs, chunked API calls, column-name normalization).

### Concrete bugs and issues

- **`02_qc/06_qc_mirtrace.sh` — corrupted command (high priority).** Line 3 mashes two invocations together:
  ```
  mirtrace qc --species hsa -c "$1" --output-dir mirtrace_resultsmirtrace --species hsa --config mirtrace_config.csv --output-dir mirtrace_results
  ```
  It should be a single clean call, e.g. `mirtrace qc --species hsa --config "$1" --output-dir mirtrace_results`.

- **Hard-coded miRge3 library path (portability).** `08`/`09` pass `-lib /home/genesis/miniconda3/envs/mirge3/miRge3_Lib`. This breaks on any machine that isn't "genesis." Make it a variable/argument with an env-var default (e.g. `${MIRGE3_LIB:?}`).

- **`09_single_end_runs_mirge3.sh` shebang is on line 2.** Line 1 is blank, so the `#!/usr/bin/env bash` shebang is ignored unless the file is invoked explicitly with `bash`. Move it to line 1.

- **`refineria.py` filtering bugs.**
  - Substring matching (`any(k in title …)`) over-filters: `'ct'` matches *function, structure, detection*; `'oral'` matches *temporal*; `'oil'` matches *boiling*; `'drug'`/`'diet'`/`'management'` are broad. Use whole-word/token matching.
  - `save_to_dumpster` uses `mode="w"` in both branches, so it **overwrites instead of appends** (the code even admits this: *"THIS SHOULD BE APPEND BUT … TOO LAZY TO FIX IT"*). Each run discards prior results.
  - Hard-coded absolute paths (`/data/papers/...`) and a nested-quote f-string (`f"{row["Title"]}"`) that only parses on Python 3.12+ make it fragile/non-portable.

- **`02_auto_get_biosample_attributes.sh` can silently drop columns.** The CSV header is built from the **first** BioSample's attributes only; attributes that later samples have but the first lacks never appear. Build the header from the union across samples.

- **Interactive prompts block automation.** `03_runinfo_postp.py`, `12_deseq2_metadata.py`, `13_run_dea.R`, and the `00_search` curators all use `input()`/`readLines`. Fine for manual use, but they prevent unattended/orchestrated runs — convert to argparse/`commandArgs` with the prompt as a fallback.

- **Stale header comments.** `16_interactive_network.R` is titled "Script: 17_interactive_network.R"; several headers reference old names. Harmless but adds to the drift.

### Style

Comments mix Spanish and English. Configuration constants (species `hsa`, adapter, RPM threshold `100`, S-score cutoffs, year cutoffs, score weights) are scattered across files rather than centralized. The `getwd()`-relative vs `ROOT_DIR`-relative split (noted in §3) is the biggest maintainability snag.

---

## 6. Documentation & organization

There **is** a substantial top-level `README.md` with a full step-by-step instructivo, environment notes, and per-stage examples — a real strength. The problem is drift: it documents an idealized naming scheme and a couple of scripts/paths that don't exist on disk, so it can't currently be followed verbatim. Bringing the README in line with the actual `00`–`16` scripts (or renaming the scripts to match the README) would by itself remove much of the confusion.

`.gitignore` exists and sensibly excludes logs, `__pycache__/`, `.venv/`, and `outputs/*`. Gaps that remain: `requirements.txt` is unpinned and covers only the `papersrch` Python dependencies — the CLI toolchain (EDirect, SRA Toolkit, FastQC, MultiQC, miRTrace, miRge3, pigz) and the R stack (DESeq2, apeglm, umap, pheatmap, tidyverse, multiMiR, clusterProfiler, org.Hs.eg.db, enrichplot, igraph, pROC, jsonlite) are documented only in prose in the README. There is no conda `environment.yml` or `renv.lock`, which matters a lot for bioinformatics reproducibility. The nested `.git` in `00_search`, the space-named `"Scripts for E-direct/"` folder, and root-level `run_test_1.csv` are organizational loose ends.

---

## 7. Prioritized recommendations

**High priority — correctness & integration**

1. Fix the corrupted `06_qc_mirtrace.sh` command line.
2. Replace the hard-coded miRge3 `-lib` path with a variable/argument (env default).
3. Reconcile inter-stage filenames so the pipeline runs without manual path-passing — most importantly, align the DEA normalized-count output with what `clinical_eval.R` reads, and the multiMiR output location with what `build_network.R` reads.
4. Fix `refineria.py` (whole-word matching; real append mode; relative paths).

**Medium priority — reproducibility**

5. Update the README to match the actual `00`–`16` scripts (or rename scripts to the README's scheme) and document the "run from inside the dataset folder" convention explicitly.
6. Pin environments: a conda `environment.yml` for the CLI/Python/R bioinformatics stack and `renv` for R; version `requirements.txt`.
7. Convert interactive prompts to arguments (keeping prompts as fallback) so stages can be orchestrated.
8. Remove/convert the nested `pipelines/00_search/.git` (make it a normal subtree or a proper submodule).

**Lower priority — robustness & polish**

9. Reconsider the global mean-RPM ≥ 100 filter in favor of a group-aware "expressed in ≥ N samples" rule.
10. Address feature-selection leakage in `clinical_eval.R` (nested CV) and the fixed ROC direction.
11. Add a thin orchestration layer (Snakemake/Nextflow) once the I/O contracts are aligned; fix stale header comments; tidy root-level scratch files.
12. **Commit your work** — the tree is far ahead of the last commit; committing now preserves all of the above before further changes.

---

*Bottom line: this is a complete and, in several places, sophisticated miRNA-seq + literature-mining pipeline. The science is essentially all present and largely sound. The remaining work is engineering polish — aligning the stages' inputs/outputs, refreshing the README, pinning the environment, and committing — to make the whole thing reproducible with minimal manual wiring.*
