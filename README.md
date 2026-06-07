**English** | [Español](README.es.md)

# miRNA-seq meta-analysis pipeline for Inflammatory Bowel Disease (IBD)

A complete end-to-end pipeline for analyzing microRNA sequencing (miRNA-seq) data in
**inflammatory bowel disease** (ulcerative colitis and Crohn's disease) in
*Homo sapiens*. It covers everything from dataset discovery in GEO/SRA to the final
interactive report and a cross-dataset comparison dashboard.

Final-year project (*proyecto terminal*) for the Biotechnology Engineering degree.

---

## Table of contents

1. [Overview](#1-overview)
2. [`datasets.yaml` — master registry](#2-datasetsyaml--master-registry)
3. [Environment setup](#3-environment-setup)
4. [Pipeline stages](#4-pipeline-stages)
   - [00_search — GEO/SRA discovery & curation](#00_search--geosra-discovery--curation)
   - [01_download — SRA run download](#01_download--sra-run-download)
   - [02_qc — Quality control & sample organization](#02_qc--quality-control--sample-organization)
   - [03_quant — miRNA quantification (miRge3)](#03_quant--mirna-quantification-mirge3)
   - [04_merge — Merge, filtering & metadata](#04_merge--merge-filtering--metadata)
   - [05_dea_r — Differential expression (DESeq2)](#05_dea_r--differential-expression-deseq2)
   - [06_targets_enrich — Target prediction (multiMiR)](#06_targets_enrich--target-prediction-multimir)
   - [07_networks_clin — Enrichment, networks & clinical evaluation](#07_networks_clin--enrichment-networks--clinical-evaluation)
   - [08_report — Report, dashboard & diagram](#08_report--report-dashboard--diagram)
5. [Shared library: prep kits](#5-shared-library-prep-kits)
6. [Parallel module: literature search engine](#6-parallel-module-literature-search-engine)
7. [Suggested full workflow](#7-suggested-full-workflow)

---

## 1. Overview

The pipeline is organized into numbered stages under `pipelines/`. Each stage produces
outputs that feed the next, all stored under `outputs/<dataset_id>/`.

| Stage | Folder | Description | Environment |
|---|---|---|---|
| 00 | `00_search` | Dataset discovery and curation in GEO/SRA | `py_env` |
| 01 | `01_download` | SRA run download (prefetch + fasterq-dump) | `py_env` |
| 02 | `02_qc` | QC with FastQC/MultiQC + miRTrace and sample organization | `py_env` |
| 03 | `03_quant` | miRNA quantification with miRge3 | `py_env` |
| 04 | `04_merge` | Count merging, RPM filtering and DESeq2 metadata | `py_env` |
| 05 | `05_dea_r` | Differential expression analysis (DESeq2) | `r_env` |
| 06 | `06_targets_enrich` | miRNA→gene target prediction (multiMiR) | `r_env` |
| 07 | `07_networks_clin` | GO/KEGG enrichment, networks, HMDD, IBD overlap and clinical ROC | `r_env` |
| 08 | `08_report` | HTML report, Streamlit dashboard and pipeline diagram | `r_env` / Python |

There is also a parallel literature-mining module
(`pipelines/search_engine_papers/`) that is independent of the main flow.

**Currently analyzed datasets** (see `datasets.yaml`): `GSE84779`, `PRJNA471862`,
`PRJNA717018`, `GSE114603`, `GSE272890`.

---

## 2. `datasets.yaml` — master registry

`datasets.yaml` (in the repo root) is the project's **single source of truth**: a
registry of every dataset with its accession, BioProject, PMID, tissue, comparisons,
sample counts, output paths and `status` (`analyzed` / `in_progress` / `suspended`).

It ensures nothing has to be inferred from folder-name prefixes. It drives both the
**dashboard** and the cross-dataset **master report**. To add a dataset, copy a block
and fill in the known fields (leave unknown ones as `null`).

---

## 3. Environment setup

### PubMed environment variable (optional)

```bash
export NCBI_API_KEY="your_api_key"
```

Raises the request rate limits when querying PubMed (used by the literature search
engine and by the stage-00 curators).

### Conda environments (WSL / Linux)

The pipeline runs under **WSL** and uses **two environments** (files are in `envs/`):

| Environment | Stages | Contents |
|---|---|---|
| **`py_env`** | 00–04 | miRge3, miRTrace, sra-tools, entrez-direct, FastQC, MultiQC, seqkit, pandas (all Python/CLI/Java tooling). Pinned to **Python 3.8** by miRge3. |
| **`r_env`** | 05–08 | R 4.x + Bioconductor: DESeq2, apeglm, multiMiR, clusterProfiler, enrichplot, org.Hs.eg.db, fgsea, tidyverse, igraph, pROC, etc. |

```bash
conda env create -f envs/py_env.yml   # or: conda env update -n py_env -f envs/py_env.yml
conda env create -f envs/r_env.yml
```

Scripts resolve their own paths (`miRge3_Lib`, `mirtrace.jar`) via `$CONDA_PREFIX`, so
you just need the correct environment active.

### Python dependencies for Windows / dashboard

`requirements.txt` lists **only** the pip packages used outside the conda environments
(e.g. the Streamlit dashboard and the search engine on Windows):

```bash
pip install -r requirements.txt
```

> ⚠️ `pandas<3` is pinned because `upsetplot 0.9.0` (the latest release) breaks under
> pandas 3.0 Copy-on-Write — see `pipelines/08_report/dashboard/upset_compat.py`.

### Tool verification

Before starting, check that the full toolchain is available:

```bash
bash pipelines/00_search/00_tool_check.sh
```

Exits 0 if everything required is present; otherwise exits 1 and writes a log listing
what is missing.

---

## 4. Pipeline stages

### 00_search — GEO/SRA discovery & curation

A set of curators that discover datasets, extract metadata and generate Markdown cards
in `outputs/docs/`, which in turn sync `datasets.yaml`.

| Script | Function |
|---|---|
| `00_tool_check.sh` | Verifies the full toolchain across both conda environments. |
| `01_geo_ibd_curator.py` | Automated discovery: searches, filters and extracts metadata for GEO/SRA datasets via a *true pivot* on the BioProject database. Uses `ignore_list.txt` as a blacklist and generates SRR target lists for download. |
| `02_bioproject_curator.py` | Targeted ("sniper") curation from one or more BioProject IDs (`PRJNA…`), bypassing GEO entirely. |
| `03_geo_id_curator.py` | Targeted curation for explicit GEO accessions (`GSE…`); same output format as `01_` but without the broad search step. |
| `04_dedup_docs.py` | Removes redundant `PRJNA*.md` cards when a `GSE*.md` already exists for the same BioProject (dry-run by default; `--confirm` to delete). |
| `05_sync_yaml.py` | Syncs `datasets.yaml` from `outputs/docs/*.md`: fills `null` fields without overwriting existing values. With `--add-new` it appends new entries. Requires `ruamel.yaml`. |

```bash
# Automated discovery
python pipelines/00_search/01_geo_ibd_curator.py

# Targeted curation by GEO accession or BioProject
python pipelines/00_search/03_geo_id_curator.py --id GSE114603
python pipelines/00_search/02_bioproject_curator.py --id PRJNA904924

# Clean up duplicates and sync the manifest
python pipelines/00_search/04_dedup_docs.py --confirm
python pipelines/00_search/05_sync_yaml.py --add-new --confirm
```

> The `complementary_scripts/` subfolder holds legacy/auxiliary utilities
> (`archive_geo_to_sra.sh`, runinfo post-processing, BioSample attribute fetching, etc.).

---

### 01_download — SRA run download

```bash
bash pipelines/01_download/06_download_runs.sh \
  [runinfo.csv|srr_list.txt] [fastq_dir] [fasterq_threads] [parallel_jobs]
```

- With no arguments it looks for `srr_list.txt` in the current directory and downloads into `sra_fastq/`.
- Uses `prefetch` → `fasterq-dump` and compresses to `.fastq.gz`; writes a local `01_download.log`.

The `complementary_scripts/` subfolder includes variants and utilities: sample-name
correction (`011_correction_name_samples.sh`) and compression (`012_compress_fastq.sh`).

---

### 02_qc — Quality control & sample organization

| Script | Function |
|---|---|
| `07_qc_fastqc.sh` | Runs FastQC over all `*.fastq.gz` and consolidates with MultiQC. |
| `08_create_config_mirtrace.sh` | Generates `mirtrace_config.csv`, picking the 3' adapter from the study's prep kit (see `kits.tsv`). |
| `09_qc_mirtrace.sh` | Runs miRTrace in batches using the kit's protocol/adapter. |
| `10_organize_samples.py` | Classifies FASTQ files into `sra_fastq/muestras_clasificadas/` based on the runinfo. |

```bash
bash pipelines/02_qc/07_qc_fastqc.sh
DATASET=GSE272890 bash pipelines/02_qc/08_create_config_mirtrace.sh
bash pipelines/02_qc/09_qc_mirtrace.sh mirtrace_config.csv
python pipelines/02_qc/10_organize_samples.py
```

---

### 03_quant — miRNA quantification (miRge3)

| Script | Function |
|---|---|
| `11_auto_make_list_subfolders.sh` | Builds the list of sample subfolders to process. |
| `12_auto_subfolders_mirge3.sh` | Safely and sequentially iterates over all subfolders of `muestras_clasificadas/`, running miRge3 one by one. Avoids WSL `/mnt/` paths (much slower). |
| `10_single_folder_mirge3.sh` | Processes a single folder / run list; `THREADS`, `PARALLEL`, `FASTQ_DIR` controls. |

```bash
# From the dataset root
bash pipelines/03_quant/11_auto_make_list_subfolders.sh
bash pipelines/03_quant/12_auto_subfolders_mirge3.sh
```

Per-sample output with `miR.Counts.csv` and `miR.RPM.csv`.

---

### 04_merge — Merge, filtering & metadata

| Script | Function |
|---|---|
| `13_mrina_cut_off.py` | Recursively finds the `miR.Counts.csv` / `miR.RPM.csv` files, unifies them into a matrix and filters by **RPM ≥ 100**. Produces `Master_Filtered_Counts_DESeq2.csv`. |
| `14_deseq2_metadata.py` | Locates the runinfo and the master matrix and generates the aligned `sample_metadata.csv` for DESeq2. |

```bash
python pipelines/04_merge/13_mrina_cut_off.py
python pipelines/04_merge/14_deseq2_metadata.py
```

---

### 05_dea_r — Differential expression (DESeq2)

```bash
Rscript pipelines/05_dea_r/15_run_dea.R
```

- Automatically finds `Master_Filtered_Counts_DESeq2.csv` and `sample_metadata.csv` from
  the dataset root.
- Interactive: shows the available metadata columns and asks which one to analyze.
- Runs DESeq2 with `apeglm` and VST; produces two models (**Default Replace** and
  **Strict No-Replace**), figures (PCA, UMAP, heatmap, volcanos) and DEA tables.

---

### 06_targets_enrich — Target prediction (multiMiR)

```bash
Rscript pipelines/06_targets_enrich/16_multimir_targets.R [dea_csv_or_dir] [padj] [lfc] [s_score_cut] [org]
# s_score_cut: optional S-score floor for the target query. Default 0 = no gate —
# query ALL DE-significant miRNAs (padj/|log2FC| only). Pass e.g. 2.0 to restrict.
```

- Fetches **validated and predicted** miRNA→gene interactions with multiMiR from the
  significant miRNAs in the DEA.
- Ranks them by a custom `S_miRNA` score, splits grouped miRNAs and queries in chunks.
- Writes the tables to `multimir_outputs/`.

---

### 07_networks_clin — Enrichment, networks & clinical evaluation

| Script | Function |
|---|---|
| `18_pathway_enrich.R` | GO (BP/CC/MF) and KEGG enrichment of the target genes with clusterProfiler; *smart file finder* for CLI path inputs. Output in `pathway_outputs/`. |
| `19_interactive_network.R` | Interactive tripartite miRNA → gene → pathway/GO dashboards (Cytoscape.js). v3: filters to strong *Functional MTI* interactions and a hub-focused view to avoid the "hairball". |
| `20_treatment_response_roc.R` | Biomarker-performance evaluation (ROC/AUC) of candidate miRNAs: glucocorticoid responders vs non-responders, and UC vs healthy controls. Output in `clinical_outputs/`. |
| `21_hmdd_prep.R` | Prepares copy-pasteable miRNA lists for **HMDD v4.0** and **miRNet 2.0** and opens both sites in the browser. |
| `22_ibd_target_overlap.R` | Assesses whether validated targets are enriched in known IBD genes via the **OpenTargets** API (GraphQL) + a per-miRNA Fisher's exact test. |

```bash
Rscript pipelines/07_networks_clin/18_pathway_enrich.R
Rscript pipelines/07_networks_clin/19_interactive_network.R
Rscript pipelines/07_networks_clin/20_treatment_response_roc.R
Rscript pipelines/07_networks_clin/21_hmdd_prep.R
Rscript pipelines/07_networks_clin/22_ibd_target_overlap.R
```

---

### 08_report — Report, dashboard & diagram

| Script | Function |
|---|---|
| `23_generate_report.py` | Generates a dataset's interactive `Final_Report.html` (Plotly). Auto-discovers the DEA subfolders (Default/Strict) and all comparisons; integrates QC, differential expression, enrichment, networks and IBD overlap. |
| `24_launch_dashboard.sh` | Launches the Streamlit dashboard (`dashboard/Dashboard.py`). Finds a Python with Streamlit among the project conda environments and serves at `http://localhost:8501`. |
| `pipeline_diagram.py` | Generates the full pipeline diagram for the thesis (`pipeline_diagram.png` at 300 DPI + vector `.pdf`). |

```bash
py -3.11 pipelines/08_report/23_generate_report.py <results_dir>
bash pipelines/08_report/24_launch_dashboard.sh
py -3.11 pipelines/08_report/pipeline_diagram.py
```

#### Comparison dashboard (`dashboard/`)

A **multi-page Streamlit** app driven by `datasets.yaml`: adding a dataset to the
manifest makes it appear automatically.

- `Dashboard.py` — home page with global metrics and the dataset registry table.
- `data.py` — data-access layer; reads the manifest and the on-disk output layout.
- `pages/1_Dataset_Explorer.py` — per-dataset explorer: QC, differential expression, enrichment and IBD overlap.
- `pages/2_Cross_Dataset_Analysis.py` — cross-dataset analysis: overlap (Venn/UpSet), direction matrix, Stouffer meta-analysis, effect-size concordance (Spearman), shared pathways and IBD-target overlap.
- `pages/3_miRNA_Lookup.py` — per-miRNA profile: differential expression across all datasets/comparisons plus IBD target-overlap relevance and the IBD genes it targets.
- `upset_compat.py` — compatibility shims for `upsetplot 0.9.0` under modern pandas/matplotlib.

---

## 5. Shared library: prep kits

The **library-prep kit** is a property of the study (every run uses the same prep), so
it is assigned per dataset in `pipelines/kits.tsv` rather than auto-detected from reads.

- `pipelines/kits.tsv` — `dataset_id → kit` map. Known kits: `illumina`,
  `nextflex_v3` (Illumina + 4,4 UMIs), `nebnext`. Anything unlisted defaults to `illumina`.
- `pipelines/lib/kit_params.sh` — library that is `source`d to resolve the kit
  (`resolve_kit`) and load its parameters (`load_kit_params`): adapters and protocol for
  miRge3 and miRTrace. Used by the stage-02 scripts.

> Note: NEXTflex v3 shares the exact 3' adapter with standard Illumina and differs only
> by the 4+4 UMIs, which adapter-sniffing cannot detect — hence the explicit per-dataset
> mapping.

---

## 6. Parallel module: literature search engine

`pipelines/search_engine_papers/` is a standalone module for mining scientific
literature (PubMed):

- `PubMed_API_0.1.py` — queries the PubMed API (uses `NCBI_API_KEY`).
- `interactive_cli.py` — interactive assistant to collect search parameters and run queries.
- `pt_search.py`, `text_utils.py`, `paper_processing/refineria.py` — TF-IDF search and refining/processing of the retrieved documents.

See its own `README.md` inside the folder for details.

---

## 7. Suggested full workflow

1. `00_search` → discover/curate datasets and sync `datasets.yaml`.
2. `01_download` → download `.fastq.gz`.
3. `02_qc` → QC (FastQC/MultiQC + miRTrace) and organize samples.
4. `03_quant` → per-miRNA counts with miRge3.
5. `04_merge` → filtered count matrix + DESeq2 metadata.
6. `05_dea_r` → differentially expressed miRNAs.
7. `06_targets_enrich` → gene targets with multiMiR.
8. `07_networks_clin` → enrichment, networks, HMDD, IBD overlap and clinical ROC.
9. `08_report` → per-dataset `Final_Report.html` + cross-dataset comparison dashboard.

Each script includes further comments and logs you can consult for additional
customization. Good luck with your analysis!
