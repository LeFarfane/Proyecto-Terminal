#!/usr/bin/env Rscript
# =============================================================================
# Script : 07_networks_clin/hmdd_prep.R
# Purpose: Prepare copy-pasteable miRNA lists for HMDD v4.0 AND miRNet 2.0,
#          then open both sites in the browser automatically.
#
#  HMDD v4.0  — https://www.cuilab.cn/hmdd
#  ─────────────────────────────────────────
#  1. Batch Search  — semicolon-delimited list (max 20 keywords).
#  2. Disease Network — two lists: upregulated / downregulated miRNAs.
#
#  miRNet 2.0  — https://www.mirnet.ca/miRNet/upload/MirUploadView.xhtml
#  ────────────────────────────────────────────────────────────────────────
#  Accepts a tab-separated file (miRNA ID <tab> log2FC) or plain ID list.
#  Use case: miRNA → target → pathway / disease network analysis.
#
#  Strategy
#  ────────
#  Candidate miRNAs = those significant in either DEA comparison
#                     (padj ≤ 0.05, |log2FC| ≥ 1).
#  Direction is determined by the most significant comparison.
#  Compound names (e.g. hsa-miR-378a-3p/378c/378d/378e) are simplified
#  to their primary isoform for external tool compatibility.
#
#  Outputs  (DEA_results_round_2/clinical_outputs/)
#  ─────────────────────────────────────────────────
#  hmdd_lists.txt          — copy-paste sections for HMDD and miRNet
#  hmdd_mirna_table.csv    — machine-readable candidate table
#  mirnet_input.tsv        — miRNA ID + log2FC, ready to upload to miRNet
# =============================================================================

suppressWarnings(suppressMessages({ library(utils); library(tools) }))

ts      <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse="")))

# ── paths ─────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(getwd())
}
SCRIPT_DIR <- tryCatch(get_script_dir(), error = function(e) normalizePath(getwd()))

# Run from the dataset root directory — finds DEA CSVs and metadata recursively
RUN_DIR  <- normalizePath(getwd())
OUT_DIR  <- file.path(RUN_DIR, "hmdd_outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dea_files <- list.files(RUN_DIR, pattern = "^DEA_.*\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(dea_files) == 0) {
  log_msg("ERROR: No DEA_*.csv files found under ", RUN_DIR, ". Run from the dataset root directory.")
  quit(status = 1)
}
log_msg(sprintf("Found %d DEA files:", length(dea_files)))
for (f in dea_files) log_msg("  ", basename(f))

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.58
MAX_BATCH   <- 20     # HMDD batch search hard limit

URL_HMDD   <- "https://www.cuilab.cn/hmdd"
URL_MIRNET <- "https://www.mirnet.ca/miRNet/upload/MirUploadView.xhtml"

# ── browser helper ────────────────────────────────────────────────────────────
# Opens a URL in the default browser.
# Handles: regular Linux/Mac, WSL (Windows Subsystem for Linux), Windows.
open_url <- function(url, label = url) {
  log_msg(sprintf("Opening: %s", label))

  # Detect WSL by checking /proc/version for "microsoft" or "WSL"
  is_wsl <- tryCatch({
    pv <- readLines("/proc/version", n = 1L, warn = FALSE)
    grepl("microsoft|WSL", pv, ignore.case = TRUE)
  }, error = function(e) FALSE)

  opened <- FALSE

  if (is_wsl) {
    # In WSL, use powershell as the single reliable method to avoid double-opens
    cmd <- sprintf("powershell.exe -Command \"Start-Process '%s'\"", url)
    ret <- tryCatch(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE),
                    error = function(e) 1L)
    opened <- (ret == 0L)
  } else {
    # Standard R browser opener (works on Linux with xdg-open, macOS, Windows)
    tryCatch({ browseURL(url); opened <- TRUE }, error = function(e) NULL)
  }

  if (!opened) {
    message(sprintf(
      "\n  Could not open browser automatically.\n  Please navigate to:\n  %s\n", url))
  }
  invisible(opened)
}

log_msg("hmdd_prep.R starting")

suppressWarnings(suppressMessages({ library(dplyr); library(readr) }))

# ── load & combine all DEA files ──────────────────────────────────────────────
dea_list <- lapply(dea_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$contrast <- tools::file_path_sans_ext(basename(f))
  df
})
dea_all <- bind_rows(dea_list)

# ── filter significant miRNAs across all contrasts ────────────────────────────
sig_all <- dea_all %>%
  filter(!is.na(padj), padj <= PADJ_THRESH, abs(log2FoldChange) >= LFC_THRESH)

if (nrow(sig_all) == 0) {
  log_msg("ERROR: No significant miRNAs found (padj <= ", PADJ_THRESH,
          ", |log2FC| >= ", LFC_THRESH, "). Exiting.")
  quit(status = 1)
}

# Per miRNA: keep the entry with the lowest padj across all contrasts
all_cand <- sig_all %>%
  group_by(miRNA) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    best_padj = padj,
    best_lfc  = log2FoldChange,
    direction = if_else(log2FoldChange > 0, "UP", "DOWN")
  ) %>%
  arrange(best_padj)

log_msg(sprintf("Total candidates: %d  (%d UP, %d DOWN)",
                nrow(all_cand),
                sum(all_cand$direction == "UP"),
                sum(all_cand$direction == "DOWN")))

# ── simplify compound miRNA names for HMDD ────────────────────────────────────
# HMDD accepts individual isoform names.
# Compound names like "hsa-miR-378a-3p/378c/378d/378e" → "hsa-miR-378a-3p"
# "hsa-let-7a-5p/7c-5p"                                 → "hsa-let-7a-5p"
simplify_name <- function(x) sub("/.*", "", x)

all_cand <- all_cand %>%
  mutate(miRNA_hmdd = simplify_name(miRNA))

# ── build HMDD lists ──────────────────────────────────────────────────────────
up_mirnas   <- all_cand %>% filter(direction == "UP")   %>% pull(miRNA_hmdd)
down_mirnas <- all_cand %>% filter(direction == "DOWN")  %>% pull(miRNA_hmdd)
all_mirnas  <- all_cand %>% pull(miRNA_hmdd)

# Batch search: top MAX_BATCH by best_padj (across both directions)
batch_list  <- head(all_mirnas, MAX_BATCH)
remainder   <- tail(all_mirnas, max(0, length(all_mirnas) - MAX_BATCH))

# ── write human-readable text file ───────────────────────────────────────────
out_txt <- file.path(OUT_DIR, "hmdd_lists.txt")
{
  sink(out_txt)

  cat("================================================================================\n")
  cat("  HMDD v4.0 miRNA Lists\n")
  cat("  Generated by hmdd_prep.R\n")
  cat(sprintf("  Date: %s\n", Sys.Date()))
  cat(sprintf("  Criteria: padj <= %.2f  |log2FC| >= %.1f  (best contrast per miRNA)\n",
              PADJ_THRESH, LFC_THRESH))
  cat(sprintf("  Total candidates: %d  (%d UP, %d DOWN)\n",
              nrow(all_cand), length(up_mirnas), length(down_mirnas)))
  cat("================================================================================\n\n")

  # ── Section 1: Batch Search ──────────────────────────────────────────────
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("SECTION 1 — BATCH SEARCH\n")
  cat("  Site : hmdd.cuilab.cn  →  Search  →  Batch Search\n")
  cat("  Limit: max 20 semicolon-separated keywords\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  cat(sprintf("List A — top %d (by padj) — COPY PASTE BELOW:\n\n", MAX_BATCH))
  cat(paste(batch_list, collapse = ";"))
  cat("\n\n")

  if (length(remainder) > 0) {
    cat(sprintf("List B — remaining %d miRNAs (second query if needed):\n\n",
                length(remainder)))
    cat(paste(remainder, collapse = ";"))
    cat("\n\n")
  }

  # ── Section 2: Disease Network ───────────────────────────────────────────
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("SECTION 2 — DISEASE NETWORK\n")
  cat("  Site : hmdd.cuilab.cn  →  Analysis  →  miRNA Disease Network\n")
  cat("  Input: two lists (upregulated & downregulated)\n")
  cat("  Note : direction = UP/DOWN relative to UC vs Healthy Controls\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  cat(sprintf("UPREGULATED in UC  (n=%d) — paste into 'Up-regulated miRNA list':\n\n",
              length(up_mirnas)))
  cat(paste(up_mirnas, collapse = "\n"))
  cat("\n\n")

  cat(sprintf("DOWNREGULATED in UC  (n=%d) — paste into 'Down-regulated miRNA list':\n\n",
              length(down_mirnas)))
  cat(paste(down_mirnas, collapse = "\n"))
  cat("\n\n")

  # ── Section 3: miRNet ────────────────────────────────────────────────────
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("SECTION 3 — miRNet 2.0\n")
  cat("  Site : https://www.mirnet.ca/miRNet/upload/MirUploadView.xhtml\n")
  cat("  Steps: (1) Select 'miRNA'  (2) ID type: 'miRBase'  (3) Organism: 'Homo sapiens'\n")
  cat("         (4) Upload the file: mirnet_input.tsv   OR   paste the list below\n")
  cat("  Note : miRNet will map targets → build interaction network → enrichment\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  cat("miRNet paste format (miRNA<tab>log2FC) — ALL candidates:\n\n")
  for (i in seq_len(nrow(all_cand))) {
    cat(sprintf("%s\t%.4f\n", all_cand$miRNA_hmdd[i], all_cand$best_lfc[i]))
  }
  cat("\n")

  # ── Section 4: Full reference table ─────────────────────────────────────
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("SECTION 4 — FULL CANDIDATE TABLE (reference)\n")
  cat("  Format: miRNA | direction | best_padj | best_log2FC | contrast\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  fmt <- "  %-35s  %-5s  %10s  %9s  %-30s\n"
  cat(sprintf(fmt, "miRNA", "dir", "best_padj", "log2FC", "contrast"))
  cat(sprintf(fmt, paste(rep("-", 35), collapse=""),
              "-----", "----------", "---------", paste(rep("-", 30), collapse="")))
  for (i in seq_len(nrow(all_cand))) {
    r <- all_cand[i, ]
    cat(sprintf(fmt,
                r$miRNA,
                r$direction,
                formatC(r$best_padj, format = "e", digits = 2),
                formatC(r$best_lfc,  format = "f", digits = 3),
                r$contrast))
  }
  cat("\n")

  sink()
}

log_msg("Saved: hmdd_lists.txt")

# ── write machine-readable CSV ────────────────────────────────────────────────
readr::write_csv(
  all_cand %>% select(miRNA, miRNA_hmdd, direction, best_padj, best_lfc, contrast),
  file.path(OUT_DIR, "hmdd_mirna_table.csv")
)
log_msg("Saved: hmdd_mirna_table.csv")

# ── write miRNet upload file (tab-separated: miRNA_ID <tab> log2FC) ───────────
mirnet_df <- all_cand %>%
  select(miRNA_ID = miRNA_hmdd, log2FC = best_lfc) %>%
  arrange(desc(abs(log2FC)))

readr::write_tsv(mirnet_df, file.path(OUT_DIR, "mirnet_input.tsv"),
                 col_names = FALSE)   # miRNet expects no header row
log_msg(sprintf("Saved: mirnet_input.tsv  (%d miRNAs)", nrow(mirnet_df)))

# ── print preview to console ──────────────────────────────────────────────────
log_msg("\n===== HMDD BATCH SEARCH (copy this) =====")
cat(paste(batch_list, collapse = ";"), "\n\n")
log_msg(sprintf("UP   (%d): %s", length(up_mirnas),   paste(up_mirnas,   collapse="; ")))
log_msg(sprintf("DOWN (%d): %s", length(down_mirnas), paste(down_mirnas, collapse="; ")))

# ── open browser windows ──────────────────────────────────────────────────────
log_msg("\n===== Opening browser tabs =====")
Sys.sleep(0.5)
open_url(URL_HMDD,   "HMDD v4.0")
Sys.sleep(1.5)   # slight delay so both tabs open cleanly
open_url(URL_MIRNET, "miRNet 2.0")

log_msg("hmdd_prep.R completed.")
log_msg(sprintf("  → HMDD   : %s", URL_HMDD))
log_msg(sprintf("  → miRNet : %s", URL_MIRNET))
