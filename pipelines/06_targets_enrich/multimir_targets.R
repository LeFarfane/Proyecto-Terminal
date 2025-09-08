#!/usr/bin/env Rscript
# ===========================================
# Script: 06_targets_enrich/multimir_targets.R
# Descripción: Fetches validated and predicted miRNA→target interactions for human
#              using multiMiR, based on significant miRNAs from DEA results.
#              Reads all DEA_DESeq2_* CSVs (or a specific file), filters by padj/log2FC,
#              queries multiMiR in chunks, and writes validated/predicted target tables.
# Fecha: 2025-09-07
# ===========================================

suppressWarnings(suppressMessages({
  library(utils)
  library(tools)
}))

# --- Simple logging ---
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- function(...) cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = "")))

# --- Paths (relative to PT/) ---
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile))
  return(normalizePath(getwd()))
}
SCRIPT_PATH <- tryCatch(get_script_path(), error = function(e) normalizePath(getwd()))
ROOT_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), "..", ".."))  # .../PT
DEA_DIR  <- file.path(ROOT_DIR, "results", "DEA")
OUT_DIR  <- file.path(ROOT_DIR, "results")
LOG_DIR  <- file.path(ROOT_DIR, "outputs", "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOG_DIR, "06_multimir.log")
sink(LOGFILE, append = FALSE, split = TRUE)

# --- Args ---
args <- commandArgs(trailingOnly = TRUE)
# usage:
#   Rscript multimir_targets.R [dea_csv_or_dir] [padj_cut] [lfc_cut] [top_n] [org]
dea_input <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else DEA_DIR
padj_cut  <- if (length(args) >= 2 && nzchar(args[[2]])) as.numeric(args[[2]]) else 0.05
lfc_cut   <- if (length(args) >= 3 && nzchar(args[[3]])) as.numeric(args[[3]]) else 1.0
top_n     <- if (length(args) >= 4 && nzchar(args[[4]])) as.integer(args[[4]]) else 30
org_code  <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else "hsa"  # human

log_line("Starting multimir_targets.R")
log_line(sprintf("dea_input=%s | padj_cut=%.3f | lfc_cut=%.2f | top_n=%d | org=%s",
                 dea_input, padj_cut, lfc_cut, top_n, org_code))

# --- Libraries needed ---
suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(multiMiR)
}))

# --- Helpers ---
list_dea_files <- function(path_or_file) {
  if (file.exists(path_or_file) && file.info(path_or_file)$isdir) {
    fs <- list.files(path_or_file, pattern = "^DEA_DESeq2_.*_vs_.*\\.csv$", full.names = TRUE)
  } else if (file.exists(path_or_file)) {
    fs <- path_or_file
  } else {
    fs <- character(0)
  }
  unique(fs)
}

read_dea_sig <- function(files, padj_cut, lfc_cut) {
  if (length(files) == 0) {
    stop("No DEA files found. Expected files like results/DEA/DEA_DESeq2_<grp>_vs_<Control>.csv")
  }
  sig <- list()
  for (f in files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    # Expect columns: miRNA, log2FoldChange, padj (as created by run_dea.R)
    needed <- c("miRNA", "log2FoldChange", "padj")
    if (!all(needed %in% names(df))) next
    df2 <- df %>%
      filter(!is.na(padj), padj <= padj_cut, !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cut) %>%
      mutate(source_file = basename(f))
    if (nrow(df2)) sig[[length(sig) + 1]] <- df2
  }
  if (!length(sig)) return(tibble(miRNA = character(), padj = numeric(), log2FoldChange = numeric(), source_file = character()))
  bind_rows(sig)
}

rank_candidates <- function(df) {
  # rank by min padj across contrasts, then by abs LFC (desc)
  df %>%
    group_by(miRNA) %>%
    summarize(min_padj = min(padj, na.rm = TRUE),
              max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
              n_contrasts = dplyr::n_distinct(source_file),
              .groups = "drop") %>%
    arrange(min_padj, desc(max_abs_lfc))
}

chunk_vec <- function(x, n) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / n))
}

# --- Load DEA and select miRNAs ---
dea_files <- list_dea_files(dea_input)
log_line(sprintf("DEA files found: %d", length(dea_files)))
if (!length(dea_files)) {
  stop(sprintf("No DEA files available under %s", dea_input))
}

dea_sig <- read_dea_sig(dea_files, padj_cut, lfc_cut)
log_line(sprintf("Significant rows across contrasts: %d", nrow(dea_sig)))

if (nrow(dea_sig) == 0) {
  stop("No significant miRNAs found with current thresholds.")
}

miR_rank <- rank_candidates(dea_sig)
miR_list <- head(miR_rank$miRNA, top_n)
log_line(sprintf("Selected miRNAs (top_n=%d): %s", length(miR_list), paste(miR_list, collapse = ", ")))

# --- Query multiMiR (validated + predicted) ---
fetch_multimir <- function(mirnas, table_type = c("validated", "predicted"), org = "hsa", chunk_size = 20) {
  table_type <- match.arg(table_type)
  chunks <- chunk_vec(mirnas, chunk_size)
  all_res <- list()
  for (i in seq_along(chunks)) {
    m <- chunks[[i]]
    log_line(sprintf("multiMiR query %s [%d/%d]: n=%d", table_type, i, length(chunks), length(m)))
    # Try with org=; keep legacy.out=FALSE so column names are modern
    res <- tryCatch(
      get_multimir(mirna = m, org = org, table = table_type, legacy.out = FALSE, summary = TRUE),
      error = function(e) {
        log_line(sprintf("WARN: get_multimir failed for chunk %d: %s", i, e$message))
        NULL
      }
    )
    if (is.null(res)) next
    df <- tryCatch(multimir_results(res), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) next
    all_res[[length(all_res) + 1]] <- df
  }
  if (!length(all_res)) return(NULL)
  bind_rows(all_res)
}

validated_df <- fetch_multimir(miR_list, table_type = "validated", org = org_code, chunk_size = 20)
predicted_df <- fetch_multimir(miR_list, table_type = "predicted", org = org_code, chunk_size = 20)

# --- Normalize columns and write outputs ---
norm_validated <- function(df) {
  if (is.null(df) || !nrow(df)) return(tibble())
  # Harmonize common columns if present
  nm <- names(df)
  tibble(
    mirna            = df[[if ("mature_mirna_id" %in% nm) "mature_mirna_id" else "mirna"]],
    target_symbol    = if ("target_symbol" %in% nm) df[["target_symbol"]] else NA_character_,
    target_entrez    = if ("target_entrez" %in% nm) df[["target_entrez"]] else NA_character_,
    target_ensembl   = if ("target_ensembl" %in% nm) df[["target_ensembl"]] else NA_character_,
    database         = if ("database" %in% nm) df[["database"]] else NA_character_,
    support_type     = if ("support_type" %in% nm) df[["support_type"]] else if ("evidence" %in% nm) df[["evidence"]] else NA_character_,
    experiment       = if ("experiment" %in% nm) df[["experiment"]] else NA_character_,
    pmid             = if ("pmid" %in% nm) df[["pmid"]] else if ("pubmed_id" %in% nm) df[["pubmed_id"]] else NA_character_
  ) %>%
    distinct() %>%
    arrange(mirna, target_symbol, database)
}

norm_predicted <- function(df) {
  if (is.null(df) || !nrow(df)) return(tibble())
  nm <- names(df)
  # Predicted can have different scoring columns; we try to keep a generic score/rank if present.
  score_col <- c("score", "predicted_score", "prediction_score", "rank", "target_rank")
  score_pick <- score_col[score_col %in% nm]
  tibble(
    mirna            = df[[if ("mature_mirna_id" %in% nm) "mature_mirna_id" else "mirna"]],
    target_symbol    = if ("target_symbol" %in% nm) df[["target_symbol"]] else NA_character_,
    target_entrez    = if ("target_entrez" %in% nm) df[["target_entrez"]] else NA_character_,
    target_ensembl   = if ("target_ensembl" %in% nm) df[["target_ensembl"]] else NA_character_,
    database         = if ("database" %in% nm) df[["database"]] else NA_character_,
    score_or_rank    = if (length(score_pick)) df[[score_pick[1]]] else NA
  ) %>%
    distinct() %>%
    arrange(mirna, target_symbol, database)
}

val_tbl  <- norm_validated(validated_df)
pred_tbl <- norm_predicted(predicted_df)

# Create out paths
out_validated <- file.path(OUT_DIR, "targets_validated.csv")
out_predicted <- file.path(OUT_DIR, "targets_predicted.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Write CSVs
readr::write_csv(val_tbl, out_validated, na = "")
readr::write_csv(pred_tbl, out_predicted, na = "")

log_line(sprintf("Saved: %s (rows=%d)", out_validated, nrow(val_tbl)))
log_line(sprintf("Saved: %s (rows=%d)", out_predicted, nrow(pred_tbl)))

# --- Summaries (optional but handy) ---
if (nrow(val_tbl)) {
  summary_miR <- val_tbl %>%
    count(mirna, name = "validated_targets") %>%
    arrange(desc(validated_targets))
  readr::write_csv(summary_miR, file.path(OUT_DIR, "targets_validated_summary_by_miRNA.csv"))
  log_line("Saved: targets_validated_summary_by_miRNA.csv")
}

if (nrow(val_tbl)) {
  genes_validated <- val_tbl %>%
    filter(!is.na(target_symbol), target_symbol != "") %>%
    distinct(target_symbol) %>%
    arrange(target_symbol)
  readr::write_csv(genes_validated, file.path(OUT_DIR, "validated_target_genes_unique.csv"))
  log_line("Saved: validated_target_genes_unique.csv")
}

log_line("Script finalizado correctamente.")
sink()  # close log tee
