#!/usr/bin/env Rscript
# ===========================================
# Script: 17_multimir_targets_baseline.R
# Description: BASELINE companion to 16_multimir_targets.R.
#              Fetches validated and predicted miRNA->target interactions for
#              *every* miRNA in the DEA output (NO padj / |log2FC| / S_score
#              cutoff), so the IBD-overlap comparison has an all-miRNA reference.
#
#              "Baseline" = all miRNAs DESeq2 actually TESTED (rows in DEA_*.csv),
#              NOT the raw miRge3.0 count matrix. Low-count miRNAs that DESeq2
#              filtered out are intentionally excluded — they were never reliably
#              measured and would only pad the reference with noise.
#
#              Outputs (multimir_baseline_outputs/):
#                targets_validated_baseline.csv   all-miRNA validated targets
#                targets_predicted_baseline.csv   all-miRNA predicted targets
#                baseline_mirna_list.csv          the miRNAs queried + provenance
# ===========================================

suppressWarnings(suppressMessages({
  library(utils)
  library(tools)
}))

# --- Simple logging ---
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- function(...) cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = "")))

# --- Paths ---
RUN_DIR <- getwd()
OUT_DIR <- file.path(RUN_DIR, "multimir_baseline_outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LOGFILE <- file.path(OUT_DIR, "multimir_baseline.log")
sink(LOGFILE, append = FALSE, split = TRUE)

# --- Args ---
# Only two args matter for the baseline: where the DEA lives, and the organism.
# There are deliberately NO significance cutoffs here.
args <- commandArgs(trailingOnly = TRUE)
dea_input <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else RUN_DIR
org_code  <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else "hsa"

log_line("Starting 17_multimir_targets_baseline.R")
log_line(sprintf("dea_input=%s | org=%s | NO significance cutoff (all DEA miRNAs)",
                 dea_input, org_code))

# --- Libraries ---
suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(multiMiR)
}))

# --- Helpers (kept identical to 16_multimir_targets.R for consistency) ---
list_dea_files <- function(path_or_file) {
  if (file.exists(path_or_file) && file.info(path_or_file)$isdir) {
    fs <- list.files(path_or_file, pattern = "^DEA_.*\\.csv$", full.names = TRUE, recursive = TRUE)
  } else if (file.exists(path_or_file)) {
    fs <- path_or_file
  } else {
    fs <- character(0)
  }
  unique(fs)
}

# Read EVERY miRNA tested in the DEA output — no filtering on padj/lfc/S_score.
read_dea_all <- function(files) {
  if (length(files) == 0) stop("No DEA files found matching the pattern.")
  all <- list()
  for (f in files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || !("miRNA" %in% names(df))) next
    df2 <- df %>%
      filter(!is.na(miRNA), nchar(trimws(miRNA)) > 0) %>%
      mutate(source_file = basename(f))
    if (nrow(df2)) all[[length(all) + 1]] <- df2[, c("miRNA", "source_file")]
  }
  if (!length(all)) return(tibble(miRNA = character(), source_file = character()))
  bind_rows(all)
}

chunk_vec <- function(x, n) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / n))
}

# Split grouped miRNAs ("hsa-miR-a/b") into separate ids (The Smart Slash Fix).
clean_mirnas <- function(mir_vec) {
  res <- character()
  for (m in mir_vec) {
    parts <- unlist(strsplit(m, "/"))
    if (length(parts) == 1) {
      res <- c(res, parts)
    } else {
      first_part <- parts[1]
      res <- c(res, first_part)
      prefix <- ""
      if (grepl("^[a-z]{3}-[a-zA-Z]+-", first_part)) {
        prefix <- sub("^([a-z]{3}-[a-zA-Z]+-).*", "\\1", first_part)
      }
      for (i in 2:length(parts)) {
        p <- parts[i]
        if (!grepl("^[a-z]{3}-", p)) {
          res <- c(res, paste0(prefix, p))
        } else {
          res <- c(res, p)
        }
      }
    }
  }
  unique(res)
}

# --- Load DEA and take ALL miRNAs ---
dea_files <- list_dea_files(dea_input)
log_line(sprintf("DEA files found: %d", length(dea_files)))
if (!length(dea_files)) stop(sprintf("No DEA files available under %s", dea_input))

dea_all <- read_dea_all(dea_files)
log_line(sprintf("Total miRNA rows across DEA files: %d", nrow(dea_all)))

miR_all_raw <- unique(dea_all$miRNA)
log_line(sprintf("Unique miRNAs in DEA output (baseline universe): %d", length(miR_all_raw)))
if (!length(miR_all_raw)) stop("No miRNAs found in DEA output.")

miR_list_clean <- clean_mirnas(miR_all_raw)
log_line(sprintf("Cleaned/Un-slashed miRNAs for multiMiR query: %d", length(miR_list_clean)))

# Record exactly which miRNAs were queried (provenance for the comparison step).
readr::write_csv(
  tibble(mirna = miR_list_clean, source = "all DEA-tested miRNAs (no cutoff)"),
  file.path(OUT_DIR, "baseline_mirna_list.csv")
)

# --- Query multiMiR (validated + predicted) ---
fetch_multimir <- function(mirnas, table_type = c("validated", "predicted"), org = "hsa", chunk_size = 20) {
  table_type <- match.arg(table_type)
  chunks <- chunk_vec(mirnas, chunk_size)
  all_res <- list()
  for (i in seq_along(chunks)) {
    m <- chunks[[i]]
    log_line(sprintf("multiMiR query %s [%d/%d]: n=%d", table_type, i, length(chunks), length(m)))
    res <- tryCatch(
      get_multimir(mirna = m, org = org, table = table_type, legacy.out = FALSE, summary = TRUE),
      error = function(e) {
        log_line(sprintf("WARN: get_multimir failed for chunk %d: %s", i, e$message))
        NULL
      }
    )
    if (is.null(res)) next
    df <- tryCatch(res@data, error = function(e) {
      log_line(sprintf("WARN: Failed to extract @data from multiMiR object: %s", e$message))
      NULL
    })
    if (is.null(df) || !nrow(df)) next
    all_res[[length(all_res) + 1]] <- df
  }
  if (!length(all_res)) return(NULL)
  bind_rows(all_res)
}

validated_df <- fetch_multimir(miR_list_clean, table_type = "validated", org = org_code, chunk_size = 20)
predicted_df <- fetch_multimir(miR_list_clean, table_type = "predicted", org = org_code, chunk_size = 20)

# --- Normalize columns (identical schema to 16_multimir_targets.R) ---
norm_validated <- function(df) {
  if (is.null(df) || !nrow(df)) return(tibble())
  nm <- names(df)
  tibble(
    mirna          = df[[if ("mature_mirna_id" %in% nm) "mature_mirna_id" else "mirna"]],
    target_symbol  = if ("target_symbol" %in% nm) df[["target_symbol"]] else NA_character_,
    target_entrez  = if ("target_entrez" %in% nm) df[["target_entrez"]] else NA_character_,
    target_ensembl = if ("target_ensembl" %in% nm) df[["target_ensembl"]] else NA_character_,
    database       = if ("database" %in% nm) df[["database"]] else NA_character_,
    support_type   = if ("support_type" %in% nm) df[["support_type"]] else if ("evidence" %in% nm) df[["evidence"]] else NA_character_,
    experiment     = if ("experiment" %in% nm) df[["experiment"]] else NA_character_,
    pmid           = if ("pmid" %in% nm) df[["pmid"]] else if ("pubmed_id" %in% nm) df[["pubmed_id"]] else NA_character_
  ) %>% distinct() %>% arrange(mirna, target_symbol, database)
}

norm_predicted <- function(df) {
  if (is.null(df) || !nrow(df)) return(tibble())
  nm <- names(df)
  score_col  <- c("score", "predicted_score", "prediction_score", "rank", "target_rank")
  score_pick <- score_col[score_col %in% nm]
  tibble(
    mirna          = df[[if ("mature_mirna_id" %in% nm) "mature_mirna_id" else "mirna"]],
    target_symbol  = if ("target_symbol" %in% nm) df[["target_symbol"]] else NA_character_,
    target_entrez  = if ("target_entrez" %in% nm) df[["target_entrez"]] else NA_character_,
    target_ensembl = if ("target_ensembl" %in% nm) df[["target_ensembl"]] else NA_character_,
    database       = if ("database" %in% nm) df[["database"]] else NA_character_,
    score_or_rank  = if (length(score_pick)) df[[score_pick[1]]] else NA
  ) %>% distinct() %>% arrange(mirna, target_symbol, database)
}

val_tbl  <- norm_validated(validated_df)
pred_tbl <- norm_predicted(predicted_df)

out_validated <- file.path(OUT_DIR, "targets_validated_baseline.csv")
out_predicted <- file.path(OUT_DIR, "targets_predicted_baseline.csv")
readr::write_csv(val_tbl, out_validated, na = "")
readr::write_csv(pred_tbl, out_predicted, na = "")

log_line(sprintf("Saved: %s (rows=%d, miRNAs=%d)",
                 out_validated, nrow(val_tbl),
                 if (nrow(val_tbl)) dplyr::n_distinct(val_tbl$mirna) else 0))
log_line(sprintf("Saved: %s (rows=%d, miRNAs=%d)",
                 out_predicted, nrow(pred_tbl),
                 if (nrow(pred_tbl)) dplyr::n_distinct(pred_tbl$mirna) else 0))

log_line("17_multimir_targets_baseline.R finished successfully.")
sink()  # close log tee
