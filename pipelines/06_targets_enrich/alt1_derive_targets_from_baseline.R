#!/usr/bin/env Rscript
# =============================================================================
# Script: derive_targets_from_baseline.R
# Purpose: Produce the same outputs as 16_multimir_targets.R WITHOUT re-querying
#          multiMiR. The baseline run (17_multimir_targets_baseline.R) already
#          fetched targets for EVERY miRNA with no filter, so the DE-significant
#          subset is just a filter of those baseline tables — instant, offline,
#          no risk of overloading / being throttled by the multiMiR server.
#
# It writes, into ./multimir_outputs/ (identical schema to 16_multimir_targets.R):
#   targets_validated.csv
#   targets_predicted.csv
#   targets_validated_summary_by_miRNA.csv
#   validated_target_genes_unique.csv
#
# Run it from a dataset root (cwd must contain DEA_results/ and
# multimir_baseline_outputs/). Same first three args as script 14:
#   Rscript derive_targets_from_baseline.R [dea_csv_or_dir] [padj] [lfc]
#   defaults: dea_csv_or_dir = ./DEA_results , padj = 0.05 , lfc = 0.58
#
# NOTE: assumes the baseline was generated from the CURRENT DEA output. If you
#       re-ran DESeq2 since the last baseline, regenerate the baseline once
#       (17_multimir_targets_baseline.R) before using this.
# =============================================================================

suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
}))

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line  <- function(...) cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = "")))

# --- Args ---
args      <- commandArgs(trailingOnly = TRUE)
RUN_DIR   <- getwd()
dea_input <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else file.path(RUN_DIR, "DEA_results")
padj_cut  <- if (length(args) >= 2 && nzchar(args[[2]])) as.numeric(args[[2]]) else 0.05
lfc_cut   <- if (length(args) >= 3 && nzchar(args[[3]])) as.numeric(args[[3]]) else 0.58

OUT_DIR      <- file.path(RUN_DIR, "multimir_outputs")
BASELINE_DIR <- file.path(RUN_DIR, "multimir_baseline_outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_line("derive_targets_from_baseline.R starting")
log_line(sprintf("dea_input=%s | padj<=%.3f | |log2FC|>=%.2f", dea_input, padj_cut, lfc_cut))

# --- Helpers (kept identical to 16_multimir_targets.R) ---
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

read_dea_sig <- function(files, padj_cut, lfc_cut) {
  if (length(files) == 0) stop("No DEA files found matching the pattern.")
  sig <- list()
  for (f in files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    needed <- c("miRNA", "log2FoldChange", "padj")
    if (!all(needed %in% names(df))) next
    df2 <- df %>% filter(!is.na(padj), padj <= padj_cut,
                         !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cut)
    if (nrow(df2)) sig[[length(sig) + 1]] <- df2
  }
  if (!length(sig)) return(character(0))
  unique(bind_rows(sig)$miRNA)
}

# The Smart Slash Fix — identical to 16_multimir_targets.R so the names match the
# baseline tables (the baseline applied the same cleaning before querying).
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
        if (!grepl("^[a-z]{3}-", p)) res <- c(res, paste0(prefix, p)) else res <- c(res, p)
      }
    }
  }
  unique(res)
}

# --- 1. DE-significant miRNAs ---
dea_files <- list_dea_files(dea_input)
log_line(sprintf("DEA files found: %d", length(dea_files)))
if (!length(dea_files)) stop(sprintf("No DEA files under %s", dea_input))

sig_raw   <- read_dea_sig(dea_files, padj_cut, lfc_cut)
if (!length(sig_raw)) stop("No significant miRNAs at the given thresholds.")
sig_clean <- clean_mirnas(sig_raw)
log_line(sprintf("DE-significant miRNAs: %d (cleaned/un-slashed: %d)", length(sig_raw), length(sig_clean)))

# --- 2. Load baseline tables ---
val_base_path  <- file.path(BASELINE_DIR, "targets_validated_baseline.csv")
pred_base_path <- file.path(BASELINE_DIR, "targets_predicted_baseline.csv")
if (!file.exists(val_base_path) && !file.exists(pred_base_path)) {
  stop(sprintf("No baseline tables in %s — run 17_multimir_targets_baseline.R first.", BASELINE_DIR))
}
read_base <- function(p) if (file.exists(p)) suppressMessages(read_csv(p, col_types = cols(.default = "c"))) else NULL
val_base  <- read_base(val_base_path)
pred_base <- read_base(pred_base_path)

# --- 3. Filter baseline to the DE-significant miRNAs ---
matched <- intersect(sig_clean, unique(c(val_base$mirna, pred_base$mirna)))
missing <- setdiff(sig_clean, matched)
log_line(sprintf("miRNAs matched in baseline: %d / %d", length(matched), length(sig_clean)))
if (length(missing)) {
  log_line(sprintf("WARN: %d DE-significant miRNAs absent from baseline (may have no targets, or baseline is stale): %s",
                   length(missing), paste(head(missing, 20), collapse = ", ")))
}

val_cols  <- c("mirna","target_symbol","target_entrez","target_ensembl","database","support_type","experiment","pmid")
pred_cols <- c("mirna","target_symbol","target_entrez","target_ensembl","database","score_or_rank")

val_tbl <- if (!is.null(val_base)) {
  val_base %>% filter(mirna %in% sig_clean) %>%
    select(any_of(val_cols)) %>% distinct() %>% arrange(mirna, target_symbol, database)
} else tibble()

pred_tbl <- if (!is.null(pred_base)) {
  pred_base %>% filter(mirna %in% sig_clean) %>%
    select(any_of(pred_cols)) %>% distinct() %>% arrange(mirna, target_symbol, database)
} else tibble()

# --- 4. Write outputs (identical names/schema to 16_multimir_targets.R) ---
write_csv(val_tbl,  file.path(OUT_DIR, "targets_validated.csv"), na = "")
write_csv(pred_tbl, file.path(OUT_DIR, "targets_predicted.csv"), na = "")
log_line(sprintf("Saved targets_validated.csv (rows=%d) and targets_predicted.csv (rows=%d)",
                 nrow(val_tbl), nrow(pred_tbl)))

if (nrow(val_tbl)) {
  val_tbl %>% count(mirna, name = "validated_targets") %>% arrange(desc(validated_targets)) %>%
    write_csv(file.path(OUT_DIR, "targets_validated_summary_by_miRNA.csv"))
  val_tbl %>% filter(!is.na(target_symbol), target_symbol != "") %>%
    distinct(target_symbol) %>% arrange(target_symbol) %>%
    write_csv(file.path(OUT_DIR, "validated_target_genes_unique.csv"))
  log_line("Saved targets_validated_summary_by_miRNA.csv and validated_target_genes_unique.csv")
}

log_line("Done — derived from baseline, no multiMiR query performed.")
