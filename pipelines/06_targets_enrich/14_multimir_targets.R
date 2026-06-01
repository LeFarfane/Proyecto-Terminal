#!/usr/bin/env Rscript
# ===========================================
# Script: multimir_targets.R
# Description: Fetches validated and predicted miRNA→target interactions for human
#              using multiMiR, based on significant miRNAs from DEA results.
#              Ranks by the custom S_miRNA score, dynamically splits grouped miRNAs,
#              queries in chunks, and writes final tables to a dedicated output folder.
# ===========================================

suppressWarnings(suppressMessages({
  library(utils)
  library(tools)
}))

# --- Simple logging ---
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- function(...) cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = "")))

# --- Paths (PORTABLE FIX + DEDICATED FOLDER) ---
RUN_DIR <- getwd()
OUT_DIR <- file.path(RUN_DIR, "multimir_outputs")

# Create the output directory automatically
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Put the log file inside the new folder too
LOGFILE <- file.path(OUT_DIR, "multimir.log")
sink(LOGFILE, append = FALSE, split = TRUE)

# --- Args ---
args <- commandArgs(trailingOnly = TRUE)
dea_input    <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else RUN_DIR
padj_cut     <- if (length(args) >= 2 && nzchar(args[[2]])) as.numeric(args[[2]]) else 0.05
lfc_cut      <- if (length(args) >= 3 && nzchar(args[[3]])) as.numeric(args[[3]]) else 0.58
s_score_cut  <- if (length(args) >= 4 && nzchar(args[[4]])) as.numeric(args[[4]]) else 3.0
org_code     <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else "hsa" 

log_line("Starting multimir_targets.R")
log_line(sprintf("dea_input=%s | padj_cut=%.3f | lfc_cut=%.2f | s_score_cut=%.2f | org=%s",
                 dea_input, padj_cut, lfc_cut, s_score_cut, org_code))

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
    fs <- list.files(path_or_file, pattern = "^DEA_.*\\.csv$", full.names = TRUE, recursive = TRUE)
  } else if (file.exists(path_or_file)) {
    fs <- path_or_file
  } else {
    fs <- character(0)
  }
  unique(fs)
}

read_dea_sig <- function(files, padj_cut, lfc_cut) {
  if (length(files) == 0) {
    stop("No DEA files found matching the pattern.")
  }
  sig <- list()
  for (f in files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    
    needed <- c("miRNA", "log2FoldChange", "padj", "S_miRNA")
    if (!all(needed %in% names(df))) next
    
    df2 <- df %>%
      filter(!is.na(padj), padj <= padj_cut, !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cut) %>%
      mutate(source_file = basename(f))
    if (nrow(df2)) sig[[length(sig) + 1]] <- df2
  }
  if (!length(sig)) return(tibble(miRNA = character(), padj = numeric(), log2FoldChange = numeric(), S_miRNA = numeric(), source_file = character()))
  bind_rows(sig)
}

rank_candidates <- function(df) {
  df %>%
    group_by(miRNA) %>%
    summarize(
      max_S_score = max(abs(S_miRNA), na.rm = TRUE), 
      min_padj = min(padj, na.rm = TRUE),
      max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
      n_contrasts = dplyr::n_distinct(source_file),
      .groups = "drop"
    ) %>%
    arrange(desc(max_S_score))
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

# Apply dynamic threshold instead of top_n
miR_filtered <- miR_rank %>% filter(max_S_score >= s_score_cut)

if (nrow(miR_filtered) == 0) {
  stop(sprintf("No miRNAs passed the dynamic threshold of S_score >= %.2f", s_score_cut))
}

miR_list <- miR_filtered$miRNA
log_line(sprintf("Selected miRNAs (S_score >= %.2f): %d candidates", s_score_cut, length(miR_list)))

# --- Split grouped miRNAs (The Smart Slash Fix) ---
clean_mirnas <- function(mir_vec) {
  res <- character()
  for (m in mir_vec) {
    parts <- unlist(strsplit(m, "/"))
    if (length(parts) == 1) {
      res <- c(res, parts)
    } else {
      first_part <- parts[1]
      res <- c(res, first_part)
      
      # Extract prefix (e.g., "hsa-miR-" or "hsa-let-")
      prefix <- ""
      if (grepl("^[a-z]{3}-[a-zA-Z]+-", first_part)) {
        prefix <- sub("^([a-z]{3}-[a-zA-Z]+-).*", "\\1", first_part)
      }
      
      # Prepend prefix to subsequent parts if they don't have one
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

miR_list_clean <- clean_mirnas(miR_list)
log_line(sprintf("Cleaned/Un-slashed miRNAs for multiMiR query: %d", length(miR_list_clean)))

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

# --- Normalize columns and write outputs ---
norm_validated <- function(df) {
  if (is.null(df) || !nrow(df)) return(tibble())
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

# Write CSVs directly to the newly created multimir_outputs folder
out_validated <- file.path(OUT_DIR, "targets_validated.csv")
out_predicted <- file.path(OUT_DIR, "targets_predicted.csv")

readr::write_csv(val_tbl, out_validated, na = "")
readr::write_csv(pred_tbl, out_predicted, na = "")

log_line(sprintf("Saved: %s (rows=%d)", out_validated, nrow(val_tbl)))
log_line(sprintf("Saved: %s (rows=%d)", out_predicted, nrow(pred_tbl)))

# --- Query multiMiR for Diseases (IBD, UC, CD) ---
fetch_multimir_disease <- function(mirnas, disease_names, org = "hsa", chunk_size = 20) {
  chunks <- chunk_vec(mirnas, chunk_size)
  all_res <- list()
  
  for (d in disease_names) {
    for (i in seq_along(chunks)) {
      m <- chunks[[i]]
      log_line(sprintf("multiMiR query disease.drug: '%s' [%d/%d]: n=%d", d, i, length(chunks), length(m)))
      res <- tryCatch(
        get_multimir(mirna = m, org = org, disease.drug = d, table = "disease.drug", legacy.out = FALSE, summary = TRUE),
        error = function(e) {
          log_line(sprintf("WARN: get_multimir failed for disease '%s' (chunk %d): %s", d, i, e$message))
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
  }
  if (!length(all_res)) return(NULL)
  dplyr::bind_rows(all_res) %>% dplyr::distinct()
}

diseases_to_query <- list(
  "IBD" = c("inflammatory bowel disease", "IBD"),
  "UC"  = c("ulcerative colitis", "UC"),
  "CD"  = c("crohn disease", "CD")
)

for (dz_key in names(diseases_to_query)) {
  dz_names <- diseases_to_query[[dz_key]]
  # Query strictly constrained to the dynamically selected miRNAs
  dz_res <- fetch_multimir_disease(mirnas = miR_list_clean, disease_names = dz_names, org = org_code, chunk_size = 20)
  
  if (!is.null(dz_res) && nrow(dz_res) > 0) {
    out_file <- file.path(OUT_DIR, sprintf("targets_disease_%s.csv", dz_key))
    readr::write_csv(dz_res, out_file, na = "")
    log_line(sprintf("Saved: %s (rows=%d)", out_file, nrow(dz_res)))
  } else {
    log_line(sprintf("No multiMiR disease/drug results found for: %s among selected miRNAs.", dz_key))
  }
}

# --- Summaries ---
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

log_line("Script finished successfully.")
sink()  # close log tee