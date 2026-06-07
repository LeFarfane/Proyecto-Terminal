#!/usr/bin/env Rscript
# =============================================================================
# Script : 07_networks_clin/clinical_roc.R
# Purpose: Biomarker-performance evaluation for candidate miRNAs.
#          Replaces clinical_eval.R — wired to the actual DEA_results_round_2
#          outputs instead of the old generic paths.
#
#  Two contrasts
#  ─────────────
#  1. R_Baseline vs NR_Baseline  ← primary clinical question:
#       can baseline miRNA expression predict glucocorticoid response?
#  2. UC_Baseline (R+NR) vs Healthy Controls
#       ← classic disease vs control
#
#  Outputs  (all under DEA_results_round_2/clinical_outputs/)
#  ──────────────────────────────────────────────────────────
#  roc_R_vs_NR.csv            per-miRNA AUC for contrast 1
#  roc_UC_vs_HC.csv           per-miRNA AUC for contrast 2
#  cv_model_R_vs_NR.csv       5-fold CV logistic (top 5 by AUC)
#  candidate_mirnas.csv       the 27 candidate miRNAs with direction
#  figures/ROC_RvsNR_*.png    top-6 ROC curves, contrast 1
#  figures/ROC_UCvsHC_*.png   top-6 ROC curves, contrast 2
#  figures/ROC_multi_RvsNR.png multi-miRNA CV model
#  figures/violin_*.png       3-group violin (HC / NR_base / R_base)
#
# Usage:  Rscript clinical_roc.R
#         (no arguments needed — all paths are hardcoded relative to ROOT_DIR)
# =============================================================================

suppressWarnings(suppressMessages({
  library(utils)
  library(tools)
}))

# ── logging ──────────────────────────────────────────────────────────────────
ts      <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) {
  m <- paste0("[", ts(), "] ", paste0(..., collapse = ""))
  cat(m, "\n")
  lf <- getOption("roc_logfile")
  if (!is.null(lf)) cat(m, "\n", file = lf, append = TRUE)
}

# ── root dir (this script lives two levels below Pipelines/) ─────────────────
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(getwd())
}
SCRIPT_DIR <- tryCatch(get_script_dir(), error = function(e) normalizePath(getwd()))
ROOT_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

# ── paths — resolved from optional arg or auto-discovered ────────────────────
# Usage:  Rscript 20_treatment_response_roc.R [DEA_results_round_2_path]
# If no argument is given the script searches under ROOT_DIR/outputs/ for a
# DEA_results_round_2 directory that contains the required files.
args_trail <- commandArgs(trailingOnly = TRUE)

find_dea_base <- function(search_root) {
  # Returns all DEA_results_round_2 dirs under search_root that contain the
  # required VST matrix (proof the full pipeline ran there).
  hits <- list.files(search_root, pattern = "^DEA_results_round_2$",
                     recursive = TRUE, full.names = TRUE, include.dirs = TRUE)
  hits <- hits[file.info(hits)$isdir %in% TRUE]
  # Keep only those with the VST matrix present
  hits[sapply(hits, function(h)
    file.exists(file.path(h, "figures_and_QC", "vst_normalized_matrix.csv")))]
}

if (length(args_trail) >= 1 && nzchar(args_trail[1])) {
  DEA_BASE <- normalizePath(args_trail[1], mustWork = FALSE)
} else {
  # 1st priority: search from the current working directory (user is already
  #               inside the right dataset folder when they run the script)
  hits <- find_dea_base(getwd())

  # 2nd priority: search the whole outputs/ tree
  if (length(hits) == 0) hits <- find_dea_base(file.path(ROOT_DIR, "outputs"))

  if (length(hits) == 0)
    stop("Cannot find a valid DEA_results_round_2 (with vst_normalized_matrix.csv). ",
         "Pass the path explicitly as argument 1.")
  if (length(hits) > 1)
    cat(sprintf("[WARN] Multiple valid DEA_results_round_2 found — using: %s\n", hits[1]))

  DEA_BASE <- hits[1]
}

# Derive sibling paths from DEA_BASE
MUESTRAS_DIR <- dirname(DEA_BASE)                     # …/muestras_clasificadas
VST_PATH     <- file.path(DEA_BASE,      "figures_and_QC", "vst_normalized_matrix.csv")
META_PATH    <- file.path(MUESTRAS_DIR,  "sample_metadata.csv")
DEA_RESP     <- file.path(DEA_BASE, "DEA_Default_Replace", "DEA_Default_UC_Responder.csv")
DEA_NONR     <- file.path(DEA_BASE, "DEA_Default_Replace", "DEA_Default_UC_NonResponder.csv")

OUT_DIR     <- file.path(DEA_BASE, "clinical_outputs")
FIG_DIR     <- file.path(OUT_DIR, "figures")
LOG_DIR     <- file.path(OUT_DIR, "logs")

for (d in c(OUT_DIR, FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOG_DIR, "clinical_roc.log")
options(roc_logfile = LOGFILE)
write("", file = LOGFILE)

# ── thresholds ────────────────────────────────────────────────────────────────
PADJ_THRESH     <- 0.05
LFC_THRESH      <- 0.58
TOP_ROC_PLOTS   <- 6
TOP_VIOLIN      <- 12
TOP_MODEL       <- 5
CV_FOLDS        <- 5
set.seed(42)

cat("\n=====================================================\n")
cat("  NOTE: This script only applies to datasets with\n")
cat("  treatment-response data (responders vs non-responders\n")
cat("  at baseline). Auto-discovers DEA_results_round_2 under\n")
cat("  outputs/, or pass the path as argument 1.\n")
cat("  Skip this step for simple case-control datasets.\n")
cat("=====================================================\n\n")

log_msg("clinical_roc.R starting")
log_msg("ROOT_DIR : ", ROOT_DIR)
log_msg("DEA_BASE : ", DEA_BASE)

# ── libraries ─────────────────────────────────────────────────────────────────
suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pROC)
}))

# ── helpers ───────────────────────────────────────────────────────────────────

# Parse the verbose group descriptions in sample_metadata into short labels.
# Groups in PRJNA471862:
#   HC          : healthy controls
#   R_Baseline  : R patients, sample taken BEFORE glucocorticoids
#   NR_Baseline : NR patients, sample taken BEFORE glucocorticoids
#   R_Day3      : R patients, sample on day 3 of treatment
#   NR_Day3     : NR patients, sample on day 3 of treatment
parse_group <- function(x) {
  x_lo <- tolower(x)
  dplyr::case_when(
    grepl("healthy", x_lo) ~ "HC",
    grepl("poor response|non-responder", x_lo) & grepl("before", x_lo) ~ "NR_Baseline",
    (grepl("good response|responders, r", x_lo) | grepl("tand who have", x_lo)) &
      grepl("before", x_lo) ~ "R_Baseline",
    grepl("3rd day|third day", x_lo) & grepl("non-responder", x_lo) ~ "NR_Day3",
    grepl("3rd day|third day", x_lo) ~ "R_Day3",
    TRUE ~ "Other"
  )
}

# Clean miRNA name for use as a filename
safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# Build a ROC plot with ggplot2 and save it
save_roc_plot <- function(roc_obj, title_str, out_path) {
  df <- data.frame(
    fpr  = 1 - roc_obj$specificities,
    tpr  = roc_obj$sensitivities
  )
  p <- ggplot(df, aes(x = fpr, y = tpr)) +
    geom_line(color = "#1f77b4", linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    coord_equal() +
    labs(title = title_str, x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 12)
  ggsave(out_path, p, width = 5, height = 5, dpi = 300)
}

# Compute per-miRNA ROC/AUC for a binary outcome vector
compute_roc_table <- function(expr_mat, y_factor, mirna_list) {
  results <- vector("list", length(mirna_list))
  for (i in seq_along(mirna_list)) {
    mir  <- mirna_list[i]
    xvec <- as.numeric(expr_mat[mir, ])
    keep <- is.finite(xvec) & !is.na(y_factor)
    if (sum(keep) < 6 || length(unique(as.character(y_factor[keep]))) < 2) next
    roc_obj <- tryCatch(
      pROC::roc(response  = y_factor[keep],
                predictor = xvec[keep],
                levels    = levels(y_factor),
                direction = "auto",
                quiet     = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_obj)) next
    results[[i]] <- data.frame(
      miRNA = mir,
      AUC   = as.numeric(pROC::auc(roc_obj)),
      n_pos = sum(as.numeric(y_factor[keep]) == 2, na.rm = TRUE),
      n_neg = sum(as.numeric(y_factor[keep]) == 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(out) || nrow(out) == 0) return(data.frame())
  out[order(-out$AUC), , drop = FALSE]
}

# ── 1. Load data ─────────────────────────────────────────────────────────────
for (p in c(VST_PATH, META_PATH, DEA_RESP, DEA_NONR)) {
  if (!file.exists(p)) { log_msg("ERROR: file not found: ", p); quit(status = 1) }
}

log_msg("Loading VST matrix ...")
expr <- tryCatch(read.csv(VST_PATH, check.names = FALSE, row.names = 1),
                 error = function(e) NULL)
if (is.null(expr)) { log_msg("ERROR: cannot read VST matrix."); quit(status = 1) }

log_msg("Loading metadata ...")
meta <- tryCatch(read.csv(META_PATH, stringsAsFactors = FALSE),
                 error = function(e) NULL)
if (is.null(meta)) { log_msg("ERROR: cannot read metadata."); quit(status = 1) }

log_msg("Loading DEA results ...")
dea_R  <- read.csv(DEA_RESP, stringsAsFactors = FALSE)
dea_NR <- read.csv(DEA_NONR, stringsAsFactors = FALSE)

# ── 2. Align expression matrix with metadata ──────────────────────────────────
common_ids <- intersect(colnames(expr), meta$sample_id)
log_msg(sprintf("Samples: expr=%d  meta=%d  common=%d",
                ncol(expr), nrow(meta), length(common_ids)))
if (length(common_ids) < 10) { log_msg("ERROR: too few overlapping samples."); quit(status = 1) }

expr <- expr[, common_ids, drop = FALSE]
meta <- meta[match(common_ids, meta$sample_id), , drop = FALSE]

# Parse group labels
meta$group_parsed <- parse_group(meta$group)
group_counts <- table(meta$group_parsed)
log_msg("Group counts: ", paste(names(group_counts), group_counts, sep = "=", collapse = "  "))

# ── 3. Candidate miRNA selection ──────────────────────────────────────────────
# Union of significant miRNAs from both comparisons (padj ≤ PADJ_THRESH, |log2FC| ≥ LFC_THRESH)
sig_R <- dea_R  %>% filter(!is.na(padj), padj <= PADJ_THRESH,
                             !is.na(log2FoldChange), abs(log2FoldChange) >= LFC_THRESH) %>%
  mutate(direction = if_else(log2FoldChange > 0, "UP", "DOWN"),
         min_padj  = padj,
         source    = "Responder")

sig_NR <- dea_NR %>% filter(!is.na(padj), padj <= PADJ_THRESH,
                              !is.na(log2FoldChange), abs(log2FoldChange) >= LFC_THRESH) %>%
  mutate(direction = if_else(log2FoldChange > 0, "UP", "DOWN"),
         min_padj  = padj,
         source    = "NonResponder")

# Combine — if a miRNA appears in both, take the entry with the lower padj for direction call
candidates <- bind_rows(sig_R, sig_NR) %>%
  group_by(miRNA) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(padj, desc(abs(log2FoldChange)))

# Keep only those present in the expression matrix
candidates <- candidates %>% filter(miRNA %in% rownames(expr))
n_cand <- nrow(candidates)
log_msg(sprintf("Candidate miRNAs: %d  (padj <= %.2f, |log2FC| >= %.1f, present in VST matrix)",
                n_cand, PADJ_THRESH, LFC_THRESH))

if (n_cand == 0) {
  log_msg("ERROR: no candidate miRNAs survive filters. Exiting.")
  quit(status = 1)
}

readr::write_csv(candidates %>%
                   select(miRNA, log2FoldChange, padj, direction, source),
                 file.path(OUT_DIR, "candidate_mirnas.csv"))
log_msg("Saved: candidate_mirnas.csv")

mirna_list <- candidates$miRNA

# ── 4. Contrast 1 — R_Baseline vs NR_Baseline ────────────────────────────────
log_msg("\n--- Contrast 1: R_Baseline vs NR_Baseline ---")

idx_c1 <- which(meta$group_parsed %in% c("R_Baseline", "NR_Baseline"))
meta_c1 <- meta[idx_c1, ]
expr_c1  <- expr[mirna_list, idx_c1, drop = FALSE]

# Ordered factor: NR = reference (0), R = positive (1)
y_c1 <- factor(meta_c1$group_parsed,
               levels = c("NR_Baseline", "R_Baseline"),
               labels = c("NR", "R"))

log_msg(sprintf("  NR_Baseline=%d  R_Baseline=%d",
                sum(y_c1 == "NR", na.rm = TRUE), sum(y_c1 == "R", na.rm = TRUE)))

roc_c1 <- compute_roc_table(as.matrix(expr_c1), y_c1, mirna_list)

if (nrow(roc_c1) > 0) {
  readr::write_csv(roc_c1, file.path(OUT_DIR, "roc_R_vs_NR.csv"))
  log_msg(sprintf("  Saved: roc_R_vs_NR.csv  (top AUC: %.3f for %s)",
                  roc_c1$AUC[1], roc_c1$miRNA[1]))

  # Top-N ROC plots
  for (mir in head(roc_c1$miRNA, TOP_ROC_PLOTS)) {
    xvec    <- as.numeric(expr_c1[mir, ])
    keep    <- is.finite(xvec) & !is.na(y_c1)
    roc_obj <- pROC::roc(response  = y_c1[keep], predictor = xvec[keep],
                         levels = levels(y_c1), direction = "auto", quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    title   <- sprintf("R vs NR  |  %s  (AUC=%.3f)", mir, auc_val)
    out_png <- file.path(FIG_DIR, sprintf("ROC_RvsNR_%s.png", safe_name(mir)))
    save_roc_plot(roc_obj, title, out_png)
  }
  log_msg(sprintf("  Saved %d ROC figures (contrast 1).", min(TOP_ROC_PLOTS, nrow(roc_c1))))

} else {
  log_msg("  WARN: could not compute any ROC for contrast 1.")
}

# ── 5. CV logistic model (Contrast 1) ─────────────────────────────────────────
if (nrow(roc_c1) >= 3) {
  log_msg(sprintf("  Building %d-fold CV logistic model (top %d miRNAs) ...",
                  CV_FOLDS, TOP_MODEL))

  top_mir <- head(roc_c1$miRNA, TOP_MODEL)
  X_raw   <- t(as.matrix(expr_c1[top_mir, , drop = FALSE]))
  X_sc    <- scale(X_raw)
  y_num   <- as.numeric(y_c1) - 1L          # 0 = NR, 1 = R
  df_cv   <- data.frame(y = y_num, X_sc, check.names = FALSE)

  keep_rows <- complete.cases(df_cv)
  df_cv     <- df_cv[keep_rows, , drop = FALSE]
  n_cv      <- nrow(df_cv)

  folds   <- sample(rep(seq_len(CV_FOLDS), length.out = n_cv))
  probs   <- rep(NA_real_, n_cv)
  for (fold in seq_len(CV_FOLDS)) {
    tr  <- which(folds != fold)
    te  <- which(folds == fold)
    fit <- tryCatch(glm(y ~ ., data = df_cv[tr, , drop = FALSE],
                        family = binomial()),
                    error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) next
    probs[te] <- suppressWarnings(
      predict(fit, newdata = df_cv[te, , drop = FALSE], type = "response"))
  }

  valid  <- !is.na(probs)
  if (sum(valid) >= 6 && length(unique(df_cv$y[valid])) == 2) {
    y_fac <- factor(df_cv$y[valid], levels = c(0, 1), labels = c("NR", "R"))
    roc_cv <- pROC::roc(response  = y_fac, predictor = probs[valid],
                        levels = c("NR", "R"), direction = "auto", quiet = TRUE)
    auc_cv <- as.numeric(pROC::auc(roc_cv))

    readr::write_csv(
      data.frame(model   = paste0("logit_top", TOP_MODEL),
                 AUC_CV  = auc_cv,
                 n       = sum(valid),
                 miRNAs  = paste(top_mir, collapse = ";"),
                 stringsAsFactors = FALSE),
      file.path(OUT_DIR, "cv_model_R_vs_NR.csv")
    )
    log_msg(sprintf("  Saved: cv_model_R_vs_NR.csv  (AUC_CV=%.3f, n=%d)", auc_cv, sum(valid)))

    save_roc_plot(
      roc_cv,
      sprintf("Multi-miRNA (top%d) CV logistic — R vs NR (AUC=%.3f)", TOP_MODEL, auc_cv),
      file.path(FIG_DIR, "ROC_multi_RvsNR.png")
    )
  } else {
    log_msg("  WARN: insufficient variation for CV ROC (contrast 1).")
  }
} else {
  log_msg("  WARN: fewer than 3 miRNAs with ROC — skipping CV model.")
}

# ── 6. Contrast 2 — UC_Baseline vs HC ─────────────────────────────────────────
log_msg("\n--- Contrast 2: UC_Baseline vs HC ---")

idx_c2  <- which(meta$group_parsed %in% c("R_Baseline", "NR_Baseline", "HC"))
meta_c2 <- meta[idx_c2, ]
expr_c2 <- expr[mirna_list, idx_c2, drop = FALSE]

y_c2 <- factor(
  if_else(meta_c2$group_parsed == "HC", "HC", "UC"),
  levels = c("HC", "UC")
)

log_msg(sprintf("  HC=%d  UC_Baseline=%d",
                sum(y_c2 == "HC"), sum(y_c2 == "UC")))

roc_c2 <- compute_roc_table(as.matrix(expr_c2), y_c2, mirna_list)

if (nrow(roc_c2) > 0) {
  readr::write_csv(roc_c2, file.path(OUT_DIR, "roc_UC_vs_HC.csv"))
  log_msg(sprintf("  Saved: roc_UC_vs_HC.csv  (top AUC: %.3f for %s)",
                  roc_c2$AUC[1], roc_c2$miRNA[1]))

  for (mir in head(roc_c2$miRNA, TOP_ROC_PLOTS)) {
    xvec    <- as.numeric(expr_c2[mir, ])
    keep    <- is.finite(xvec) & !is.na(y_c2)
    roc_obj <- pROC::roc(response  = y_c2[keep], predictor = xvec[keep],
                         levels = levels(y_c2), direction = "auto", quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    title   <- sprintf("UC vs HC  |  %s  (AUC=%.3f)", mir, auc_val)
    out_png <- file.path(FIG_DIR, sprintf("ROC_UCvsHC_%s.png", safe_name(mir)))
    save_roc_plot(roc_obj, title, out_png)
  }
  log_msg(sprintf("  Saved %d ROC figures (contrast 2).", min(TOP_ROC_PLOTS, nrow(roc_c2))))
} else {
  log_msg("  WARN: could not compute any ROC for contrast 2.")
}

# ── 7. Violin plots — top miRNAs, three-group comparison ─────────────────────
log_msg("\n--- Violin plots (3-group: HC / NR_Baseline / R_Baseline) ---")

# Rank miRNAs: if contrast 1 available use that, else fall back to contrast 2
rank_source <- if (nrow(roc_c1) > 0) roc_c1 else roc_c2
top_violin_list <- head(rank_source$miRNA, TOP_VIOLIN)

idx_vio  <- which(meta$group_parsed %in% c("HC", "NR_Baseline", "R_Baseline"))
meta_vio <- meta[idx_vio, ]
expr_vio <- as.matrix(expr[top_violin_list, idx_vio, drop = FALSE])

group_palette <- c(HC = "#4daf4a", NR_Baseline = "#e41a1c", R_Baseline = "#377eb8")
group_labels  <- c(HC = "Healthy\nControl",
                   NR_Baseline = "UC\nNon-Responder",
                   R_Baseline  = "UC\nResponder")

for (mir in top_violin_list) {
  df_v <- data.frame(
    expr  = as.numeric(expr_vio[mir, ]),
    group = meta_vio$group_parsed,
    stringsAsFactors = FALSE
  ) %>%
    filter(group %in% c("HC", "NR_Baseline", "R_Baseline"),
           is.finite(expr)) %>%
    mutate(group = factor(group, levels = c("HC", "NR_Baseline", "R_Baseline")))

  if (nrow(df_v) < 6) next

  p <- ggplot(df_v, aes(x = group, y = expr, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.12, outlier.size = 0.7, outlier.alpha = 0.5,
                 fill = "white", alpha = 0.8) +
    scale_fill_manual(values  = group_palette,
                      labels  = group_labels,
                      name    = NULL) +
    scale_x_discrete(labels  = group_labels) +
    labs(title = mir,
         x     = NULL,
         y     = "VST expression") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x     = element_text(size = 10))

  ggsave(
    file.path(FIG_DIR, sprintf("violin_%s.png", safe_name(mir))),
    p, width = 5.5, height = 4.5, dpi = 300
  )
}
log_msg(sprintf("Saved violin plots for %d miRNAs.", length(top_violin_list)))

# ── 8. Summary to stdout ──────────────────────────────────────────────────────
log_msg("\n========== SUMMARY ==========")
log_msg(sprintf("Candidate miRNAs evaluated  : %d", n_cand))
if (nrow(roc_c1) > 0) {
  log_msg(sprintf("Contrast 1 (R vs NR)        : top AUC = %.3f  (%s)",
                  roc_c1$AUC[1], roc_c1$miRNA[1]))
}
if (nrow(roc_c2) > 0) {
  log_msg(sprintf("Contrast 2 (UC vs HC)       : top AUC = %.3f  (%s)",
                  roc_c2$AUC[1], roc_c2$miRNA[1]))
}
cv_file <- file.path(OUT_DIR, "cv_model_R_vs_NR.csv")
if (file.exists(cv_file)) {
  cv_row <- readr::read_csv(cv_file, show_col_types = FALSE)
  log_msg(sprintf("Multi-miRNA CV-AUC (R vs NR): %.3f", cv_row$AUC_CV[1]))
}
log_msg(sprintf("Outputs saved to: %s", OUT_DIR))
log_msg("Script completed successfully.")
