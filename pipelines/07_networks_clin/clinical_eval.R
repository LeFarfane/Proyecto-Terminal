#!/usr/bin/env Rscript
# ===========================================
# Script: 07_networks_clin/clinical_eval.R
# Descripción: Clinical evaluation of candidate miRNAs:
#              (1) Correlates miRNA expression with clinical numeric variables (Spearman, BH FDR).
#              (2) Case/control (or responder) classification per-miRNA (ROC/AUC) + simple multi-miRNA model.
#              (3) Saves tables and figures under results/clinical/.
# Fecha: 2025-09-07
# ===========================================

suppressWarnings(suppressMessages({
  library(utils)
  library(tools)
}))

# --- Logging (stdout + archivo) ---
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- function(...) {
  msg <- sprintf("[%s] %s", timestamp(), paste0(..., collapse = ""))
  cat(msg, "\n")
  if (!is.null(getOption("clin_logfile"))) {
    cat(msg, "\n", file = getOption("clin_logfile"), append = TRUE)
  }
}

# --- Paths base (relativos a PT/) ---
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile))
  return(normalizePath(getwd()))
}
SCRIPT_PATH <- tryCatch(get_script_path(), error = function(e) normalizePath(getwd()))
ROOT_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), "..", ".."))  # .../PT

CLIN_DIR <- file.path(ROOT_DIR, "results", "clinical")
FIG_DIR  <- file.path(CLIN_DIR, "figures")
LOG_DIR  <- file.path(ROOT_DIR, "outputs", "logs")
dir.create(CLIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,  recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOG_DIR, "07_clinical_eval.log")
options(clin_logfile = LOGFILE)
write("", file = LOGFILE)

# --- Args ---
# uso:
#   Rscript clinical_eval.R [norm_counts_csv] [meta_csv] [mirna_list_file] [outcome_col] [positive_label] [top_n_fallback] [numeric_vars_csv]
# defaults:
#   norm_counts_csv = results/DEA/normalized_counts_deseq2.csv
#   meta_csv        = metadata/sample_metadata.csv
#   mirna_list_file = "" (si vacío, toma miRNAs sig de DEA o top var = 30)
#   outcome_col     = "" (se infiere: 'response'/'responder' o 'disease')
#   positive_label  = "" (se infiere: 'Responder' o 'Case' (no-Control))
#   top_n_fallback  = 30
#   numeric_vars_csv= "" (si vacío, detecta numeric en meta; prioriza Mayo, CDAI, CRP, calprotectina)
args <- commandArgs(trailingOnly = TRUE)
default_counts <- file.path(ROOT_DIR, "results", "DEA", "normalized_counts_deseq2.csv")
default_meta   <- file.path(ROOT_DIR, "metadata", "sample_metadata.csv")
counts_csv     <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_counts
meta_csv       <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else default_meta
mirna_list_fn  <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else ""
outcome_col    <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else ""
positive_label <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else ""
top_n_fallback <- if (length(args) >= 6 && nzchar(args[[6]])) as.integer(args[[6]]) else 30
numeric_vars_fn<- if (length(args) >= 7 && nzchar(args[[7]])) args[[7]] else ""

log_line("Starting clinical_eval.R")
log_line(sprintf("counts=%s | meta=%s | mirna_list=%s | outcome_col=%s | positive_label=%s | topN=%d | num_vars=%s",
                 counts_csv, meta_csv, ifelse(nzchar(mirna_list_fn), mirna_list_fn, "auto"),
                 ifelse(nzchar(outcome_col), outcome_col, "auto"),
                 ifelse(nzchar(positive_label), positive_label, "auto"),
                 top_n_fallback,
                 ifelse(nzchar(numeric_vars_fn), numeric_vars_fn, "auto")))

# --- Librerías necesarias ---
suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pROC)
}))

# --- Lectura de datos ---
if (!file.exists(counts_csv)) { log_line(sprintf("ERROR: counts file not found: %s", counts_csv)); quit(status = 1) }
if (!file.exists(meta_csv))   { log_line(sprintf("ERROR: metadata file not found: %s", meta_csv)); quit(status = 1) }

expr <- tryCatch(read.csv(counts_csv, check.names = FALSE, row.names = 1), error = function(e) NULL)
meta <- tryCatch(read.csv(meta_csv, stringsAsFactors = FALSE), error = function(e) NULL)
if (is.null(expr) || is.null(meta)) { log_line("ERROR: failed to read inputs."); quit(status = 1) }

if (!("sample_id" %in% names(meta))) { log_line("ERROR: metadata must include 'sample_id'."); quit(status = 1) }

# Alinear muestras
common <- intersect(colnames(expr), meta$sample_id)
if (length(common) < 4) {
  log_line(sprintf("ERROR: too few overlapping samples (<4). expr(%d) ∩ meta(%d) = %d",
                   ncol(expr), nrow(meta), length(common)))
  quit(status = 1)
}
expr <- expr[, common, drop = FALSE]
meta <- meta[match(common, meta$sample_id), , drop = FALSE]
stopifnot(identical(colnames(expr), meta$sample_id))

# Transformación para correlaciones/clasificación: log2(1+normalized)
expr_log <- log2(1 + as.matrix(expr))

# --- Selección de miRNAs a evaluar ---
read_list <- function(path) {
  if (!file.exists(path) || !nzchar(path)) return(character(0))
  ext <- tolower(file_ext(path))
  if (ext %in% c("csv", "tsv")) {
    df <- tryCatch({
      if (ext == "csv") readr::read_csv(path, show_col_types = FALSE) else readr::read_tsv(path, show_col_types = FALSE)
    }, error = function(e) NULL)
    if (is.null(df)) return(character(0))
    cols <- tolower(names(df))
    for (cand in c("mirna","id","name","feature")) {
      if (cand %in% cols) return(unique(df[[names(df)[match(cand, cols)]]]))
    }
    # si una sola columna
    if (ncol(df) == 1) return(unique(df[[1]]))
    return(character(0))
  } else {
    v <- readLines(path, warn = FALSE)
    v <- trimws(v); v <- v[nchar(v) > 0]
    unique(v)
  }
}

miR_list <- character(0)
if (nzchar(mirna_list_fn)) {
  miR_list <- read_list(mirna_list_fn)
}

if (!length(miR_list)) {
  # Intentar desde DEA: tomar todos los DEA_DESeq2_* sig y elegir top_n_fallback por frecuencia/padj
  dea_dir <- file.path(ROOT_DIR, "results", "DEA")
  dea_files <- list.files(dea_dir, pattern = "^DEA_DESeq2_.*_vs_.*\\.csv$", full.names = TRUE)
  sig_all <- list()
  for (f in dea_files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    if (!all(c("miRNA","padj","log2FoldChange") %in% names(df))) next
    df <- df %>% filter(!is.na(padj), padj <= 0.05, !is.na(log2FoldChange), abs(log2FoldChange) >= 1)
    if (nrow(df)) sig_all[[length(sig_all)+1]] <- df
  }
  if (length(sig_all)) {
    sig <- bind_rows(sig_all)
    miR_list <- sig %>%
      group_by(miRNA) %>%
      summarize(min_padj = min(padj, na.rm = TRUE),
                max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
                .groups = "drop") %>%
      arrange(min_padj, desc(max_abs_lfc)) %>%
      slice_head(n = top_n_fallback) %>%
      pull(miRNA)
    log_line(sprintf("Selected miRNAs from DEA (n=%d).", length(miR_list)))
  }
}

if (!length(miR_list)) {
  # Fallback: top variable miRNAs en expr_log
  vrs <- apply(expr_log, 1, stats::var, na.rm = TRUE)
  miR_list <- names(sort(vrs, decreasing = TRUE))[seq_len(min(top_n_fallback, length(vrs)))]
  log_line(sprintf("Selected top variable miRNAs (n=%d).", length(miR_list)))
}

# Filtrar expr a los miRNAs de interés (si existen en la matriz)
miR_list <- intersect(miR_list, rownames(expr_log))
if (!length(miR_list)) { log_line("ERROR: None of the selected miRNAs are present in expression matrix."); quit(status = 1) }

# --- Determinar variables clínicas numéricas ---
read_numeric_vars <- function(meta, provided_file = "") {
  # prioridad por nombres comunes; luego todas numéricas
  preferred <- c("Mayo", "mayo", "CDAI", "cdai", "CRP", "crp",
                 "calprotectina", "Calprotectina", "calprotectin", "FCP", "fcp")
  if (nzchar(provided_file) && file.exists(provided_file)) {
    v <- read_list(provided_file)
    v <- v[v %in% names(meta)]
    # coaccionar a numeric si es char
    for (nm in v) {
      if (!is.numeric(meta[[nm]])) suppressWarnings(meta[[nm]] <- as.numeric(meta[[nm]]))
    }
    v <- v[sapply(meta[v], function(x) is.numeric(x) && sum(!is.na(x)) >= 4)]
    return(v)
  }
  # Detectar
  cand <- intersect(preferred, names(meta))
  # coercionar
  for (nm in cand) if (!is.numeric(meta[[nm]])) suppressWarnings(meta[[nm]] <- as.numeric(meta[[nm]]))
  cand <- cand[sapply(meta[cand], function(x) is.numeric(x) && sum(!is.na(x)) >= 4)]
  # si no hay, tomar todas numéricas razonables
  if (!length(cand)) {
    nums <- names(Filter(function(x) is.numeric(x) && sum(!is.na(x)) >= 4, meta))
    cand <- setdiff(nums, c("age_if_not_relevant"))  # placeholder para excluir si quieres
  }
  cand
}

num_vars <- read_numeric_vars(meta, numeric_vars_fn)
if (!length(num_vars)) {
  log_line("WARN: No numeric clinical variables detected; correlations will be skipped.")
}

# --- Outcome binario ---
infer_outcome <- function(meta, desired_col = "", pos_label = "") {
  if (nzchar(desired_col) && desired_col %in% names(meta)) {
    y <- meta[[desired_col]]
    # inferir positivo si no viene
    if (!nzchar(pos_label)) {
      # Heurística
      if (grepl("response|respon", desired_col, ignore.case = TRUE)) {
        pos_label <- "Responder"
      } else {
        pos_label <- "Case"
      }
    }
    return(list(col = desired_col, pos = pos_label))
  }
  # prefer 'response'/'responder'
  for (nm in names(meta)) {
    if (grepl("^response$|respon", nm, ignore.case = TRUE)) return(list(col = nm, pos = ifelse(nzchar(pos_label), pos_label, "Responder")))
  }
  # fallback: disease
  if ("disease" %in% names(meta)) return(list(col = "disease", pos = ifelse(nzchar(pos_label), pos_label, "Case")))
  # not found
  return(list(col = "", pos = ""))
}
oc <- infer_outcome(meta, outcome_col, positive_label)
if (!nzchar(oc$col)) {
  log_line("WARN: No suitable outcome column found; ROC/AUC will be skipped.")
}

# Construir vector binario si posible
make_binary_outcome <- function(vec, positive_label = "Case") {
  v <- as.character(vec)
  # Detectar 'Control'
  neg_label <- "Control"
  # Si positive_label no aparece, asumir todo lo no-Control es positivo
  if (!any(v == positive_label, na.rm = TRUE)) {
    pos_idx <- which(!is.na(v) & !grepl("^control$", v, ignore.case = TRUE))
    y <- rep(NA_integer_, length(v))
    y[!is.na(v) & grepl("^control$", v, ignore.case = TRUE)] <- 0L
    y[pos_idx] <- 1L
    return(factor(y, levels = c(0,1), labels = c("Control","Case")))
  } else {
    y <- ifelse(v == positive_label, 1L,
                ifelse(grepl("^control$", v, ignore.case = TRUE), 0L, NA_integer_))
    # Si aún quedan NA, conviértelos a 1 (positivo) si no son Control
    y[is.na(y) & !grepl("^control$", v, ignore.case = TRUE) & !is.na(v)] <- 1L
    pos_name <- ifelse(positive_label == "Case", "Case", positive_label)
    return(factor(y, levels = c(0,1), labels = c("Control", pos_name)))
  }
}

y_bin <- NULL
if (nzchar(oc$col)) {
  y_bin <- make_binary_outcome(meta[[oc$col]], oc$pos)
  if (sum(!is.na(y_bin)) < 4 || length(unique(na.omit(y_bin))) < 2) {
    log_line("WARN: Outcome has insufficient data/variation; ROC/AUC will be skipped.")
    y_bin <- NULL
  } else {
    log_line(sprintf("Outcome: %s (positive='%s') | n=%d (non-NA)", oc$col, oc$pos, sum(!is.na(y_bin))))
  }
}

# --- (1) Correlaciones Spearman ---
corr_out <- data.frame()
if (length(num_vars)) {
  log_line(sprintf("Running Spearman correlations for %d miRNAs x %d clinical vars ...", length(miR_list), length(num_vars)))
  for (mir in miR_list) {
    x <- as.numeric(expr_log[mir, ])
    for (v in num_vars) {
      y <- meta[[v]]
      # coercionar a numeric si hace falta
      if (!is.numeric(y)) suppressWarnings(y <- as.numeric(y))
      if (sum(!is.na(x) & !is.na(y)) >= 4) {
        ct <- suppressWarnings(suppressMessages(cor.test(x, y, method = "spearman", exact = FALSE)))
        rho <- unname(ct$estimate)
        p   <- unname(ct$p.value)
        corr_out <- rbind(corr_out, data.frame(miRNA = mir, variable = v, rho = rho, p = p, n = sum(!is.na(x) & !is.na(y)), stringsAsFactors = FALSE))
      }
    }
  }
  if (nrow(corr_out)) {
    # Ajuste FDR por variable (o global)
    corr_out$padj_global <- p.adjust(corr_out$p, method = "BH")
    # por variable (útil si interpretas por clínico)
    corr_out <- corr_out %>%
      group_by(variable) %>%
      mutate(padj_by_var = p.adjust(p, method = "BH")) %>%
      ungroup() %>%
      arrange(padj_global, desc(abs(rho)))
    readr::write_csv(corr_out, file.path(CLIN_DIR, "correlations.csv"))
    log_line(sprintf("Saved: correlations.csv (rows=%d)", nrow(corr_out)))
  } else {
    log_line("No valid pairs for correlations.")
  }
}

# --- (2) ROC/AUC por miRNA ---
roc_tbl <- data.frame()
plot_roc <- function(roc_obj, title, outfile) {
  p <- ggplot() +
    geom_line(aes(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal() +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity")
  ggsave(outfile, p, width = 5, height = 5, dpi = 300)
}

if (!is.null(y_bin)) {
  log_line("Computing ROC/AUC per miRNA ...")
  y_num <- as.numeric(y_bin) - 1L  # 0/1
  for (mir in miR_list) {
    x <- as.numeric(expr_log[mir, ])
    keep <- is.finite(x) & !is.na(y_num)
    if (length(unique(y_num[keep])) < 2 || sum(keep) < 4) next
    roc_obj <- tryCatch(pROC::roc(response = y_bin[keep], predictor = x[keep], levels = levels(y_bin), direction = "<"), error = function(e) NULL)
    if (is.null(roc_obj)) next
    auc_val <- as.numeric(pROC::auc(roc_obj))
    roc_tbl <- rbind(roc_tbl, data.frame(miRNA = mir, AUC = auc_val, n = sum(keep), stringsAsFactors = FALSE))
  }
  if (nrow(roc_tbl)) {
    roc_tbl <- roc_tbl[order(-roc_tbl$AUC), , drop = FALSE]
    readr::write_csv(roc_tbl, file.path(CLIN_DIR, "roc_auc.csv"))
    log_line(sprintf("Saved: roc_auc.csv (rows=%d)", nrow(roc_tbl)))
    # Graficar top 6 ROC
    top_plot <- head(roc_tbl$miRNA, 6)
    for (mir in top_plot) {
      x <- as.numeric(expr_log[mir, ])
      keep <- is.finite(x) & !is.na(y_bin)
      roc_obj <- pROC::roc(response = y_bin[keep], predictor = x[keep], levels = levels(y_bin), direction = "<")
      outpng <- file.path(FIG_DIR, sprintf("ROC_%s.png", gsub("[^A-Za-z0-9._-]", "_", mir)))
      plot_roc(roc_obj, sprintf("ROC: %s (AUC=%.3f)", mir, as.numeric(pROC::auc(roc_obj))), outpng)
    }
  } else {
    log_line("No ROC/AUC could be computed.")
  }
}

# --- (3) Figuras de expresión por grupos (opcional, top 12 miRNAs) ---
if ("disease" %in% names(meta)) {
  log_line("Plotting expression by disease groups (top 12 miRNAs) ...")
  to_plot <- head(miR_list, 12)
  df_long <- data.frame(
    miRNA = rep(to_plot, each = ncol(expr_log)),
    expr  = as.numeric(unlist(lapply(to_plot, function(m) expr_log[m, ]))),
    sample_id = rep(colnames(expr_log), times = length(to_plot)),
    disease = rep(meta$disease, times = length(to_plot)),
    stringsAsFactors = FALSE
  )
  for (mir in to_plot) {
    d <- df_long[df_long$miRNA == mir, , drop = FALSE]
    p <- ggplot(d, aes(x = disease, y = expr)) +
      geom_violin(trim = TRUE) + geom_boxplot(width = 0.15, outlier.size = 0.8) +
      labs(title = sprintf("%s expression by group", mir), x = "Group", y = "log2(1+normalized)") +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
    ggsave(file.path(FIG_DIR, sprintf("expr_by_group_%s.png", gsub("[^A-Za-z0-9._-]", "_", mir))),
           p, width = 6.5, height = 4.5, dpi = 300)
  }
}

# --- (4) Modelo simple multi-miRNA con validación cruzada (5-fold) ---
if (!is.null(y_bin) && nrow(roc_tbl) >= 3) {
  log_line("Building simple multi-miRNA logistic model with 5-fold CV (top 5 miRNAs by AUC) ...")
  top_mir <- head(roc_tbl$miRNA, 5)
  X <- t(expr_log[top_mir, , drop = FALSE])
  X <- scale(X)  # estandarizar
  dat <- data.frame(y = as.numeric(y_bin) - 1L, X, check.names = FALSE)
  set.seed(123)
  k <- 5
  n <- nrow(dat)
  folds <- sample(rep(1:k, length.out = n))
  probs <- rep(NA_real_, n)
  for (fold in 1:k) {
    train_idx <- which(folds != fold)
    test_idx  <- which(folds == fold)
    fit <- tryCatch(glm(y ~ ., data = dat[train_idx, , drop = FALSE], family = binomial()), error = function(e) NULL)
    if (is.null(fit)) next
    pr <- suppressWarnings(predict(fit, newdata = dat[test_idx, , drop = FALSE], type = "response"))
    probs[test_idx] <- pr
  }
  keep <- !is.na(probs)
  if (sum(keep) >= 4 && length(unique(dat$y[keep])) == 2) {
    roc_obj <- pROC::roc(response = factor(dat$y[keep], levels = c(0,1), labels = levels(y_bin)), predictor = probs[keep], direction = "<")
    auc_cv <- as.numeric(pROC::auc(roc_obj))
    readr::write_csv(
      data.frame(model = "logit_top5", AUC_CV = auc_cv, n = sum(keep), miRNAs = paste(top_mir, collapse = ",")),
      file.path(CLIN_DIR, "multimirna_cv_auc.csv")
    )
    log_line(sprintf("Saved: multimirna_cv_auc.csv (AUC_CV=%.3f, n=%d)", auc_cv, sum(keep)))
    # Figura ROC del modelo
    p <- ggplot() +
      geom_line(aes(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      coord_equal() +
      labs(title = sprintf("Multi-miRNA logistic (AUC=%.3f)", auc_cv),
           x = "1 - Specificity", y = "Sensitivity")
    ggsave(file.path(FIG_DIR, "ROC_multimirna_top5.png"), p, width = 5, height = 5, dpi = 300)
  } else {
    log_line("WARN: Could not compute CV ROC for multi-miRNA model (insufficient variation).")
  }
}

# --- Fin ---
log_line("Script finalizado correctamente.")
