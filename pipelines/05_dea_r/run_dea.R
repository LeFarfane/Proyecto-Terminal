#!/usr/bin/env Rscript
# ===========================================
# Script: 05_dea_r/run_dea.R
# Descripción: Differential Expression Analysis (DEA) for miRNA counts using DESeq2 (primary)
#              and optional edgeR/limma-voom. Aligns counts↔metadata, filters low abundance,
#              builds a dynamic design (~ covariates + disease), exports tables and figures.
# Fecha: 2025-09-07
# ===========================================

# --- Configuración segura / Logging ---
suppressWarnings(suppressMessages({
  library(utils)
}))

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- function(...) {
  msg <- paste0("[", timestamp(), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  if (!is.null(getOption("dea_logfile"))) {
    cat(msg, "\n", file = getOption("dea_logfile"), append = TRUE)
  }
}

# --- Paths de proyecto (relativos a PT/) ---
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg)))
  # Fallback (interactive or sourced)
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile))
  return(normalizePath(getwd()))
}
SCRIPT_PATH <- tryCatch(get_script_path(), error = function(e) normalizePath(getwd()))
ROOT_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), "..", ".."))  # .../PT

# --- Defaults (sobre-escribibles por args) ---
default_counts <- file.path(ROOT_DIR, "data_processed", "miRNA_counts_raw.csv")
default_meta   <- file.path(ROOT_DIR, "metadata", "sample_metadata.csv")
DEA_DIR   <- file.path(ROOT_DIR, "results", "DEA")
FIG_DIR   <- file.path(DEA_DIR, "figures")
LOG_DIR   <- file.path(ROOT_DIR, "outputs", "logs")
dir.create(DEA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOG_DIR, "05_dea_r.log")
options(dea_logfile = LOGFILE)

# --- CLI Args (muy simples) ---
args <- commandArgs(trailingOnly = TRUE)
# Posicionales: counts_csv, meta_csv, control_label, padj_cut, lfc_cut
counts_csv   <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_counts
meta_csv     <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else default_meta
control_lab  <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else "Control"
padj_cut     <- if (length(args) >= 4 && nzchar(args[[4]])) as.numeric(args[[4]]) else 0.05
lfc_cut      <- if (length(args) >= 5 && nzchar(args[[5]])) as.numeric(args[[5]]) else 1.0
run_voom     <- if (length(args) >= 6 && nzchar(args[[6]])) as.logical(args[[6]]) else FALSE

# --- Inicio ---
write("", file = LOGFILE)  # reset log
log_line("Iniciando script...")
log_line(sprintf("ROOT_DIR: %s", ROOT_DIR))
log_line(sprintf("counts_csv: %s", counts_csv))
log_line(sprintf("meta_csv  : %s", meta_csv))
log_line(sprintf("control   : %s | padj_cut: %.3f | lfc_cut: %.2f | voom: %s", control_lab, padj_cut, lfc_cut, run_voom))

# --- Librerías necesarias ---
suppressWarnings(suppressMessages({
  library(readr)
  library(tibble)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(apeglm)
  library(pheatmap)
  library(matrixStats)
}))

# edgeR/limma opcional
if (isTRUE(run_voom)) {
  suppressWarnings(suppressMessages({
    library(edgeR)
    library(limma)
  }))
}

# --- Lectura de datos ---
if (!file.exists(counts_csv)) {
  log_line(sprintf("ERROR: counts file not found: %s", counts_csv)); quit(status = 1)
}
if (!file.exists(meta_csv)) {
  log_line(sprintf("ERROR: metadata file not found: %s", meta_csv)); quit(status = 1)
}

counts <- tryCatch({
  read.csv(counts_csv, check.names = FALSE, row.names = 1)
}, error = function(e) { log_line(sprintf("ERROR reading counts: %s", e$message)); quit(status = 1) })

meta <- tryCatch({
  read.csv(meta_csv, stringsAsFactors = FALSE)
}, error = function(e) { log_line(sprintf("ERROR reading metadata: %s", e$message)); quit(status = 1) })

# --- Alineación de muestras ---
if (!("sample_id" %in% colnames(meta))) {
  log_line("ERROR: metadata must include 'sample_id' column."); quit(status = 1)
}
if (!("disease" %in% colnames(meta))) {
  log_line("ERROR: metadata must include 'disease' column."); quit(status = 1)
}

a <- colnames(counts)
b <- meta$sample_id
common <- intersect(a, b)
if (length(common) < 4) {
  log_line(sprintf("ERROR: too few overlapping samples (<4). counts(%d) ∩ meta(%d) = %d",
                   length(a), length(b), length(common)))
  quit(status = 1)
}
counts <- counts[, common, drop = FALSE]
meta   <- meta[match(common, meta$sample_id), , drop = FALSE]
stopifnot(identical(colnames(counts), meta$sample_id))

# --- Sanitizar tipos y covariables ---
meta$disease <- as.factor(meta$disease)
covars <- c("batch", "age", "sex", "tissue_or_biofluid", "platform")
avail_covars <- covars[covars %in% colnames(meta) & sapply(covars, function(x) dplyr::n_distinct(meta[[x]]) > 1)]
# convertir covars categóricas a factor (age puede ser numérica)
for (cv in avail_covars) {
  if (cv != "age") {
    if (!is.factor(meta[[cv]])) meta[[cv]] <- as.factor(meta[[cv]])
  } else {
    # age: si es char, intenta numeric
    if (!is.numeric(meta[[cv]])) {
      suppressWarnings(meta[[cv]] <- as.numeric(meta[[cv]]))
    }
  }
}

# --- Control como referencia (autodetección si no está) ---
if (!(control_lab %in% levels(meta$disease))) {
  candidates <- levels(meta$disease)[grepl("control|healthy", levels(meta$disease), ignore.case = TRUE)]
  if (length(candidates) >= 1) {
    log_line(sprintf("Aviso: control '%s' no encontrado. Usando detectado: '%s'.", control_lab, candidates[1]))
    control_lab <- candidates[1]
  } else {
    log_line(sprintf("ERROR: control label '%s' not in disease levels: %s",
                     control_lab, paste(levels(meta$disease), collapse = ", ")))
    quit(status = 1)
  }
}
meta$disease <- relevel(meta$disease, ref = control_lab)

# --- Filtrado básico de baja abundancia ---
keep <- rowSums(counts >= 10) >= max(2, floor(0.2 * ncol(counts)))
counts_f <- counts[keep, , drop = FALSE]
log_line(sprintf("Filtrado: %d -> %d miRNAs (>=10 counts in >=%d samples).",
                 nrow(counts), nrow(counts_f), max(2, floor(0.2 * ncol(counts)))))

# --- Fórmula de diseño dinámica ---
form_terms <- c(avail_covars, "disease")
design_formula <- as.formula(paste("~", paste(form_terms, collapse = " + ")))
log_line(sprintf("Design formula: %s", paste(deparse(design_formula), collapse = "")))

# --- DESeq2 pipeline ---
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts_f)),
                              colData   = meta,
                              design    = design_formula)

dds <- DESeq(dds, minReplicatesForReplace = Inf, parallel = FALSE)
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(DEA_DIR, "normalized_counts_deseq2.csv"))
log_line("Saved: normalized_counts_deseq2.csv")

# VST para visualización
vst_mat <- assay(vst(dds, blind = TRUE))
write.csv(vst_mat, file.path(DEA_DIR, "vst_matrix.csv"))
log_line("Saved: vst_matrix.csv")

# Contrastes (cada nivel de disease contra control)
contrast_levels <- setdiff(levels(meta$disease), control_lab)

# Función de volcano
plot_volcano <- function(res_tbl, title, outfile, lfc_cut = 1, padj_cut = 0.05) {
  df <- res_tbl
  df$neglog10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
  df$sig <- ifelse(!is.na(df$padj) & df$padj <= padj_cut & abs(df$log2FoldChange) >= lfc_cut, "sig", "ns")
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(alpha = sig), size = 1.3) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    labs(title = title, x = "log2 fold change", y = "-log10(FDR)") +
    theme_minimal()
  ggsave(outfile, p, width = 6, height = 4, dpi = 300)
}

# Resultados por contraste con shrinkage (apeglm)
suppressWarnings(suppressMessages(library(apeglm)))
for (grp in contrast_levels) {
  log_line(sprintf("DESeq2: %s vs %s", grp, control_lab))
  res <- results(dds, contrast = c("disease", grp, control_lab))
  # shrink
  res_shr <- tryCatch(
    lfcShrink(dds, contrast = c("disease", grp, control_lab), type = "apeglm"),
    error = function(e) {
      log_line(sprintf("WARN: apeglm shrinkage failed for %s vs %s (%s). Using unshrunken LFC.",
                       grp, control_lab, e$message))
      res
    }
  )
  res_tbl <- as.data.frame(res_shr)
  res_tbl$miRNA <- rownames(res_tbl)
  res_tbl <- res_tbl %>% select(miRNA, log2FoldChange, lfcSE, stat, pvalue, padj) %>% arrange(padj)
  out_csv <- file.path(DEA_DIR, sprintf("DEA_DESeq2_%s_vs_%s.csv", grp, control_lab))
  write.csv(res_tbl, out_csv, row.names = FALSE)
  log_line(sprintf("Saved: %s (n=%d)", basename(out_csv), nrow(res_tbl)))

  # volcano
  volcano_png <- file.path(FIG_DIR, sprintf("volcano_DESeq2_%s_vs_%s.png", grp, control_lab))
  try(plot_volcano(res_tbl, sprintf("DESeq2: %s vs %s", grp, control_lab), volcano_png, lfc_cut, padj_cut), silent = TRUE)
  if (file.exists(volcano_png)) log_line(paste("Saved figure:", basename(volcano_png)))
}

# Heatmap top-variance (hasta 50)
suppressWarnings(suppressMessages(library(pheatmap)))
var_order <- order(matrixStats::rowVars(vst_mat), decreasing = TRUE)
sel <- head(var_order, n = min(50, nrow(vst_mat)))
ann <- meta %>% select(disease) %>% as.data.frame()
rownames(ann) <- meta$sample_id
heatmap_png <- file.path(FIG_DIR, "heatmap_top50_vst.png")
try(
  pheatmap::pheatmap(
    vst_mat[sel, , drop = FALSE],
    show_rownames = TRUE, show_colnames = FALSE,
    annotation_col = ann,
    filename = heatmap_png, width = 7, height = 9
  ),
  silent = TRUE
)
if (file.exists(heatmap_png)) log_line(paste("Saved figure:", basename(heatmap_png)))

# --- edgeR + limma-voom (opcional) ---
if (isTRUE(run_voom)) {
  log_line("Running edgeR + limma-voom (optional)...")
  y <- DGEList(counts = as.matrix(counts_f))
  y <- calcNormFactors(y, method = "TMM")
  # construir matriz de diseño
  form_str <- paste("~", paste(form_terms, collapse = " + "))
  meta_design <- meta
  for (cv in avail_covars) {
    if (cv != "age") {
      if (!is.factor(meta_design[[cv]])) meta_design[[cv]] <- as.factor(meta_design[[cv]])
    } else {
      if (!is.numeric(meta_design[[cv]])) { suppressWarnings(meta_design[[cv]] <- as.numeric(meta_design[[cv]])) }
    }
  }
  X <- model.matrix(as.formula(form_str), data = meta_design)
  v <- voom(y, design = X, plot = FALSE)
  fit <- lmFit(v, X)
  fit <- eBayes(fit)
  coef_names <- colnames(coef(fit))
  for (grp in contrast_levels) {
    coef_name <- paste0("disease", grp)
    if (!(coef_name %in% coef_names)) {
      log_line(sprintf("WARN: coeficiente %s no encontrado en voom; disponibles: %s", coef_name, paste(coef_names, collapse = ", ")))
      next
    }
    tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
    tt <- tt %>% tibble::rownames_to_column("miRNA") %>%
      dplyr::rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)
    out_csv <- file.path(DEA_DIR, sprintf("DEA_edgeR_voom_%s_vs_%s.csv", grp, control_lab))
    write.csv(tt, out_csv, row.names = FALSE)
    log_line(sprintf("Saved: %s (n=%d)", basename(out_csv), nrow(tt)))
    # volcano
    volcano_png <- file.path(FIG_DIR, sprintf("volcano_edgeR_voom_%s_vs_%s.png", grp, control_lab))
    try(plot_volcano(tt, sprintf("edgeR/voom: %s vs %s", grp, control_lab), volcano_png, lfc_cut, padj_cut), silent = TRUE)
    if (file.exists(volcano_png)) log_line(paste("Saved figure:", basename(volcano_png)))
  }
}

# --- Session info ---
sink(file.path(DEA_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

# --- Fin ---
log_line("Script finalizado correctamente.")
