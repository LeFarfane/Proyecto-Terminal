#!/usr/bin/env Rscript

# ---------------------------
# 0) Setup & Interactive Prompt
# ---------------------------
cat("\n=====================================================\n")
cat("      DESeq2 Automated miRNA Pipeline\n")
cat("=====================================================\n")

find_file <- function(name) {
  hits <- list.files(".", pattern = paste0("^", name, "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop(paste("Cannot find", name, "— run from the dataset root directory."))
  if (length(hits) > 1) cat(paste("[WARN] Multiple matches for", name, "— using first:", hits[1], "\n"))
  hits[1]
}
counts_csv  <- find_file("Master_Filtered_Counts_DESeq2.csv")
meta_csv    <- find_file("sample_metadata.csv")

# Read metadata early so we can show real column options
meta_preview <- read.csv(meta_csv, stringsAsFactors = FALSE, nrows = 0)
candidate_cols <- colnames(meta_preview)[colnames(meta_preview) != "sample_id"]

cat("Available metadata columns:\n")
for (i in seq_along(candidate_cols)) {
  cat(sprintf("  [%d] %s\n", i, candidate_cols[i]))
}
cat("\nWhich column would you like to analyze? (type name or number): ")

input_col <- trimws(readLines(file("stdin"), n = 1))

# Accept number or name
if (grepl("^[0-9]+$", input_col)) {
  idx <- as.integer(input_col)
  if (idx >= 1 && idx <= length(candidate_cols)) {
    condition_col <- candidate_cols[idx]
  } else {
    stop(paste("Invalid selection:", input_col))
  }
} else if (input_col == "") {
  condition_col <- candidate_cols[1]
  cat(paste(" -> No input detected. Defaulting to:", condition_col, "\n"))
} else {
  condition_col <- input_col
}
cat(paste(" -> Selected:", condition_col, "\n"))

results_dir_strict  <- "./DEA_results/DEA_No_Outlier_Replacement"
results_dir_default <- "./DEA_results/DEA_Outlier_Replacement_On"
fig_dir             <- "./DEA_results/figures_and_QC"

for (d in c(results_dir_strict, results_dir_default, fig_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---------------------------
# 1) Load data & Prep
# ---------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, stringsAsFactors = FALSE)

if (!(condition_col %in% colnames(meta))) {
  stop(paste("Column '", condition_col, "' not found. Available:", paste(colnames(meta), collapse = ", ")))
}

# --- STANDARDIZE & SHORTEN METADATA LABELS ---
raw_conditions <- meta[[condition_col]]

clean_conditions <- sapply(raw_conditions, function(x) {
  if (grepl("^NC$|healthy.control|normal.control", x, ignore.case = TRUE)) {
    return("HealthyControl")
  } else if (grepl("^UC$|ulcerative.colitis", x, ignore.case = TRUE)) {
    return("UlcerativeColitis")
  } else if (grepl("non-responder", x, ignore.case = TRUE)) {
    return("UC_NonResponder")
  } else if (grepl("good.response|treatment.responder", x, ignore.case = TRUE)) {
    return("UC_Responder")
  } else if (grepl("^IlealCD$|ileal.crohn|crohn", x, ignore.case = TRUE)) {
    return("IlealCD")
  } else {
    return(make.names(substr(x, 1, 20)))
  }
})

meta$Condition <- clean_conditions
control_label  <- "HealthyControl"

# Align columns and rows
a <- colnames(counts)
b <- meta$sample_id
common <- intersect(a, b)
counts <- counts[, common, drop = FALSE]
meta   <- meta[match(common, meta$sample_id), , drop = FALSE]

rownames(meta) <- meta$sample_id

# Set factors and reference levels
meta$Condition <- as.factor(meta$Condition)
if (!(control_label %in% levels(meta$Condition))) {
  stop(paste("Control label", control_label, "not found after text cleaning!"))
}
meta$Condition <- relevel(meta$Condition, ref = control_label)
contrast_levels <- setdiff(levels(meta$Condition), control_label)

covars <- c("batch", "age", "sex", "tissue_or_biofluid", "platform")
avail_covars <- covars[covars %in% colnames(meta) & sapply(covars, function(x) dplyr::n_distinct(meta[[x]]) > 1)]
design_formula <- as.formula(paste("~", paste(c(avail_covars, "Condition"), collapse = " + ")))

cat(paste("\n[INFO] Design formula:", deparse(design_formula), "\n"))

# ---------------------------
# 2) Base DESeq2 Object & QC Plots
# ---------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(umap)
  library(pheatmap)
})

dds_base <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                                   colData   = meta,
                                   design    = design_formula)

cat("[INFO] Calculating Variance Stabilizing Transformation...\n")
vst_obj <- varianceStabilizingTransformation(dds_base, blind = TRUE, fitType = "local")
vst_mat <- assay(vst_obj)
write.csv(vst_mat, file.path(fig_dir, "vst_normalized_matrix.csv"))

anno_df <- data.frame(Condition = meta$Condition)
rownames(anno_df) <- rownames(meta)

pheatmap::pheatmap(vst_mat, show_rownames = TRUE, show_colnames = FALSE,
                   fontsize_row = 5, 
                   annotation_col = anno_df,
                   filename = file.path(fig_dir, "heatmap_ALL_miRNAs.png"), width = 8, height = 10)

pca_data <- plotPCA(vst_obj, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() + labs(title = paste("PCA Plot colored by", condition_col))
ggsave(file.path(fig_dir, "PCA_plot.png"), p_pca, width = 7, height = 5)

set.seed(42)
umap_res <- umap(t(vst_mat))
umap_df <- data.frame(UMAP1 = umap_res$layout[,1], UMAP2 = umap_res$layout[,2], Condition = meta$Condition)
p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() + labs(title = paste("UMAP colored by", condition_col))
ggsave(file.path(fig_dir, "UMAP_plot.png"), p_umap, width = 7, height = 5)


# ---------------------------
# 3) The Dual DESeq2 Runs & S_miRNA Score
# ---------------------------
suppressPackageStartupMessages({ library(apeglm) })

run_dea_pipeline <- function(dds, out_dir, run_name) {
  
  png(file.path(out_dir, paste0("dispersion_plot_", run_name, ".png")), width = 800, height = 600)
  plotDispEsts(dds, main = paste("Dispersion Estimates -", run_name))
  dev.off()
  
  res_names <- resultsNames(dds)
  
  for (grp in contrast_levels) {
    coef_name <- paste0("Condition_", grp, "_vs_", control_label)
    if (!(coef_name %in% res_names)) next
    
    # 1. Grab original results to keep the 'stat' column
    res_orig <- results(dds, name = coef_name)
    
    # 2. Grab shrunken results for accurate log2FoldChange
    res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    
    # 3. Stitch them together and calculate the S_miRNA score
    res_tbl <- as.data.frame(res_orig) %>%
      rownames_to_column(var = "miRNA") %>%
      mutate(
        # Replace original LFC and lfcSE with the shrunken values
        log2FoldChange = res_shr$log2FoldChange,
        lfcSE = res_shr$lfcSE,
        
        # Math Safety: Prevent -log10(0) or -log10(NA) crashes
        padj_safe = ifelse(is.na(padj) | padj == 0, 1e-300, padj),
        
        # FIXED: Removed the abs() so the score keeps its direction!
        S_miRNA = log2FoldChange * -log10(padj_safe)
      ) %>%
      # Select exactly the columns you requested
      select(miRNA, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, S_miRNA) %>%
      # FIXED: Sort by the absolute magnitude so the most significant 
      # up AND down-regulated miRNAs stay at the top of the CSV
      arrange(desc(abs(S_miRNA)))
    
    # Export CSV
    write.csv(res_tbl, file.path(out_dir, paste0("DEA_", run_name, "_", grp, ".csv")), row.names = FALSE)
    
    # Volcano Plot
    df <- res_tbl %>% mutate(sig = ifelse(!is.na(padj) & padj <= 0.05 & abs(log2FoldChange) >= 1, "sig", "ns"))
    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(alpha = sig, color = sig), size = 1.5) +
      scale_color_manual(values = c("ns" = "grey", "sig" = "red")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_minimal() + labs(title = paste(run_name, "-", grp, "vs Control"))
    ggsave(file.path(out_dir, paste0("volcano_", run_name, "_", grp, ".png")), p, width = 7, height = 5)
  }
}

cat("\nRunning Strict Model (minReplicatesForReplace = Inf)...\n")
dds_strict <- DESeq(dds_base, minReplicatesForReplace = Inf, fitType = "local", parallel = FALSE)
run_dea_pipeline(dds_strict, results_dir_strict, "Strict")

cat("\nRunning Default Model (Cook's distance replacement enabled)...\n")
dds_default <- DESeq(dds_base, fitType = "local", parallel = FALSE)
run_dea_pipeline(dds_default, results_dir_default, "Default")

cat("\n[SUCCESS] Both runs complete! Check the 'results' folder.\n")