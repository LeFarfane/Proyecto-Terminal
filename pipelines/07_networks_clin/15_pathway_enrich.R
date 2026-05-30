#!/usr/bin/env Rscript
# ===========================================
# Script: 15_pathway_enrich.R
# Description: Performs GO (BP, CC, MF) and KEGG pathway enrichment 
#              on a list of target genes using clusterProfiler.
#              Features a Smart File Finder for CLI path inputs.
# ===========================================

suppressWarnings(suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
}))

# --- Setup Args & Paths ---
RUN_DIR <- getwd()
args <- commandArgs(trailingOnly = TRUE)
# Argument 1: Input directory or exact file path. Defaults to current directory.
input_arg <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else RUN_DIR

# --- Smart File Finder ---
if (dir.exists(input_arg)) {
  # If given a directory, scan for the file in common subfolders
  potential_paths <- c(
    file.path(input_arg, "validated_target_genes_unique.csv"),
    file.path(input_arg, "multimir_outputs", "validated_target_genes_unique.csv"),
    file.path(input_arg, "multimir_output", "validated_target_genes_unique.csv") # Handles the missing 's' typo!
  )
  input_file <- Find(file.exists, potential_paths)
  
  if (is.null(input_file)) {
    stop(sprintf("Cannot find 'validated_target_genes_unique.csv' inside %s or its subfolders.", input_arg))
  }
} else if (file.exists(input_arg)) {
  # If given the exact file directly
  input_file <- input_arg
} else {
  stop(sprintf("Invalid input path: %s", input_arg))
}

# --- Output Setup ---
OUT_DIR <- file.path(RUN_DIR, "pathway_outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LOGFILE <- file.path(OUT_DIR, "pathway_enrichment.log")
sink(LOGFILE, append = FALSE, split = TRUE)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- function(...) cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = "")))

log_line("Starting 15_pathway_enrich.R")
log_line(sprintf("Using input file: %s", input_file))

# --- 1. Load Gene List ---
genes_df <- read_csv(input_file, show_col_types = FALSE)
gene_symbols <- unique(genes_df$target_symbol)
log_line(sprintf("Loaded %d unique target genes for enrichment.", length(gene_symbols)))

# --- 2. Gene ID Conversion ---
log_line("Translating Gene Symbols to Entrez IDs...")
gene_conversion <- tryCatch({
  bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) {
  log_line("WARN: Issue during ID conversion: ", e$message)
  return(NULL)
})

if (is.null(gene_conversion) || nrow(gene_conversion) == 0) {
   stop("FATAL: Failed to translate any gene symbols to Entrez IDs.")
}

entrez_genes <- gene_conversion$ENTREZID
log_line(sprintf("Successfully translated %d genes to Entrez IDs.", length(entrez_genes)))

# --- 3. GO Enrichment (Gene Ontology) ---
log_line("Running GO Enrichment (BP, CC, MF)...")
ego <- enrichGO(gene          = entrez_genes,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",  
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)   

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  write_csv(as.data.frame(ego), file.path(OUT_DIR, "GO_Enrichment_Results.csv"))
  log_line("Saved GO_Enrichment_Results.csv")
  
  p_go <- dotplot(ego, split="ONTOLOGY", showCategory=10) + 
          facet_grid(ONTOLOGY~., scale="free") +
          ggtitle("Top 10 GO Terms per Ontology")
  ggsave(file.path(OUT_DIR, "GO_Dotplot.png"), p_go, width = 10, height = 12)
} else {
  log_line("WARN: No significant GO terms found.")
}

# --- 4. KEGG Pathway Enrichment ---
log_line("Running KEGG Pathway Enrichment...")
ekegg <- enrichKEGG(gene         = entrez_genes,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write_csv(as.data.frame(ekegg_readable), file.path(OUT_DIR, "KEGG_Enrichment_Results.csv"))
  log_line("Saved KEGG_Enrichment_Results.csv")
  
  p_kegg <- dotplot(ekegg, showCategory=20) + 
            ggtitle("Top 20 Enriched KEGG Pathways")
  ggsave(file.path(OUT_DIR, "KEGG_Dotplot.png"), p_kegg, width = 9, height = 7)
} else {
  log_line("WARN: No significant KEGG pathways found.")
}

log_line("Pathway enrichment pipeline finished successfully!")
sink()