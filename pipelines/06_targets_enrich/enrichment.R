#!/usr/bin/env Rscript
# ===========================================
# Script: 06_targets_enrich/enrichment.R
# Descripción: Enrichment analysis (KEGG, GO:BP, opcional Reactome) for validated target genes.
#              Reads target symbols (CSV/TXT), maps to ENTREZ, runs clusterProfiler,
#              and writes CSV tables + optional dotplot PNGs.
# Fecha: 2025-09-07
# ===========================================

suppressWarnings(suppressMessages({
  library(utils)
  library(tools)
}))

# --- Logging simple (con archivo) ---
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- function(...) {
  msg <- sprintf("[%s] %s", timestamp(), paste0(..., collapse = ""))
  cat(msg, "\n")
  if (!is.null(getOption("enrich_logfile"))) {
    cat(msg, "\n", file = getOption("enrich_logfile"), append = TRUE)
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

ENRICH_DIR <- file.path(ROOT_DIR, "results", "enrichment")
FIG_DIR    <- file.path(ENRICH_DIR, "figures")
LOG_DIR    <- file.path(ROOT_DIR, "outputs", "logs")
dir.create(ENRICH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,    recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOG_DIR, "06_enrichment.log")
options(enrich_logfile = LOGFILE)
write("", file = LOGFILE)

# --- Argumentos ---
# uso:
#   Rscript enrichment.R [targets_file] [organism_kegg] [p_cut] [q_cut] [minGS] [maxGS] [background_file]
# defaults:
#   targets_file: results/validated_target_genes_unique.csv (o fallback a results/targets_validated.csv)
#   organism_kegg: "hsa"
#   p_cut: 0.05, q_cut: 0.05, minGS: 10, maxGS: 500
args <- commandArgs(trailingOnly = TRUE)
default_targets <- {
  f1 <- file.path(ROOT_DIR, "results", "validated_target_genes_unique.csv")
  f2 <- file.path(ROOT_DIR, "results", "targets_validated.csv")
  if (file.exists(f1)) f1 else f2
}
targets_file <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_targets
org_kegg     <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else "hsa"
p_cut        <- if (length(args) >= 3 && nzchar(args[[3]])) as.numeric(args[[3]]) else 0.05
q_cut        <- if (length(args) >= 4 && nzchar(args[[4]])) as.numeric(args[[4]]) else 0.05
minGS        <- if (length(args) >= 5 && nzchar(args[[5]])) as.integer(args[[5]]) else 10
maxGS        <- if (length(args) >= 6 && nzchar(args[[6]])) as.integer(args[[6]]) else 500
bg_file      <- if (length(args) >= 7 && nzchar(args[[7]])) args[[7]] else ""

log_line("Starting enrichment.R")
log_line(sprintf("targets=%s | org_kegg=%s | p=%.3f | q=%.3f | minGS=%d | maxGS=%d | bg=%s",
                 targets_file, org_kegg, p_cut, q_cut, minGS, maxGS, ifelse(nzchar(bg_file), bg_file, "none")))

# --- Librerías de análisis ---
suppressWarnings(suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(enrichplot)
}))

have_reactome <- requireNamespace("ReactomePA", quietly = TRUE)

# --- Lectura de targets ---
read_targets <- function(path) {
  if (!file.exists(path)) stop(sprintf("Targets file not found: %s", path))
  ext <- tolower(file_ext(path))
  if (ext %in% c("csv", "tsv", "txt")) {
    # Heurística de columnas
    df <- tryCatch({
      if (ext == "csv") readr::read_csv(path, show_col_types = FALSE)
      else readr::read_tsv(path, show_col_types = FALSE)
    }, error = function(e) {
      # fallback texto simple por línea
      tibble::tibble(target_symbol = readLines(path, warn = FALSE))
    })
    cols <- tolower(names(df))
    pick <- NULL
    for (cand in c("target_symbol", "symbol", "gene", "genes", "gene_symbol")) {
      if (cand %in% cols) { pick <- names(df)[match(cand, cols)]; break }
    }
    if (is.null(pick)) {
      # si el archivo validated_target_genes_unique.csv viene con una sola columna
      if (ncol(df) == 1) {
        pick <- names(df)[1]
      } else {
        stop("Could not infer a symbol column; expected 'target_symbol' or similar.")
      }
    }
    vec <- df[[pick]] %>% as.character() %>% trimws()
    vec <- vec[vec != "" & !is.na(vec)]
    unique(vec)
  } else {
    # texto plano
    vec <- readLines(path, warn = FALSE)
    vec <- trimws(vec)
    vec <- vec[vec != ""]
    unique(vec)
  }
}

symbols <- read_targets(targets_file)
if (length(symbols) == 0) stop("No symbols found in targets file.")
log_line(sprintf("Unique target symbols: %d", length(symbols)))

# --- Background (opcional) ---
bg_symbols <- character(0)
if (nzchar(bg_file)) {
  try({
    bg_symbols <- read_targets(bg_file)
    log_line(sprintf("Background symbols: %d", length(bg_symbols)))
  }, silent = TRUE)
}

# --- Mapear a ENTREZ ---
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))
map_symbols <- function(syms) {
  if (!length(syms)) return(tibble::tibble(SYMBOL = character(), ENTREZID = character()))
  bitr(syms, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    distinct()
}

gene_map <- map_symbols(symbols)
if (!nrow(gene_map)) stop("None of the symbols could be mapped to ENTREZ IDs.")
entrez <- unique(gene_map$ENTREZID)
log_line(sprintf("Mapped to ENTREZ IDs: %d", length(entrez)))

bg_entrez <- NULL
if (length(bg_symbols)) {
  bg_map <- map_symbols(bg_symbols)
  if (nrow(bg_map)) {
    bg_entrez <- unique(bg_map$ENTREZID)
    log_line(sprintf("Mapped background ENTREZ IDs: %d", length(bg_entrez)))
  } else {
    log_line("WARN: Could not map any background symbols; proceeding without universe.")
  }
}

# --- Enrichment: KEGG ---
log_line("Running enrichKEGG ...")
ekegg <- tryCatch(
  enrichKEGG(
    gene            = entrez,
    organism        = org_kegg,         # e.g., 'hsa'
    keyType         = "kegg",
    pAdjustMethod   = "BH",
    pvalueCutoff    = p_cut,
    qvalueCutoff    = q_cut,
    universe        = bg_entrez,
    minGSSize       = minGS,
    maxGSSize       = maxGS
  ),
  error = function(e) { log_line(sprintf("WARN: enrichKEGG failed: %s", e$message)); NULL }
)

# --- Enrichment: GO Biological Process ---
log_line("Running enrichGO (BP) ...")
ego <- tryCatch(
  enrichGO(
    gene            = entrez,
    OrgDb           = org.Hs.eg.db,
    keyType         = "ENTREZID",
    ont             = "BP",
    pAdjustMethod   = "BH",
    pvalueCutoff    = p_cut,
    qvalueCutoff    = q_cut,
    universe        = bg_entrez,
    minGSSize       = minGS,
    maxGSSize       = maxGS,
    readable        = TRUE
  ),
  error = function(e) { log_line(sprintf("WARN: enrichGO failed: %s", e$message)); NULL }
)

# --- Enrichment: Reactome (opcional) ---
ereact <- NULL
if (have_reactome) {
  log_line("Running ReactomePA::enrichPathway ...")
  ereact <- tryCatch(
    ReactomePA::enrichPathway(
      gene            = entrez,
      organism        = "human",
      pAdjustMethod   = "BH",
      pvalueCutoff    = p_cut,
      qvalueCutoff    = q_cut,
      universe        = bg_entrez,
      minGSSize       = minGS,
      maxGSSize       = maxGS,
      readable        = TRUE
    ),
    error = function(e) { log_line(sprintf("WARN: Reactome enrichment failed: %s", e$message)); NULL }
  )
} else {
  log_line("ReactomePA not installed; skipping Reactome enrichment.")
}

# --- Guardar tablas ---
save_table <- function(obj, path) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    log_line(sprintf("No results to save for %s", basename(path)))
    return(invisible(FALSE))
  }
  df <- as.data.frame(obj)
  readr::write_csv(df, path)
  log_line(sprintf("Saved: %s (rows=%d)", path, nrow(df)))
  invisible(TRUE)
}

kegg_csv <- file.path(ENRICH_DIR, "kegg.csv")
go_csv   <- file.path(ENRICH_DIR, "go_bp.csv")
react_csv<- file.path(ENRICH_DIR, "reactome.csv")

save_table(ekegg, kegg_csv)
save_table(ego,   go_csv)
if (!is.null(ereact)) save_table(ereact, react_csv)

# --- Figuras (dotplot) ---
plot_dot <- function(obj, outfile, title = NULL, show = 20) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(invisible(FALSE))
  p <- try(enrichplot::dotplot(obj, showCategory = show) + ggtitle(title %||% ""), silent = TRUE)
  if (inherits(p, "try-error")) return(invisible(FALSE))
  try(ggsave(outfile, p, width = 7, height = 5, dpi = 300), silent = TRUE)
  if (file.exists(outfile)) log_line(sprintf("Saved figure: %s", outfile))
  invisible(TRUE)
}
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

if (!is.null(ekegg))  plot_dot(ekegg,  file.path(FIG_DIR, "kegg_dotplot.png"),  "KEGG Enrichment")
if (!is.null(ego))    plot_dot(ego,    file.path(FIG_DIR, "go_bp_dotplot.png"), "GO Biological Process")
if (!is.null(ereact)) plot_dot(ereact, file.path(FIG_DIR, "reactome_dotplot.png"), "Reactome Pathways")

# --- Resumen breve ---
summarize_top <- function(obj, top = 5) {
  if (is.null(obj)) return(invisible(NULL))
  df <- as.data.frame(obj)
  if (!nrow(df)) return(invisible(NULL))
  head(df[, c("ID", "Description", "GeneRatio", "p.adjust", "Count")], top) %>% print(row.names = FALSE)
}
log_line("Top KEGG terms:")
try(summarize_top(ekegg), silent = TRUE)
log_line("Top GO:BP terms:")
try(summarize_top(ego), silent = TRUE)
if (!is.null(ereact)) {
  log_line("Top Reactome terms:")
  try(summarize_top(ereact), silent = TRUE)
}

# --- Fin ---
log_line("Script finalizado correctamente.")
