#!/usr/bin/env Rscript
# =============================================================================
# Script : 07_networks_clin/ibd_target_overlap.R
# Purpose: Assess whether candidate miRNA validated targets are enriched in
#          known IBD (UC / Crohn's / IBD) disease genes via:
#
#   1. OpenTargets Platform API (GraphQL, free, no auth required).
#      Queries three EFO terms: UC, Crohn's disease, IBD.
#      Falls back to a curated embedded gene list if the API is unavailable.
#
#   2. Fisher's exact test per miRNA — are this miRNA's validated targets
#      more IBD-associated than expected by chance (background = all
#      validated targets in the multiMiR output)?
#
#  Outputs  (DEA_results_round_2/clinical_outputs/)
#  ─────────────────────────────────────────────────
#  ibd_target_overlap.csv       per-miRNA: IBD targets, enrichment OR, FDR
#  ibd_target_genes.csv         miRNA × IBD gene membership (long format)
#  ibd_gene_list.csv            the IBD gene list actually used + source
#  figures/ibd_enrichment.png   dot/bar chart of significant miRNAs
#
# =============================================================================

suppressWarnings(suppressMessages({ library(utils); library(tools) }))

ts      <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse="")))

# ── paths ─────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(getwd())
}
SCRIPT_DIR <- tryCatch(get_script_dir(), error = function(e) normalizePath(getwd()))

# Run from the dataset root directory — finds files recursively
RUN_DIR  <- normalizePath(getwd())

find_file <- function(name) {
  hits <- list.files(RUN_DIR, pattern = paste0("^", name, "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop(paste("Cannot find", name, "— run from the dataset root directory."))
  if (length(hits) > 1) log_msg("WARN: multiple matches for ", name, " — using first: ", hits[1])
  hits[1]
}

VAL_TARGETS <- find_file("targets_validated.csv")
OUT_DIR     <- file.path(RUN_DIR, "ibd_overlap_outputs")
FIG_DIR     <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ_THRESH    <- 0.05       # Fisher's FDR threshold for reporting
OT_PAGE_SIZE   <- 2500       # OpenTargets max per-page result
OT_SCORE_MIN   <- 0.0        # keep all associations (including weak evidence)
OT_TIMEOUT_S   <- 30         # seconds before giving up on API call

log_msg("ibd_target_overlap.R starting")

suppressWarnings(suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
}))

# ── 1. Load validated targets ─────────────────────────────────────────────────
log_msg("Reading validated targets: ", VAL_TARGETS)
val <- read.csv(VAL_TARGETS, stringsAsFactors = FALSE)

# Keep only columns we need
val <- val %>%
  filter(!is.na(mirna), nchar(trimws(mirna)) > 0,
         !is.na(target_symbol), nchar(trimws(target_symbol)) > 0) %>%
  mutate(mirna         = trimws(mirna),
         target_symbol = trimws(target_symbol)) %>%
  select(mirna, target_symbol) %>%
  distinct()

mirna_targets <- split(val$target_symbol, val$mirna)  # list: mirna → character vector
universe_genes <- unique(val$target_symbol)
n_universe     <- length(universe_genes)
log_msg(sprintf("  Loaded %d unique mirna-gene pairs  (%d miRNAs, %d unique genes)",
                nrow(val), length(mirna_targets), n_universe))

# ── 2. Fetch IBD gene list from OpenTargets API ────────────────────────────────

# EFO identifiers:
#   EFO_0003767  Inflammatory bowel disease (IBD)
#   EFO_0000729  Ulcerative colitis (UC)
#   EFO_0000384  Crohn's disease (CD)
OPENTARGETS_URL <- "https://api.platform.opentargets.org/api/v4/graphql"

fetch_ot_genes <- function(efo_id, label = efo_id) {
  api_available <- requireNamespace("httr",     quietly = TRUE) &&
                   requireNamespace("jsonlite", quietly = TRUE)
  if (!api_available) {
    log_msg(sprintf("  WARN: httr/jsonlite not installed — cannot query OpenTargets for %s", label))
    return(character(0))
  }

  query <- sprintf(
    '{ disease(efoId: "%s") {
        name
        associatedTargets(page: {index: 0, size: %d}) {
          rows { target { approvedSymbol } score }
        }
      }
    }', efo_id, OT_PAGE_SIZE)

  resp <- tryCatch(
    httr::POST(
      url     = OPENTARGETS_URL,
      httr::content_type_json(),
      body    = jsonlite::toJSON(list(query = query), auto_unbox = TRUE),
      httr::timeout(OT_TIMEOUT_S)
    ),
    error = function(e) NULL
  )

  if (is.null(resp) || httr::status_code(resp) != 200) {
    log_msg(sprintf("  WARN: OpenTargets API request failed for %s (HTTP %s)",
                    label,
                    if (is.null(resp)) "no response" else httr::status_code(resp)))
    return(character(0))
  }

  parsed  <- tryCatch(
    jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                       simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(parsed)) return(character(0))

  rows <- tryCatch(
    parsed$data$disease$associatedTargets$rows,
    error = function(e) NULL
  )
  if (is.null(rows) || !is.data.frame(rows)) return(character(0))

  genes <- rows$target$approvedSymbol
  genes <- genes[!is.na(genes) & nchar(genes) > 0]
  log_msg(sprintf("  OpenTargets %s ('%s'): %d gene associations retrieved",
                  efo_id, label, length(genes)))
  genes
}

log_msg("Fetching IBD gene associations from OpenTargets Platform ...")
genes_ibd <- fetch_ot_genes("EFO_0003767", "IBD")
genes_uc  <- fetch_ot_genes("EFO_0000729", "Ulcerative colitis")
genes_cd  <- fetch_ot_genes("EFO_0000384", "Crohn's disease")

ibd_genes_api <- unique(c(genes_ibd, genes_uc, genes_cd))
log_msg(sprintf("  Combined unique IBD genes from OpenTargets: %d", length(ibd_genes_api)))

# ── 3. Curated fallback gene list ─────────────────────────────────────────────
# Core IBD susceptibility genes from landmark GWAS meta-analyses
# (Jostins et al. 2012 Nature; de Lange et al. 2017 Nat Genet;
#  Liu et al. 2015 Nat Genet; Anderson et al. 2011 Nat Genet)
curated_ibd_genes <- c(
  # Innate immunity / autophagy / barrier
  "NOD2", "ATG16L1", "IRGM", "CARD9", "LRRK2", "RIPK2",
  # IL-23 / Th17 axis
  "IL23R", "IL12B", "JAK2", "STAT3", "TYK2", "RORC",
  # IL-10 pathway
  "IL10", "IL10RA", "IL10RB",
  # TNF / NF-κB
  "TNF", "TNFSF15", "TNFRSF14", "TRAF3IP2",
  # Epithelial barrier
  "CDH1", "ECM1", "LAMB1", "ITLN1", "GNA12", "HNF4A",
  # Adaptive immunity
  "PTPN22", "PTPN2", "TAGAP", "IL2", "IL21", "IL2RA",
  # MST1 / complement
  "MST1", "C1orf106",
  # Metabolism / cytokine signalling
  "IKZF1", "IRF5", "SMAD3", "PRDM1", "NKX2-3", "KIF21B",
  # UC-specific loci
  "TNFRSF9", "IL16", "PIM1", "PIAS3", "DAP", "GNA12", "BTNL2",
  # CD-specific loci
  "CDKAL1", "BST1", "PTGER4", "LYRM4", "IBD5",
  # Additional well-replicated hits
  "C11orf30", "ZMIZ1", "SLC9A4", "SP140", "ICOSLG",
  "TNFRSF6B", "FCGR2A", "FCGR3A", "GPR35"
)

# ── 4. Decide which gene list to use ──────────────────────────────────────────
if (length(ibd_genes_api) >= 50) {
  ibd_genes <- ibd_genes_api
  source_label <- "OpenTargets Platform (EFO_0003767 + EFO_0000729 + EFO_0000384)"
  log_msg(sprintf("Using OpenTargets gene list (n=%d genes).", length(ibd_genes)))
} else {
  ibd_genes    <- curated_ibd_genes
  source_label <- "Curated GWAS meta-analysis list (Jostins 2012, de Lange 2017)"
  log_msg(sprintf("API returned < 50 genes — falling back to curated list (n=%d genes).",
                  length(ibd_genes)))
}

readr::write_csv(
  data.frame(gene = ibd_genes, source = source_label),
  file.path(OUT_DIR, "ibd_gene_list.csv")
)
log_msg("Saved: ibd_gene_list.csv")

# ── 5. Per-miRNA Fisher's exact test ──────────────────────────────────────────
log_msg("Running Fisher's exact tests ...")

ibd_in_universe <- intersect(ibd_genes, universe_genes)
n_ibd_univ      <- length(ibd_in_universe)
n_non_ibd_univ  <- n_universe - n_ibd_univ

log_msg(sprintf("  IBD genes in multiMiR universe: %d / %d (%.1f%%)",
                n_ibd_univ, n_universe, 100 * n_ibd_univ / n_universe))

results <- lapply(names(mirna_targets), function(mir) {
  tgts  <- unique(mirna_targets[[mir]])
  n_mir <- length(tgts)
  if (n_mir == 0) return(NULL)

  overlap <- intersect(tgts, ibd_in_universe)
  a       <- length(overlap)              # mir targets ∩ IBD genes
  b       <- n_mir - a                    # mir targets, not IBD
  c_      <- n_ibd_univ - a              # IBD genes, not mir targets
  d       <- n_non_ibd_univ - b          # neither

  if (a == 0) {
    return(data.frame(
      miRNA             = mir,
      n_validated       = n_mir,
      n_ibd_overlap     = 0L,
      pct_ibd           = 0,
      OR                = 0,
      p_fisher          = 1,
      padj              = NA_real_,
      overlap_genes     = NA_character_,
      stringsAsFactors  = FALSE
    ))
  }

  ft <- tryCatch(
    fisher.test(matrix(c(a, b, c_, d), nrow = 2, ncol = 2),
                alternative = "greater"),
    error = function(e) NULL
  )
  if (is.null(ft)) return(NULL)

  data.frame(
    miRNA            = mir,
    n_validated      = n_mir,
    n_ibd_overlap    = a,
    pct_ibd          = round(100 * a / n_mir, 2),
    OR               = round(unname(ft$estimate), 3),
    p_fisher         = ft$p.value,
    padj             = NA_real_,
    overlap_genes    = paste(sort(overlap), collapse = ";"),
    stringsAsFactors = FALSE
  )
})

fisher_tbl <- do.call(rbind, Filter(Negate(is.null), results))
fisher_tbl <- fisher_tbl %>%
  mutate(padj = p.adjust(p_fisher, method = "BH")) %>%
  arrange(padj, desc(OR))

readr::write_csv(fisher_tbl, file.path(OUT_DIR, "ibd_target_overlap.csv"))
log_msg(sprintf("Saved: ibd_target_overlap.csv (%d miRNAs, %d with padj <= %.2f)",
                nrow(fisher_tbl),
                sum(fisher_tbl$padj <= PADJ_THRESH, na.rm = TRUE),
                PADJ_THRESH))

# ── 6. Long-format gene table ─────────────────────────────────────────────────
sig_mirnas   <- fisher_tbl %>% filter(padj <= PADJ_THRESH | n_ibd_overlap >= 3)

if (nrow(sig_mirnas) > 0) {
  gene_long <- lapply(seq_len(nrow(sig_mirnas)), function(i) {
    mir   <- sig_mirnas$miRNA[i]
    genes <- strsplit(sig_mirnas$overlap_genes[i], ";")[[1]]
    genes <- genes[!is.na(genes) & nchar(genes) > 0]
    if (!length(genes)) return(NULL)
    data.frame(miRNA = mir, ibd_gene = genes, stringsAsFactors = FALSE)
  })
  gene_long_df <- do.call(rbind, Filter(Negate(is.null), gene_long))
  if (!is.null(gene_long_df) && nrow(gene_long_df) > 0) {
    readr::write_csv(gene_long_df, file.path(OUT_DIR, "ibd_target_genes.csv"))
    log_msg(sprintf("Saved: ibd_target_genes.csv (%d rows)", nrow(gene_long_df)))
  }
}

# ── 7. Figure — enrichment bar chart ─────────────────────────────────────────
plot_tbl <- fisher_tbl %>%
  filter(n_ibd_overlap > 0) %>%
  arrange(padj, desc(OR)) %>%
  head(20) %>%
  mutate(
    miRNA      = factor(miRNA, levels = rev(miRNA)),
    sig_label  = case_when(
      padj <= 0.01 ~ "**",
      padj <= 0.05 ~ "*",
      TRUE         ~ ""
    ),
    bar_fill   = if_else(padj <= PADJ_THRESH, "#d62728", "#aec7e8")
  )

if (nrow(plot_tbl) > 0) {
  p <- ggplot(plot_tbl, aes(x = miRNA, y = n_ibd_overlap, fill = bar_fill)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig_label), hjust = -0.3, size = 4.5) +
    scale_fill_identity(
      guide  = "legend",
      labels = c("#d62728" = "FDR ≤ 0.05", "#aec7e8" = "FDR > 0.05"),
      breaks = c("#d62728", "#aec7e8"),
      name   = NULL
    ) +
    coord_flip() +
    labs(
      title    = sprintf("Validated targets overlapping IBD genes\n(source: %s)", source_label),
      subtitle = sprintf("* FDR ≤ 0.05   ** FDR ≤ 0.01   (Fisher's exact, BH correction)\nUniverse: %d validated target genes; IBD gene set: %d genes",
                         n_universe, length(ibd_genes)),
      x        = NULL,
      y        = "# IBD-associated validated targets"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title      = element_text(size = 11, face = "bold"))

  fig_h <- max(4, 0.4 * nrow(plot_tbl) + 2)
  ggsave(file.path(FIG_DIR, "ibd_enrichment.png"), p,
         width = 8, height = fig_h, dpi = 300)
  log_msg("Saved: figures/ibd_enrichment.png")
}

# ── 8b. Validated vs Predicted vs Baseline comparison ─────────────────────────
# Re-run the IDENTICAL Fisher test (same IBD gene list) over the ALL-miRNA
# baseline target tables produced by 17_multimir_targets_baseline.R, for both
# validated and predicted evidence, each with its own pooled-target universe.
# Significant miRNAs (those that fed the canonical analysis above) are flagged so
# the dashboard can slice the three layers the thesis compares:
#   validated      = is_significant & evidence == "validated"
#   non-validated  = is_significant & evidence == "predicted"
#   baseline       = all miRNAs (the full distribution), either evidence
# The canonical ibd_target_overlap.csv above is left untouched; this writes a
# separate ibd_target_overlap_comparison.csv. If the baseline files are missing
# (script 15 not run yet) the comparison is skipped without failing the run.

find_file_opt <- function(name) {
  hits <- list.files(RUN_DIR, pattern = paste0("^", name, "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) return(NA_character_)
  if (length(hits) > 1) log_msg("  comparison: multiple matches for ", name, " — using first: ", hits[1])
  hits[1]
}

run_fisher_layer <- function(target_csv, evidence_label, ibd_genes) {
  if (is.na(target_csv) || !file.exists(target_csv)) {
    log_msg(sprintf("  comparison: %s baseline targets not found — skipping that layer", evidence_label))
    return(NULL)
  }
  df <- tryCatch(read.csv(target_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !nrow(df) || !all(c("mirna", "target_symbol") %in% names(df))) {
    log_msg(sprintf("  comparison: %s file unusable — skipping that layer", evidence_label))
    return(NULL)
  }
  df <- df %>%
    filter(!is.na(mirna), nchar(trimws(mirna)) > 0,
           !is.na(target_symbol), nchar(trimws(target_symbol)) > 0) %>%
    mutate(mirna = trimws(mirna), target_symbol = trimws(target_symbol)) %>%
    select(mirna, target_symbol) %>%
    distinct()

  mt       <- split(df$target_symbol, df$mirna)
  universe <- unique(df$target_symbol)
  n_uni    <- length(universe)
  ibd_uni  <- intersect(ibd_genes, universe)
  n_ibd    <- length(ibd_uni)
  n_non    <- n_uni - n_ibd
  log_msg(sprintf("  comparison [%s]: %d miRNAs, universe=%d genes, IBD-in-universe=%d",
                  evidence_label, length(mt), n_uni, n_ibd))

  rows <- lapply(names(mt), function(mir) {
    tg <- unique(mt[[mir]]); nm <- length(tg)
    if (nm == 0) return(NULL)
    ov <- intersect(tg, ibd_uni)
    a  <- length(ov); b <- nm - a; c_ <- n_ibd - a; d <- n_non - b
    if (a == 0) {
      return(data.frame(miRNA = mir, evidence = evidence_label, n_targets = nm,
                        n_ibd_overlap = 0L, pct_ibd = 0, OR = 0, p_fisher = 1,
                        padj = NA_real_, overlap_genes = NA_character_,
                        stringsAsFactors = FALSE))
    }
    ft <- tryCatch(fisher.test(matrix(c(a, b, c_, d), nrow = 2, ncol = 2),
                               alternative = "greater"),
                   error = function(e) NULL)
    if (is.null(ft)) return(NULL)
    data.frame(miRNA = mir, evidence = evidence_label, n_targets = nm,
               n_ibd_overlap = a, pct_ibd = round(100 * a / nm, 2),
               OR = round(unname(ft$estimate), 3), p_fisher = ft$p.value,
               padj = NA_real_, overlap_genes = paste(sort(ov), collapse = ";"),
               stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  tbl$universe_size   <- n_uni
  tbl$ibd_in_universe <- n_ibd
  tbl %>% mutate(padj = p.adjust(p_fisher, method = "BH"))
}

log_msg("Building validated/predicted/baseline comparison ...")
base_val  <- find_file_opt("targets_validated_baseline.csv")
base_pred <- find_file_opt("targets_predicted_baseline.csv")

if (is.na(base_val) && is.na(base_pred)) {
  log_msg("Baseline target files not found — run 06_targets_enrich/17_multimir_targets_baseline.R "
          , "first to enable the comparison. Skipping ibd_target_overlap_comparison.csv.")
} else {
  cmp <- bind_rows(
    run_fisher_layer(base_val,  "validated", ibd_genes),
    run_fisher_layer(base_pred, "predicted", ibd_genes)
  )
  if (!is.null(cmp) && nrow(cmp)) {
    # `val` (loaded at step 1) holds the significant validated targets — its
    # miRNAs are exactly the candidates that fed the canonical analysis.
    sig_set <- unique(val$mirna)
    cmp <- cmp %>%
      mutate(is_significant = miRNA %in% sig_set) %>%
      arrange(evidence, padj, desc(OR))
    out_cmp <- file.path(OUT_DIR, "ibd_target_overlap_comparison.csv")
    readr::write_csv(cmp, out_cmp)
    log_msg(sprintf("Saved: ibd_target_overlap_comparison.csv (%d rows; %d significant miRNA-rows; evidence: %s)",
                    nrow(cmp), sum(cmp$is_significant),
                    paste(sort(unique(cmp$evidence)), collapse = ", ")))
  } else {
    log_msg("Comparison produced no rows — check the baseline target files.")
  }
}

# ── 8. Summary ────────────────────────────────────────────────────────────────
log_msg("\n========== SUMMARY ==========")
log_msg(sprintf("IBD gene list used    : %d genes  (%s)", length(ibd_genes), source_label))
log_msg(sprintf("miRNAs tested         : %d", nrow(fisher_tbl)))
log_msg(sprintf("With any IBD overlap  : %d", sum(fisher_tbl$n_ibd_overlap > 0)))
log_msg(sprintf("FDR <= %.2f           : %d", PADJ_THRESH,
                sum(fisher_tbl$padj <= PADJ_THRESH, na.rm = TRUE)))
if (nrow(fisher_tbl) > 0 && any(fisher_tbl$n_ibd_overlap > 0)) {
  top3 <- fisher_tbl %>% filter(n_ibd_overlap > 0) %>% head(3)
  for (i in seq_len(nrow(top3))) {
    r <- top3[i, ]
    log_msg(sprintf("  %-30s  overlap=%d  OR=%.2f  padj=%.3g",
                    r$miRNA, r$n_ibd_overlap, r$OR, r$padj))
  }
}
log_msg(sprintf("Outputs saved to: %s", OUT_DIR))
log_msg("ibd_target_overlap.R completed.")
