#!/usr/bin/env Rscript
# ===========================================
# Script: 19_interactive_network.R
# Description: Production tripartite miRNA -> Gene -> Pathway/GO dashboards (Cytoscape.js).
#
#  v3 — performance + biology + aesthetics overhaul
#  ------------------------------------------------
#  WHY THE REWRITE: v2 built ~2,500 nodes / ~24,000 edges per dashboard because it kept
#  every multiMiR "validated" row (including weak / TarBase-negative evidence). That is a
#  hairball: slow to render and biologically unreadable. v3 fixes this on three axes.
#
#  BIOLOGY (signal, not noise):
#   - Keep only STRONG experimentally-validated interactions (support_type == "Functional MTI").
#   - Drop rows with no gene symbol.
#   - Hub-focused default view: genes targeted by >= MIN_MIRNAS_PER_GENE significant miRNAs
#     (live-adjustable slider in the UI, down to 1 to reveal the full strong-evidence network).
#   - Edge weight encodes evidence depth (# distinct PMIDs supporting the interaction).
#
#  OPTIMIZATION:
#   - Cytoscape.js with a ONE-SHOT fcose layout — no continuous physics simulation.
#     Filtering only toggles node visibility (edges auto-hide), so it never re-simulates.
#   - Straight edges, deterministic layout, compact embedded JSON.
#
#  AESTHETICS / PRODUCTION:
#   - Tripartite-aware fcose + concentric (ringed) layouts, live stats header, node search,
#     degree slider, neighborhood focus-on-click, group filters, PNG export, polished legend.
#   - GO pathway nodes tinted by ontology (BP / CC / MF).
# ===========================================

suppressWarnings(suppressMessages({
  library(readr)
  library(jsonlite)
}))

RUN_DIR   <- getwd()
args      <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else RUN_DIR
OUT_DIR   <- file.path(RUN_DIR, "network_outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── tunable parameters ──────────────────────────────────────────────────────────
STRONG_SUPPORT      <- c("Functional MTI")  # multiMiR support_type(s) treated as strong
MAX_PATHWAYS        <- 20                    # cap on enriched terms per dashboard
MIN_MIRNAS_PER_GENE <- 2                     # default UI slider value (hub view); data keeps >= 1
HUB_THRESHOLD       <- 3                     # gene targeted by >= this many miRNAs -> "hub" highlight
LFC_THRESH          <- 0.58
PADJ_THRESH         <- 0.05

log_line <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse="")))
log_line("Starting tripartite dashboard generation (v3, Cytoscape.js)...")

# ── helpers ───────────────────────────────────────────────────────────────────
find_file <- function(pattern, dir, required = TRUE) {
  files <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    if (required) stop(sprintf("Cannot find file matching: %s", pattern))
    return(NULL)
  }
  files[1]
}

rescale_range <- function(x, lo, hi) {
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || r[1] == r[2]) return(rep((lo + hi) / 2, length(x)))
  lo + (x - r[1]) / (r[2] - r[1]) * (hi - lo)
}

# Expand "gene1/gene2/..." rows into long Pathway/Gene pairs (base-R separate_rows)
expand_terms <- function(sig_terms) {
  do.call(rbind, lapply(seq_len(nrow(sig_terms)), function(i) {
    genes <- strsplit(sig_terms$geneID[i], "/", fixed = TRUE)[[1]]
    genes <- genes[nzchar(genes)]
    if (length(genes) == 0) return(NULL)
    data.frame(Pathway = sig_terms$Description[i], Gene = genes, stringsAsFactors = FALSE)
  }))
}

# ── load data ─────────────────────────────────────────────────────────────────
# DEA: find ALL matching files — handles 1 or N contrasts automatically
dea_files <- list.files(input_dir, pattern = "^DEA_Default.*\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(dea_files) == 0) stop(sprintf("No DEA_Default*.csv found under: %s", input_dir))
log_line(sprintf("  DEA file(s) found: %d — %s",
                 length(dea_files), paste(basename(dea_files), collapse = ", ")))

dea_list <- lapply(dea_files, function(f) {
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) { log_line(sprintf("  WARN: could not read %s", basename(f))); return(NULL) }
  df$contrast <- tools::file_path_sans_ext(basename(f))
  df
})
dea_list <- Filter(Negate(is.null), dea_list)
if (length(dea_list) == 0) stop("All DEA files failed to load.")
dea_all <- as.data.frame(do.call(rbind, dea_list), stringsAsFactors = FALSE)

val_df  <- read_csv(find_file("^targets_validated\\.csv$", input_dir), show_col_types = FALSE)
kegg_df <- read_csv(find_file("^KEGG_Enrichment_Results\\.csv$", input_dir), show_col_types = FALSE)
go_path <- find_file("^GO_Enrichment_Results\\.csv$", input_dir, required = FALSE)
go_df   <- if (!is.null(go_path)) read_csv(go_path, show_col_types = FALSE) else NULL

val_df <- as.data.frame(val_df, stringsAsFactors = FALSE)

# Strong, gene-annotated interactions only
strong <- val_df[val_df$support_type %in% STRONG_SUPPORT &
                 !is.na(val_df$target_symbol) & nzchar(val_df$target_symbol), ]
strong$pmid <- as.character(strong$pmid)
log_line(sprintf("  Strong-evidence interactions (%s): %d rows / %d unique pairs",
                 paste(STRONG_SUPPORT, collapse=","), nrow(strong),
                 nrow(unique(strong[, c("mirna","target_symbol")]))))

# Evidence depth: # distinct PMIDs per (miRNA, gene)
strong$key <- paste(strong$mirna, strong$target_symbol, sep = "||")
ev_pmids   <- tapply(strong$pmid, strong$key,
                     function(x) length(unique(x[!is.na(x) & nzchar(x)])))

# ── significant miRNAs — works for 1 or N contrasts ──────────────────────────
# Strip "DEA_Default_" prefix from contrast names for clean display
clean_contrast <- function(x) sub("^DEA_Default_", "", x)

sig_all <- dea_all[!is.na(dea_all$padj) & dea_all$padj <= PADJ_THRESH &
                   !is.na(dea_all$log2FoldChange) & abs(dea_all$log2FoldChange) >= LFC_THRESH, ]
if (nrow(sig_all) == 0) {
  log_line("WARN: No miRNAs pass thresholds — falling back to top 30 by |S_miRNA| across all contrasts.")
  ord <- order(-abs(ifelse(is.na(dea_all$S_miRNA), 0, dea_all$S_miRNA)))
  sig_all <- dea_all[head(ord, 30), ]
}

# Best row per miRNA (lowest padj = most significant contrast wins)
sig_all       <- sig_all[order(sig_all$padj), ]
sig_mirnas_df <- sig_all[!duplicated(sig_all$miRNA), ]

# Contrast metadata — used by build_elements and write_dashboard
n_dea          <- length(dea_files)
contrast_names <- sort(unique(sig_mirnas_df$contrast))           # primary contrasts present
# All contrasts each miRNA is significant in (may be >1 when n_dea > 1)
contrast_all_map <- tapply(sig_all$contrast, sig_all$miRNA,
                           function(x) paste(clean_contrast(sort(unique(x))), collapse = " + "))
# 0-based index into palette for JS
contrast_idx_map <- setNames(seq_along(contrast_names) - 1L, contrast_names)
sig_mirnas_df$contrastIdx <- contrast_idx_map[sig_mirnas_df$contrast]

log_line(sprintf("  Significant miRNAs: %d across %d contrast(s): %s",
                 nrow(sig_mirnas_df), n_dea,
                 paste(clean_contrast(contrast_names), collapse = ", ")))

# ── build Cytoscape element payload ─────────────────────────────────────────────
build_elements <- function(enrich_df, label, onto_col = NULL) {
  enrich_df <- as.data.frame(enrich_df, stringsAsFactors = FALSE)

  sig_terms <- enrich_df[!is.na(enrich_df$p.adjust) & enrich_df$p.adjust <= PADJ_THRESH, ]
  sig_terms <- sig_terms[order(sig_terms$p.adjust), ]
  if (nrow(sig_terms) == 0) {
    log_line(sprintf("  WARN (%s): no significant terms — using top 10.", label))
    sig_terms <- head(enrich_df[order(enrich_df$p.adjust), ], 10)
  } else if (nrow(sig_terms) > MAX_PATHWAYS) {
    sig_terms <- sig_terms[seq_len(MAX_PATHWAYS), ]
  }
  log_line(sprintf("  %s: %d significant terms", label, nrow(sig_terms)))

  te <- expand_terms(sig_terms)
  if (is.null(te) || nrow(te) == 0) return(NULL)

  # miRNA -> gene pairs: strong evidence, gene in an enriched pathway, miRNA significant
  pairs <- strong[strong$target_symbol %in% unique(te$Gene) &
                  strong$mirna %in% sig_mirnas_df$miRNA, ]
  if (nrow(pairs) == 0) return(NULL)
  tf <- unique(pairs[, c("mirna", "target_symbol")])
  names(tf) <- c("miRNA", "Gene")

  # Keep ALL strong-evidence genes in the data; the UI slider hides singletons by default.
  n_mir_per_gene <- table(tf$Gene)
  tf$n_mirnas    <- as.integer(n_mir_per_gene[tf$Gene])

  genes  <- sort(unique(tf$Gene))
  mirnas <- sort(unique(tf$miRNA))
  te     <- te[te$Gene %in% genes, ]
  paths  <- unique(te$Pathway)
  if (length(paths) == 0) return(NULL)

  # ── per-entity attributes ──────────────────────────────────────────────────
  lfc_map   <- setNames(sig_mirnas_df$log2FoldChange, sig_mirnas_df$miRNA)
  padj_map  <- setNames(sig_mirnas_df$padj,           sig_mirnas_df$miRNA)
  score_map <- setNames(ifelse(is.na(sig_mirnas_df$S_miRNA), 0, sig_mirnas_df$S_miRNA),
                        sig_mirnas_df$miRNA)

  gene_nmir <- setNames(as.integer(n_mir_per_gene[genes]), genes)
  pinfo <- sig_terms[match(paths, sig_terms$Description), ]
  pinfo$neg_log_p <- -log10(pinfo$p.adjust)
  onto <- if (!is.null(onto_col) && onto_col %in% names(pinfo)) as.character(pinfo[[onto_col]])
          else if ("subcategory" %in% names(pinfo)) as.character(pinfo$subcategory)
          else rep("", length(paths))
  onto[is.na(onto)] <- ""

  # Top-25% significance thresholds (lowest padj = most significant)
  gold_mirna <- as.numeric(quantile(padj_map[mirnas], 0.25, na.rm = TRUE))
  gold_path  <- as.numeric(quantile(pinfo$p.adjust,   0.25, na.rm = TRUE))

  # Pre-computed display sizes
  mir_size  <- setNames(rescale_range(abs(score_map[mirnas]), 22, 46), mirnas)
  gene_size <- setNames(rescale_range(gene_nmir, 16, 34), genes)
  path_size <- setNames(rescale_range(pinfo$neg_log_p, 30, 52), paths)

  # ── nodes ──────────────────────────────────────────────────────────────────
  node_list <- list()

  for (m in mirnas) {
    lfc <- lfc_map[[m]]; pj <- padj_map[[m]]; sc <- score_map[[m]]
    dir <- if (is.na(lfc)) "na" else if (lfc > 0) "up" else "down"
    top <- !is.na(pj) && !is.na(gold_mirna) && pj <= gold_mirna
    # Contrast info — primary contrast and (when multi-DEA) all contrasts
    m_row       <- sig_mirnas_df[sig_mirnas_df$miRNA == m, ]
    m_contrast  <- if (nrow(m_row) > 0) m_row$contrast[1] else ""
    m_cidx      <- if (nrow(m_row) > 0) m_row$contrastIdx[1] else 0L
    m_call      <- if (!is.null(contrast_all_map[[m]])) contrast_all_map[[m]] else clean_contrast(m_contrast)
    contrast_tip <- if (n_dea > 1) sprintf("<br>Contrast(s): %s", m_call) else
                    sprintf("<br>Contrast: %s", clean_contrast(m_contrast))
    tip <- sprintf("<b>%s</b><br>LFC: %.3f (%s)<br>padj: %.2e<br>S_miRNA: %.2f<br>Targets shown: %d%s%s",
                   m, lfc, if (identical(dir,"up")) "UP &#9650;" else "DOWN &#9660;",
                   pj, sc, sum(tf$miRNA == m), contrast_tip,
                   if (top) "<br><b>&#9733; Top 25% significance</b>" else "")
    node_list[[length(node_list)+1]] <- list(data = list(
      id = m, label = m, group = "miRNA", dir = dir,
      lfc = round(lfc, 4), padj = pj, size = round(mir_size[[m]], 1),
      topSig = as.integer(top), ndeg = sum(tf$miRNA == m),
      contrast = m_contrast, contrastIdx = as.integer(m_cidx), tip = tip))
  }

  for (g in genes) {
    nm  <- gene_nmir[[g]]
    hub <- nm >= HUB_THRESHOLD
    tip <- sprintf("<b>%s</b><br>Targeted by %d significant miRNA(s)%s",
                   g, nm, if (hub) "<br><b>&#9733; Hub gene</b>" else "")
    node_list[[length(node_list)+1]] <- list(data = list(
      id = g, label = g, group = "Gene", isHub = as.integer(hub),
      nmir = nm, size = round(gene_size[[g]], 1), tip = tip))
  }

  for (i in seq_along(paths)) {
    p <- paths[i]
    pj <- pinfo$p.adjust[i]
    top <- !is.na(pj) && !is.na(gold_path) && pj <= gold_path
    gr  <- if ("GeneRatio" %in% names(pinfo)) as.character(pinfo$GeneRatio[i]) else ""
    cnt <- if ("Count" %in% names(pinfo)) pinfo$Count[i] else NA
    n_in <- length(unique(te$Gene[te$Pathway == p]))
    tip <- sprintf("<b>%s</b><br>p.adjust: %.2e<br>Gene ratio: %s<br>Targets in view: %d%s%s",
                   p, pj, gr, n_in,
                   if (nzchar(onto[i])) sprintf("<br>Category: %s", onto[i]) else "",
                   if (top) "<br><b>&#9733; Top 25% significance</b>" else "")
    node_list[[length(node_list)+1]] <- list(data = list(
      id = p, label = p, group = "Pathway", onto = onto[i],
      size = round(path_size[[p]], 1), topSig = as.integer(top),
      ndeg = n_in, tip = tip))
  }

  # ── edges ──────────────────────────────────────────────────────────────────
  edge_list <- list()
  ev_counts <- pmax(1L, as.integer(ev_pmids[paste(tf$miRNA, tf$Gene, sep = "||")]))
  ev_counts[is.na(ev_counts)] <- 1L
  ew <- rescale_range(log1p(ev_counts), 1.4, 5.5)

  for (i in seq_len(nrow(tf))) {
    m <- tf$miRNA[i]; g <- tf$Gene[i]
    lfc <- lfc_map[[m]]
    dir <- if (is.na(lfc)) "na" else if (lfc > 0) "up" else "down"
    tip <- sprintf("<b>%s &rarr; %s</b><br>miRNA: %s<br>Evidence: %d PMID(s)",
                   m, g, if (identical(dir,"up")) "UP" else "DOWN", ev_counts[i])
    edge_list[[length(edge_list)+1]] <- list(data = list(
      id = sprintf("mg%d", i), source = m, target = g,
      etype = "mg", dir = dir, w = round(ew[i], 2), tip = tip))
  }

  gp <- unique(te[, c("Gene", "Pathway")])
  for (i in seq_len(nrow(gp))) {
    edge_list[[length(edge_list)+1]] <- list(data = list(
      id = sprintf("gp%d", i), source = gp$Gene[i], target = gp$Pathway[i],
      etype = "gp", dir = "na", w = 1,
      tip = sprintf("<b>%s &rarr; %s</b><br>gene member of pathway", gp$Gene[i], gp$Pathway[i])))
  }

  stats <- list(mirna = length(mirnas), gene = length(genes), pathway = length(paths),
                edges = length(edge_list), label = label,
                contrastNames = contrast_names, nDea = n_dea)
  log_line(sprintf("  %s graph: %d miRNAs, %d genes, %d pathways, %d edges",
                   label, stats$mirna, stats$gene, stats$pathway, stats$edges))

  list(json  = toJSON(list(nodes = node_list, edges = edge_list),
                      auto_unbox = TRUE, na = "null", pretty = FALSE),
       stats = stats)
}

# ── HTML dashboard writer (Cytoscape.js) ────────────────────────────────────────
write_dashboard <- function(payload, out_file, page_title, pathway_label, is_go = FALSE) {
  if (is.null(payload)) {
    log_line(sprintf("SKIP: no connected nodes for '%s'.", page_title))
    return(invisible(NULL))
  }
  # Contrast palette — colorblind-safe, distinct from node fill colors
  CONTRAST_PALETTE <- c("#e07b00","#2ca02c","#9467bd","#17becf","#bcbd22","#8c564b")
  cnames_r  <- payload$stats$contrastNames
  n_dea_r   <- payload$stats$nDea
  used_pal  <- CONTRAST_PALETTE[seq_len(min(length(cnames_r), length(CONTRAST_PALETTE)))]

  cfg <- toJSON(list(pageTitle = page_title, pathwayLabel = pathway_label,
                     minMir = MIN_MIRNAS_PER_GENE, hubThreshold = HUB_THRESHOLD,
                     isGO = is_go, nDea = n_dea_r,
                     contrastNames = cnames_r, contrastPalette = used_pal,
                     stats = payload$stats), auto_unbox = TRUE)

  # ── contrast-specific HTML snippets (only rendered when n_dea > 1) ──────────
  if (n_dea_r > 1) {
    contrast_checks <- paste(sapply(seq_along(cnames_r), function(i) {
      col  <- used_pal[i]
      name <- clean_contrast(cnames_r[i])
      sprintf('<label><input type="checkbox" class="chk-contrast" data-contrast="%s" checked> <span style="display:inline-block;width:9px;height:9px;border-radius:50%%;background:%s;border:1.5px solid rgba(0,0,0,0.25);vertical-align:middle;margin-right:3px;"></span>%s</label>',
              cnames_r[i], col, name)
    }), collapse = "\n      ")
    contrast_filter_html <- sprintf(
      '\n      <div class="filter-group">
        <div class="filter-header"><span>Contrasts</span>
          <span><button class="btn-mini" data-all-contrast="1">All</button>
                <button class="btn-mini" data-all-contrast="0">None</button></span></div>
        <div class="checkbox-list">%s</div>
      </div>', contrast_checks)
    contrast_legend_html <- paste(sapply(seq_along(cnames_r), function(i) {
      sprintf('        <div class="legend-row"><span class="sw" style="background:transparent;border:2.5px solid %s"></span> %s</div>',
              used_pal[i], clean_contrast(cnames_r[i]))
    }), collapse = "\n")
  } else {
    contrast_filter_html <- ""
    contrast_legend_html <- ""
  }

  html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>', page_title, '</title>
<script src="https://unpkg.com/cytoscape@3.30.2/dist/cytoscape.min.js"></script>
<script src="https://unpkg.com/layout-base@2.0.1/layout-base.js"></script>
<script src="https://unpkg.com/cose-base@2.2.0/cose-base.js"></script>
<script src="https://unpkg.com/cytoscape-fcose@2.2.0/cytoscape-fcose.js"></script>
<style>
  :root{ --mir-up:#d62728; --mir-down:#1f77b4; --gene:#4DBBD5; --hub:#FF7F00;
         --path:#7B2D8B; --gold:#FFD700; --ink:#2c3e50; --panel:rgba(255,255,255,0.97); }
  *{ box-sizing:border-box; }
  html,body{ margin:0; padding:0; height:100%; font-family:"Segoe UI",Tahoma,sans-serif;
             background:#eef1f4; color:var(--ink); overflow:hidden; }
  #cy{ position:absolute; inset:0; z-index:1; background:
       radial-gradient(circle at 30% 20%, #ffffff 0%, #e7ecf1 100%); }

  #topbar{ position:absolute; top:0; left:0; right:0; height:52px; z-index:20;
           display:flex; align-items:center; gap:14px; padding:0 16px;
           background:linear-gradient(90deg,#1f2a44,#2c3e50); color:#fff;
           box-shadow:0 2px 10px rgba(0,0,0,0.18); }
  #topbar h1{ font-size:15px; font-weight:600; margin:0; white-space:nowrap; }
  #topbar .sub{ font-size:11px; opacity:0.7; font-weight:400; }
  .stat-chips{ display:flex; gap:8px; margin-left:auto; flex-wrap:wrap; }
  .chip{ background:rgba(255,255,255,0.12); border-radius:12px; padding:3px 10px;
         font-size:11px; display:flex; gap:5px; align-items:center; }
  .chip b{ font-size:12px; }
  .chip .dot{ width:8px; height:8px; border-radius:50%; display:inline-block; }

  #panel{ position:absolute; top:64px; left:12px; z-index:15; width:288px;
          max-height:calc(100vh - 76px); background:var(--panel); border-radius:10px;
          box-shadow:0 6px 22px rgba(0,0,0,0.14); display:flex; flex-direction:column;
          overflow:hidden; }
  #panel .scroll{ overflow-y:auto; padding:12px; display:flex; flex-direction:column; gap:12px; }
  .toolbar{ display:flex; flex-wrap:wrap; gap:6px; }
  .toolbar button, .toolbar select{ font-size:11px; padding:5px 8px; cursor:pointer;
       border:1px solid #cdd4da; background:#fff; border-radius:5px; color:var(--ink); }
  .toolbar button:hover{ background:#eef2f6; }
  .search{ width:100%; padding:7px 9px; border:1px solid #cdd4da; border-radius:6px; font-size:12px; }
  .slider-row{ font-size:12px; display:flex; flex-direction:column; gap:4px; }
  .slider-row input[type=range]{ width:100%; }
  .slider-row .val{ font-weight:600; color:var(--mir-up); }
  .forces{ border:1px solid #e0e5ea; border-radius:6px; padding:8px 10px; background:#fbfcfd; }
  .forces summary{ font-size:12px; font-weight:700; cursor:pointer; color:var(--ink); }
  .forces .slider-row{ margin-top:8px; }

  .filter-group{ display:flex; flex-direction:column; }
  .filter-header{ display:flex; justify-content:space-between; align-items:center;
                  margin-bottom:4px; font-size:12px; font-weight:700; }
  .btn-mini{ font-size:10px; padding:1px 6px; cursor:pointer; border:1px solid #cdd4da;
             background:#fff; border-radius:3px; margin-left:3px; }
  .checkbox-list{ border:1px solid #d7dde2; border-radius:5px; padding:5px; max-height:120px;
                  overflow-y:auto; background:#fff; font-size:12px; }
  .checkbox-list label{ display:block; margin-bottom:2px; cursor:pointer; white-space:nowrap;
                        overflow:hidden; text-overflow:ellipsis; }
  .checkbox-list input{ margin-right:6px; }

  .legend{ font-size:11px; border-top:1px solid #e3e8ec; padding-top:8px; }
  .legend b{ font-size:12px; }
  .legend-row{ display:flex; align-items:center; gap:7px; margin:3px 0; }
  .sw{ width:13px; height:13px; border-radius:50%; flex-shrink:0; }
  .sw-hex{ width:14px; height:14px; background:var(--path); flex-shrink:0;
           clip-path:polygon(25% 0,75% 0,100% 50%,75% 100%,25% 100%,0 50%); }
  .sw-line{ width:20px; height:3px; border-radius:2px; flex-shrink:0; }
  .muted{ color:#8a94a0; font-size:10px; margin-top:4px; }

  #tooltip{ position:absolute; z-index:40; pointer-events:none; display:none; max-width:280px;
            background:rgba(28,36,54,0.96); color:#fff; font-size:12px; line-height:1.45;
            padding:8px 11px; border-radius:7px; box-shadow:0 4px 16px rgba(0,0,0,0.3); }
  #tooltip b{ color:#ffe27a; }
  #hint{ position:absolute; bottom:10px; right:14px; z-index:15; font-size:11px; color:#7a838f;
         background:rgba(255,255,255,0.8); padding:4px 10px; border-radius:6px; }
</style>
</head>
<body>
  <div id="topbar">
    <div>
      <h1>', page_title, '</h1>
      <div class="sub">Tripartite regulatory network &middot; strong-evidence miRNA targets</div>
    </div>
    <div class="stat-chips" id="stat-chips"></div>
  </div>

  <div id="panel">
    <div class="scroll">
      <input id="search" class="search" type="text" placeholder="Search miRNA / gene / pathway...">
      <div class="toolbar">
        <select id="layout-sel" title="Layout">
          <option value="fcose">Layout: organic (fCoSE)</option>
          <option value="concentric">Layout: tripartite rings</option>
          <option value="cose">Layout: spring (CoSE)</option>
        </select>
        <button id="btn-relayout" title="Re-run layout on visible nodes">Re-layout</button>
        <button id="btn-fit">Fit</button>
        <button id="btn-reset">Reset view</button>
        <button id="btn-labels">Labels</button>
        <button id="btn-png">Export PNG</button>
      </div>

      <div class="slider-row">
        <label>Min. miRNAs per gene shown: <span class="val" id="deg-val">', MIN_MIRNAS_PER_GENE, '</span></label>
        <input id="deg-slider" type="range" min="1" max="8" step="1" value="', MIN_MIRNAS_PER_GENE, '">
      </div>

      <details class="forces" open>
        <summary>Forces (spread / cluster)</summary>
        <div class="slider-row">
          <label>Repulsion (spread apart): <span class="val" id="rep-val">9000</span></label>
          <input id="f-rep" type="range" min="2000" max="30000" step="500" value="9000">
        </div>
        <div class="slider-row">
          <label>Link distance: <span class="val" id="lnk-val">110</span></label>
          <input id="f-lnk" type="range" min="40" max="320" step="10" value="110">
        </div>
        <div class="slider-row">
          <label>Gravity (pull to centre): <span class="val" id="grv-val">0.18</span></label>
          <input id="f-grv" type="range" min="0" max="0.8" step="0.02" value="0.18">
        </div>
      </details>

      ', contrast_filter_html, '
      <div class="filter-group">
        <div class="filter-header"><span style="color:var(--mir-up)">miRNAs</span>
          <span><button class="btn-mini" data-all="miRNA" data-on="1">All</button>
                <button class="btn-mini" data-all="miRNA" data-on="0">None</button></span></div>
        <div class="checkbox-list" id="list-miRNA"></div>
      </div>
      <div class="filter-group">
        <div class="filter-header"><span style="color:var(--gene)">Genes</span>
          <span><button class="btn-mini" data-all="Gene" data-on="1">All</button>
                <button class="btn-mini" data-all="Gene" data-on="0">None</button></span></div>
        <div class="checkbox-list" id="list-Gene"></div>
      </div>
      <div class="filter-group">
        <div class="filter-header"><span style="color:var(--path)">', pathway_label, '</span>
          <span><button class="btn-mini" data-all="Pathway" data-on="1">All</button>
                <button class="btn-mini" data-all="Pathway" data-on="0">None</button></span></div>
        <div class="checkbox-list" id="list-Pathway"></div>
      </div>

      <div class="legend">
        <b>Legend</b>
        <div class="legend-row"><span class="sw" style="background:var(--mir-up)"></span> miRNA upregulated (LFC &gt; 0)</div>
        <div class="legend-row"><span class="sw" style="background:var(--mir-down)"></span> miRNA downregulated (LFC &lt; 0)</div>
        <div class="legend-row"><span class="sw" style="background:var(--gene)"></span> Target gene</div>
        <div class="legend-row"><span class="sw" style="background:var(--hub)"></span> Hub gene (&#8805;', HUB_THRESHOLD, ' miRNAs)</div>
        <div class="legend-row"><span class="sw-hex"></span> ', pathway_label, ' term</div>
        <div class="legend-row"><span class="sw-line" style="background:var(--mir-up)"></span> Edge: miRNA UP &rarr; gene</div>
        <div class="legend-row"><span class="sw-line" style="background:var(--mir-down)"></span> Edge: miRNA DOWN &rarr; gene</div>
        <div class="legend-row"><span class="sw" style="background:transparent;border:2px solid var(--gold)"></span> Top 25% most significant</div>
', contrast_legend_html, '
        <div class="muted">Node size = significance &middot; edge width = evidence (PMIDs) &middot; click a node to focus its neighborhood &middot; hover for details</div>
      </div>
    </div>
  </div>

  <div id="cy"></div>
  <div id="tooltip"></div>
  <div id="hint">Drag to pan &middot; scroll to zoom &middot; click background to clear focus</div>

<script>
const DATA = ', payload$json, ';
const CFG  = ', cfg, ';

cytoscape.use(window.cytoscapeFcose);

const stylesheet = [
  { selector: "node", style: {
      "width":"data(size)", "height":"data(size)", "label":"data(label)",
      "font-size":9, "color":"#243047", "text-valign":"center", "text-halign":"center",
      "text-margin-y":-1, "text-outline-color":"#ffffff", "text-outline-width":2.2,
      "border-width":1.5, "border-color":"#5a6472", "min-zoomed-font-size":7 } },

  { selector: "node[group=\'miRNA\'][dir=\'up\']",   style:{ "background-color":"var(--x)" } },
  { selector: "node[group=\'miRNA\']", style:{ "shape":"ellipse" } },
  { selector: "node[group=\'miRNA\'][dir=\'up\']",   style:{ "background-color":"#d62728" } },
  { selector: "node[group=\'miRNA\'][dir=\'down\']", style:{ "background-color":"#1f77b4" } },
  { selector: "node[group=\'miRNA\'][dir=\'na\']",   style:{ "background-color":"#888888" } },

  { selector: "node[group=\'Gene\']", style:{ "shape":"ellipse", "background-color":"#4DBBD5" } },
  { selector: "node[group=\'Gene\'][isHub=1]", style:{ "background-color":"#FF7F00",
      "font-size":10, "border-color":"#b35a00" } },

  { selector: "node[group=\'Pathway\']", style:{ "shape":"hexagon", "background-color":"#7B2D8B",
      "color":"#ffffff", "text-outline-color":"#4a1a55", "text-outline-width":2, "font-size":10,
      "font-weight":"bold" } },
  { selector: "node[group=\'Pathway\'][onto=\'CC\']", style:{ "background-color":"#1B9E77", "text-outline-color":"#0c4a38" } },
  { selector: "node[group=\'Pathway\'][onto=\'MF\']", style:{ "background-color":"#D95F02", "text-outline-color":"#7a3500" } },
  { selector: "node[group=\'Pathway\'][onto=\'BP\']", style:{ "background-color":"#6A3D9A", "text-outline-color":"#3a1f59" } },

  { selector: "node[topSig=1]", style:{ "border-width":3.5, "border-color":"#FFD700" } },

  { selector: "edge", style:{ "width":"data(w)", "curve-style":"straight",
      "line-color":"#c2c8cf", "target-arrow-shape":"none", "opacity":0.55 } },
  { selector: "edge[etype=\'mg\'][dir=\'up\']",   style:{ "line-color":"#d62728", "opacity":0.6 } },
  { selector: "edge[etype=\'mg\'][dir=\'down\']", style:{ "line-color":"#1f77b4", "opacity":0.6 } },
  { selector: "edge[etype=\'gp\']", style:{ "line-color":"#b9c0c8", "line-style":"dashed", "width":1, "opacity":0.4 } },

  { selector: ".faded", style:{ "opacity":0.07, "text-opacity":0 } },
  { selector: ".hl-node", style:{ "border-width":4, "border-color":"#111", "z-index":99 } },
  { selector: ".hl-edge", style:{ "opacity":0.95, "width":3.2, "z-index":98 } },
  { selector: ".nolabel node, node.nolabel", style:{ "text-opacity":0 } },
  { selector: ".dimlabel", style:{ "text-opacity":0 } }
];

// Obsidian-style live force controls. Higher repulsion / link distance = more spread;
// higher gravity = tighter clustering toward the centre.
const FORCES = { repulsion: 9000, linkDist: 110, gravity: 0.18 };
function layoutOpts(name){
  if (name === "concentric")
    return { name:"concentric", animate:false, equidistant:false,
             minNodeSpacing: 10 + FORCES.repulsion/400,
             concentric: n => ({ Pathway:3, Gene:2, miRNA:1 })[n.data("group")] || 1,
             levelWidth: () => 1 };
  if (name === "cose")
    return { name:"cose", animate:false, idealEdgeLength: FORCES.linkDist,
             nodeRepulsion: FORCES.repulsion * 100, nodeOverlap: 12,
             gravity: FORCES.gravity, numIter:1500 };
  return { name:"fcose", quality:"proof", animate:false, randomize:true, packComponents:true,
           nodeSeparation: 30 + FORCES.repulsion/120,
           nodeRepulsion: () => FORCES.repulsion,
           idealEdgeLength: e => e.data("etype")==="gp" ? FORCES.linkDist*1.5 : FORCES.linkDist,
           edgeElasticity: 0.4, gravity: FORCES.gravity, gravityRange: 4.5,
           numIter: 3000, tile: true, tilingPaddingVertical: 18, tilingPaddingHorizontal: 18 };
}

const cy = cytoscape({
  container: document.getElementById("cy"),
  elements: [...DATA.nodes, ...DATA.edges],
  style: stylesheet,
  layout: layoutOpts("fcose"),
  wheelSensitivity: 0.22,
  minZoom: 0.08, maxZoom: 4
});

// ── side-panel checkbox lists ────────────────────────────────────────────────
const groupState = { miRNA:true, Gene:true, Pathway:true };
function buildList(group){
  const box = document.getElementById("list-"+group);
  const ns = DATA.nodes.filter(n => n.data.group===group)
                       .sort((a,b)=>a.data.id.localeCompare(b.data.id));
  box.innerHTML = ns.map(n => {
    // When multiple contrasts: add a colored dot in the miRNA list matching the border color
    let dot = "";
    if (group === "miRNA" && CFG.nDea > 1 && n.data.contrastIdx !== undefined) {
      const col = CFG.contrastPalette[n.data.contrastIdx] || "#aaa";
      dot = `<span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:${col};border:1px solid rgba(0,0,0,0.2);margin-right:4px;vertical-align:middle;flex-shrink:0"></span>`;
    }
    return `<label><input type="checkbox" class="chk" data-id="${n.data.id}" checked>${dot}${n.data.label}</label>`;
  }).join("");
}
["miRNA","Gene","Pathway"].forEach(buildList);

// ── filtering (visibility only — no re-layout, edges auto-hide with their nodes) ──
const checkedIds = () => new Set([...document.querySelectorAll(".chk:checked")].map(c=>c.dataset.id));
function applyFilters(){
  const minMir = +document.getElementById("deg-slider").value;
  const checked = checkedIds();
  // Contrast filter (only active when multiple DEA files were loaded)
  const hiddenContrasts = new Set(
    [...document.querySelectorAll(".chk-contrast:not(:checked)")].map(c => c.dataset.contrast)
  );
  cy.batch(() => {
    cy.nodes().forEach(n => {
      const g = n.data("group");
      let vis = checked.has(n.id());
      if (vis && g==="Gene"  && n.data("nmir") < minMir) vis = false;
      if (vis && g==="miRNA" && hiddenContrasts.has(n.data("contrast"))) vis = false;
      n.style("display", vis ? "element" : "none");
    });
    // prune now-orphaned pathways and miRNAs (no visible gene neighbour)
    cy.nodes("[group=\'Pathway\'], [group=\'miRNA\']").forEach(n => {
      if (n.style("display")==="none") return;
      const liveGene = n.connectedEdges().connectedNodes("[group=\'Gene\']")
                        .some(g => g.style("display")!=="none");
      if (!liveGene) n.style("display","none");
    });
  });
  updateStats();
}

// ── live stat chips ──────────────────────────────────────────────────────────
function updateStats(){
  const vis = g => cy.nodes("[group=\'"+g+"\']").filter(n=>n.style("display")!=="none").length;
  const ve  = cy.edges().filter(e=>e.style("display")!=="none"
                && e.source().style("display")!=="none" && e.target().style("display")!=="none").length;
  const chips = [
    ["miRNAs", vis("miRNA"), "#d62728"],
    ["Genes",  vis("Gene"),  "#4DBBD5"],
    [CFG.pathwayLabel, vis("Pathway"), "#7B2D8B"],
    ["Edges",  ve, "#9aa3ad"]
  ];
  document.getElementById("stat-chips").innerHTML = chips.map(c =>
    `<div class="chip"><span class="dot" style="background:${c[2]}"></span>${c[0]} <b>${c[1]}</b></div>`
  ).join("");
}

// ── neighborhood focus on click ──────────────────────────────────────────────
function clearFocus(){ cy.elements().removeClass("faded hl-node hl-edge"); }
cy.on("tap","node", evt => {
  const n = evt.target;
  const hood = n.closedNeighborhood();
  cy.elements().addClass("faded");
  hood.removeClass("faded");
  n.addClass("hl-node");
  n.connectedEdges().removeClass("faded").addClass("hl-edge");
});
cy.on("tap", evt => { if (evt.target === cy) clearFocus(); });

// ── tooltip ──────────────────────────────────────────────────────────────────
const tip = document.getElementById("tooltip");
function showTip(evt){
  const t = evt.target.data("tip"); if (!t) return;
  tip.innerHTML = t; tip.style.display = "block";
  const p = evt.renderedPosition || evt.position;
  const x = Math.min(p.x + 16, window.innerWidth - 300);
  const y = Math.min(p.y + 14, window.innerHeight - 120);
  tip.style.left = x + "px"; tip.style.top = y + "px";
}
cy.on("mouseover","node,edge", showTip);
cy.on("mousemove","node,edge", showTip);
cy.on("mouseout","node,edge", () => tip.style.display="none");
cy.on("pan zoom drag", () => tip.style.display="none");

// ── controls ─────────────────────────────────────────────────────────────────
document.querySelectorAll(".chk").forEach(c => c.addEventListener("change", applyFilters));
document.querySelectorAll("[data-all]").forEach(b => b.addEventListener("click", () => {
  const g = b.dataset.all, on = b.dataset.on==="1";
  document.querySelectorAll("#list-"+g+" .chk").forEach(c => c.checked = on);
  applyFilters();
}));
const degSlider = document.getElementById("deg-slider");
degSlider.addEventListener("input", () => {
  document.getElementById("deg-val").textContent = degSlider.value;
  applyFilters();
});

// Obsidian-style force sliders — re-run the layout (debounced) on the visible subgraph.
let forceTimer = null;
function onForce(){
  FORCES.repulsion = +document.getElementById("f-rep").value;
  FORCES.linkDist  = +document.getElementById("f-lnk").value;
  FORCES.gravity   = +document.getElementById("f-grv").value;
  document.getElementById("rep-val").textContent = FORCES.repulsion;
  document.getElementById("lnk-val").textContent = FORCES.linkDist;
  document.getElementById("grv-val").textContent = FORCES.gravity.toFixed(2);
  clearTimeout(forceTimer);
  forceTimer = setTimeout(() => runLayout(document.getElementById("layout-sel").value), 240);
}
["f-rep","f-lnk","f-grv"].forEach(id => document.getElementById(id).addEventListener("input", onForce));
document.getElementById("btn-fit").addEventListener("click", () => cy.fit(cy.elements(":visible"), 40));
document.getElementById("btn-reset").addEventListener("click", () => {
  clearFocus(); cy.fit(cy.elements(":visible"), 40);
});
document.getElementById("btn-labels").addEventListener("click", () => {
  document.body.classList.toggle("nolabel");
  cy.nodes().toggleClass("dimlabel", document.body.classList.contains("nolabel"));
});
document.getElementById("btn-png").addEventListener("click", () => {
  const png = cy.png({ scale:2, full:true, bg:"#ffffff" });
  const a = document.createElement("a");
  a.href = png; a.download = CFG.pageTitle.replace(/[^a-z0-9]+/gi,"_") + ".png"; a.click();
});
let layoutRunning = false, layoutPending = null;
function runLayout(name){
  if (layoutRunning){ layoutPending = name; return; }
  layoutRunning = true;
  const eles = cy.elements().filter(e => e.style("display")!=="none");
  const lay = eles.layout(layoutOpts(name));
  lay.one("layoutstop", () => {
    layoutRunning = false;
    cy.fit(eles, 40);
    if (layoutPending){ const n = layoutPending; layoutPending = null; runLayout(n); }
  });
  lay.run();
}
document.getElementById("layout-sel").addEventListener("change", e => runLayout(e.target.value));
document.getElementById("btn-relayout").addEventListener("click",
  () => runLayout(document.getElementById("layout-sel").value));

const search = document.getElementById("search");
search.addEventListener("keydown", e => {
  if (e.key !== "Enter") return;
  const q = search.value.trim().toLowerCase(); if (!q) { clearFocus(); return; }
  const hit = cy.nodes().filter(n => n.style("display")!=="none"
                && n.data("label").toLowerCase().includes(q));
  if (hit.length === 0) return;
  clearFocus();
  cy.elements().addClass("faded");
  hit.removeClass("faded").addClass("hl-node");
  hit.connectedEdges().removeClass("faded").addClass("hl-edge");
  hit.connectedEdges().connectedNodes().removeClass("faded");
  cy.animate({ fit:{ eles:hit.closedNeighborhood(), padding:80 } }, { duration:400 });
});

// ── multi-contrast wiring (no-op when nDea === 1) ────────────────────────────
if (CFG.nDea > 1) {
  // Apply per-contrast border color to miRNA nodes via programmatic style
  CFG.contrastNames.forEach((name, i) => {
    const col = CFG.contrastPalette[i] || "#aaa";
    cy.nodes(`[group="miRNA"][contrastIdx=${i}]`).style({ "border-color": col, "border-width": 2.5 });
  });
  // Gold top-sig border overrides contrast border (re-apply so it wins)
  cy.nodes("[group=\'miRNA\'][topSig=1]").style({ "border-color": "#FFD700", "border-width": 3.5 });

  // Wire contrast checkboxes
  document.querySelectorAll(".chk-contrast").forEach(c =>
    c.addEventListener("change", applyFilters));
  document.querySelectorAll("[data-all-contrast]").forEach(b =>
    b.addEventListener("click", () => {
      const on = b.dataset.allContrast === "1";
      document.querySelectorAll(".chk-contrast").forEach(c => c.checked = on);
      applyFilters();
    }));
}

// ── init ─────────────────────────────────────────────────────────────────────
cy.ready(() => { applyFilters(); cy.fit(cy.elements(":visible"), 40); updateStats(); });
</script>
</body>
</html>')

  writeLines(html, out_file)
  log_line(sprintf("Saved: %s  (%.0f KB)", basename(out_file), file.info(out_file)$size/1024))
}

# ── generate dashboards ───────────────────────────────────────────────────────
log_line("Building KEGG dashboard...")
kegg_payload <- build_elements(kegg_df, "KEGG")
write_dashboard(kegg_payload, file.path(OUT_DIR, "Automated_Dashboard_KEGG.html"),
                "miRNA Tripartite Network — KEGG Pathways", "KEGG Pathways", is_go = FALSE)

legacy <- file.path(OUT_DIR, "Automated_Dashboard.html")
if (file.exists(legacy)) file.remove(legacy)

if (!is.null(go_df)) {
  log_line("Building GO dashboard...")
  go_payload <- build_elements(go_df, "GO", onto_col = "ONTOLOGY")
  write_dashboard(go_payload, file.path(OUT_DIR, "Automated_Dashboard_GO.html"),
                  "miRNA Tripartite Network — GO Terms", "GO Terms", is_go = TRUE)
} else {
  log_line("SKIP: GO_Enrichment_Results.csv not found.")
}

log_line("Done.")
