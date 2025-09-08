#!/usr/bin/env Rscript
# ===========================================
# Script: 07_networks_clin/build_network.R
# Descripción: Construye la red miRNA–mRNA a partir de interacciones validadas (y opcionalmente
#              predichas), genera archivos de aristas/nodos, exporta SIF para Cytoscape y
#              calcula métricas de centralidad y módulos (comunidades).
# Fecha: 2025-09-07
# ===========================================

suppressWarnings(suppressMessages({
  library(utils)
  library(tools)
}))

# --- Logging simple (archivo + stdout) ---
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- function(...) {
  msg <- sprintf("[%s] %s", timestamp(), paste0(..., collapse = ""))
  cat(msg, "\n")
  if (!is.null(getOption("net_logfile"))) {
    cat(msg, "\n", file = getOption("net_logfile"), append = TRUE)
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

OUT_DIR  <- file.path(ROOT_DIR, "results", "networks")
FIG_DIR  <- file.path(OUT_DIR, "figures")
LOG_DIR  <- file.path(ROOT_DIR, "outputs", "logs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOG_DIR, "07_build_network.log")
options(net_logfile = LOGFILE)
write("", file = LOGFILE)

# --- Args ---
# uso:
#   Rscript build_network.R [validated_csv] [predicted_csv] [include_predicted] [min_evidence] [min_degree] [mirna_list_file] [string_ppi_file]
# defaults:
#   validated_csv    = results/targets_validated.csv
#   predicted_csv    = results/targets_predicted.csv
#   include_predicted= FALSE
#   min_evidence     = 1   (número mínimo de fuentes/records por arista miRNA–gene)
#   min_degree       = 1   (filtra nodos con grado < min_degree del grafo final)
#   mirna_list_file  = ""  (si se provee, restringe la red a estos miRNA)
#   string_ppi_file  = ""  (opcional TSV/CSV con PPI gene-gene; columnas: geneA,geneB o preferredName_A/B)
args <- commandArgs(trailingOnly = TRUE)
validated_csv     <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else file.path(ROOT_DIR, "results", "targets_validated.csv")
predicted_csv     <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else file.path(ROOT_DIR, "results", "targets_predicted.csv")
include_predicted <- if (length(args) >= 3 && nzchar(args[[3]])) as.logical(args[[3]]) else FALSE
min_evidence      <- if (length(args) >= 4 && nzchar(args[[4]])) as.integer(args[[4]]) else 1
min_degree        <- if (length(args) >= 5 && nzchar(args[[5]])) as.integer(args[[5]]) else 1
mirna_list_file   <- if (length(args) >= 6 && nzchar(args[[6]])) args[[6]] else ""
string_ppi_file   <- if (length(args) >= 7 && nzchar(args[[7]])) args[[7]] else ""

log_line("Starting build_network.R")
log_line(sprintf("validated=%s | predicted=%s | include_predicted=%s | min_evidence=%d | min_degree=%d",
                 validated_csv, predicted_csv, include_predicted, min_evidence, min_degree))
if (nzchar(mirna_list_file)) log_line(sprintf("miRNA filter list: %s", mirna_list_file))
if (nzchar(string_ppi_file)) log_line(sprintf("STRING/other PPI file: %s", string_ppi_file))

# --- Librerías necesarias ---
suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
}))

# --- Helpers ---
read_targets_table <- function(path, kind = c("validated","predicted")) {
  kind <- match.arg(kind)
  if (!file.exists(path)) {
    log_line(sprintf("WARN: %s targets file not found: %s", kind, path))
    return(tibble())
  }
  df <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) {
    log_line(sprintf("WARN: Empty or unreadable %s table: %s", kind, path))
    return(tibble())
  }
  # Normalizar nombres de columnas comunes
  nm <- names(df)
  col_mir <- if ("mirna" %in% nm) "mirna" else if ("mature_mirna_id" %in% nm) "mature_mirna_id" else NA
  col_sym <- if ("target_symbol" %in% nm) "target_symbol" else if ("symbol" %in% nm) "symbol" else NA
  col_db  <- if ("database" %in% nm) "database" else NA
  col_sup <- if ("support_type" %in% nm) "support_type" else if ("evidence" %in% nm) "evidence" else NA
  if (is.na(col_mir) || is.na(col_sym)) {
    log_line(sprintf("WARN: Could not find miRNA/target columns in %s; skipping.", basename(path)))
    return(tibble())
  }
  df2 <- df %>%
    transmute(
      mirna  = .data[[col_mir]] %>% as.character(),
      target = .data[[col_sym]] %>% as.character(),
      database = if (!is.na(col_db)) .data[[col_db]] else NA_character_,
      support  = if (!is.na(col_sup)) .data[[col_sup]] else NA_character_,
      evidence_kind = kind
    ) %>%
    filter(!is.na(mirna), !is.na(target), nchar(mirna) > 0, nchar(target) > 0) %>%
    distinct()
  df2
}

read_list_vec <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(character(0))
  vec <- readLines(path, warn = FALSE)
  vec <- trimws(vec)
  unique(vec[nchar(vec) > 0])
}

read_ppi <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(tibble())
  df <- tryCatch({
    ext <- tolower(file_ext(path))
    if (ext == "csv") readr::read_csv(path, show_col_types = FALSE) else readr::read_tsv(path, show_col_types = FALSE)
  }, error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(tibble())
  cols <- tolower(names(df))
  ga <- NULL; gb <- NULL
  for (cand in list(c("genea","geneb"), c("preferredname_a","preferredname_b"), c("protein1","protein2"))) {
    if (all(cand %in% cols)) { ga <- names(df)[match(cand[1], cols)]; gb <- names(df)[match(cand[2], cols)]; break }
  }
  if (is.null(ga) || is.null(gb)) return(tibble())
  df %>% transmute(geneA = .data[[ga]] %>% as.character(), geneB = .data[[gb]] %>% as.character()) %>%
    filter(!is.na(geneA), !is.na(geneB), nchar(geneA) > 0, nchar(geneB) > 0) %>%
    distinct()
}

# --- Cargar tablas ---
val_tbl  <- read_targets_table(validated_csv, kind = "validated")
pred_tbl <- if (isTRUE(include_predicted)) read_targets_table(predicted_csv, kind = "predicted") else tibble()

if (!nrow(val_tbl) && !nrow(pred_tbl)) stop("No target tables available to build the network.")

all_edges <- bind_rows(val_tbl, pred_tbl)

# Filtrar por lista de miRNAs si se aporta
miR_filter <- read_list_vec(mirna_list_file)
if (length(miR_filter)) {
  before <- nrow(all_edges)
  all_edges <- all_edges %>% filter(mirna %in% miR_filter)
  log_line(sprintf("Applied miRNA filter: kept %d/%d edges", nrow(all_edges), before))
}

# Agregar evidencias por arista (miRNA–gene)
edge_collapsed <- all_edges %>%
  group_by(mirna, target) %>%
  summarize(
    evidence_count = n(),
    sources = paste(sort(unique(na.omit(database))), collapse = ";"),
    kinds   = paste(sort(unique(evidence_kind)), collapse = ";"),
    support_types = paste(sort(unique(na.omit(support))), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(evidence_count), mirna, target)

# Filtrar por mínimo de evidencia
before_e <- nrow(edge_collapsed)
edge_collapsed <- edge_collapsed %>% filter(evidence_count >= min_evidence)
log_line(sprintf("Edges after evidence filter (>= %d): %d/%d", min_evidence, nrow(edge_collapsed), before_e))
if (!nrow(edge_collapsed)) stop("No edges left after applying evidence filter.")

# Construir grafo bipartito miRNA–gene
g_edges <- edge_collapsed %>%
  transmute(from = mirna, to = target, interaction = "targets", evidence_count, sources, kinds, support_types)

# Crear tabla de nodos con tipo
miRNA_nodes <- tibble(node = unique(g_edges$from), type = "miRNA")
gene_nodes  <- tibble(node = unique(g_edges$to),   type = "gene")
nodes_tbl   <- bind_rows(miRNA_nodes, gene_nodes) %>% distinct()

# Grafo inicial
g <- graph_from_data_frame(d = g_edges, vertices = nodes_tbl, directed = FALSE)
V(g)$type <- ifelse(V(g)$name %in% miRNA_nodes$node, "miRNA", "gene")

# Integrar PPI (opcional)
ppi <- read_ppi(string_ppi_file)
if (nrow(ppi)) {
  # Mantener solo genes ya presentes (target nodes)
  ppi2 <- ppi %>% filter(geneA %in% gene_nodes$node, geneB %in% gene_nodes$node) %>% distinct()
  if (nrow(ppi2)) {
    ppi_edges <- ppi2 %>% transmute(from = geneA, to = geneB, interaction = "PPI", evidence_count = NA_integer_,
                                    sources = "PPI", kinds = "PPI", support_types = NA_character_)
    # Unir grafo
    g <- graph_from_data_frame(d = bind_rows(g_edges, ppi_edges), vertices = nodes_tbl, directed = FALSE)
    log_line(sprintf("Added PPI edges: %d", nrow(ppi2)))
  } else {
    log_line("PPI file provided but no overlapping gene-gene edges matched target set.")
  }
}

# Filtrar por grado mínimo (si se solicita)
if (min_degree > 1) {
  deg <- degree(g)
  keep_nodes <- names(deg)[deg >= min_degree]
  g <- induced_subgraph(g, vids = keep_nodes)
  log_line(sprintf("Applied degree filter (>= %d). Nodes kept: %d", min_degree, vcount(g)))
}

# Calcular métricas de centralidad (solo sobre el grafo final)
deg <- degree(g, mode = "all")
btw <- betweenness(g, directed = FALSE, normalized = TRUE)
clo <- closeness(g, normalized = TRUE)
cent_tbl <- tibble(
  node = names(deg),
  type = V(g)$type[match(names(deg), V(g)$name)],
  degree = as.integer(deg),
  betweenness = as.numeric(btw),
  closeness = as.numeric(clo)
) %>% arrange(desc(degree), desc(betweenness))

# Comunidades (Louvain si el grafo es conexo/adecuado)
comm <- tryCatch(cluster_louvain(g), error = function(e) NULL)
if (!is.null(comm)) {
  membership <- membership(comm)
  cent_tbl$module <- membership[match(cent_tbl$node, names(membership))]
  log_line(sprintf("Detected communities (Louvain): %d modules", length(unique(membership))))
} else {
  cent_tbl$module <- NA_integer_
  log_line("WARN: Community detection failed or graph too small.")
}

# --- Salidas ---
edges_csv <- file.path(OUT_DIR, "miRNA_mRNA_edges.csv")
nodes_csv <- file.path(OUT_DIR, "miRNA_mRNA_nodes.csv")
cent_csv  <- file.path(OUT_DIR, "centrality_modules.csv")
sif_path  <- file.path(OUT_DIR, "miRNA_mRNA_network.sif")

# Edges finales desde el grafo (garantiza consistencia con filtros aplicados)
edges_final <- as_data_frame(g, what = "edges") %>%
  rename(source = from, target = to) %>%
  mutate(interaction = ifelse(is.na(interaction), "targets", interaction))
# Completar metadatos de evidencia si faltan (por PPI u otros)
edges_final <- edges_final %>%
  left_join(g_edges %>% rename(source = from, target = to), by = c("source","target","interaction")) %>%
  mutate(
    evidence_count = coalesce(evidence_count, 1L),
    sources        = coalesce(sources, interaction),
    kinds          = coalesce(kinds, interaction),
    support_types  = coalesce(support_types, NA_character_)
  ) %>% distinct()

readr::write_csv(edges_final, edges_csv)
readr::write_csv(nodes_tbl %>% filter(node %in% V(g)$name), nodes_csv)
readr::write_csv(cent_tbl, cent_csv)

# Exportar SIF (source \t interaction \t target)
sif_df <- edges_final %>% transmute(V1 = source, V2 = interaction, V3 = target)
write.table(sif_df, file = sif_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

log_line(sprintf("Saved edges: %s (rows=%d)", edges_csv, nrow(edges_final)))
log_line(sprintf("Saved nodes: %s (rows=%d)", nodes_csv, nrow(nodes_tbl %>% filter(node %in% V(g)$name))))
log_line(sprintf("Saved centrality/modules: %s (rows=%d)", cent_csv, nrow(cent_tbl)))
log_line(sprintf("Saved SIF: %s", sif_path))
log_line(sprintf("Graph summary — nodes: %d, edges: %d", vcount(g), ecount(g)))

# Pequeña figura opcional: distribución de grados (PNG)
suppressWarnings(suppressMessages(library(ggplot2)))
deg_df <- data.frame(degree = degree(g))
png(file.path(FIG_DIR, "degree_distribution.png"), width = 800, height = 600, res = 120)
print(ggplot(deg_df, aes(x = degree)) + geom_histogram(binwidth = 1) +
        labs(title = "Degree distribution", x = "Degree", y = "Count"))
dev.off()
log_line("Saved figure: degree_distribution.png")

# --- Fin ---
log_line("Script finalizado correctamente.")
