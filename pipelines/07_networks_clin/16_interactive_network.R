#!/usr/bin/env Rscript
# ===========================================
# Script: 16_interactive_network.R
# Description: Generates Ultra-Lightweight HTML dashboards with "Auto-Freeze" physics.
#              Each dashboard features 3 multi-select menus that automatically reveal
#              biological connections (miRNA -> Gene -> Pathway/GO term) when selected.
#
#              Two dashboards are produced:
#                Automated_Dashboard_KEGG.html  (miRNA -> Gene -> KEGG pathway)
#                Automated_Dashboard_GO.html    (miRNA -> Gene -> GO term)
#              Automated_Dashboard.html is kept as a copy of the KEGG dashboard
#              for backward compatibility with older report links.
# ===========================================

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(jsonlite)
}))

RUN_DIR <- getwd()
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else RUN_DIR

OUT_DIR <- file.path(RUN_DIR, "network_outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_line <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
log_line("Starting Auto-Freeze Tripartite Dashboard Generation...")

find_file <- function(pattern, dir, required = TRUE) {
  files <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    if (required) stop(sprintf("Could not find file matching: %s", pattern))
    return(NULL)
  }
  return(files[1])
}

dea_df  <- read_csv(find_file("^DEA_.*\\.csv$", input_dir), show_col_types = FALSE)
val_df  <- read_csv(find_file("^targets_validated\\.csv$", input_dir), show_col_types = FALSE)
kegg_df <- read_csv(find_file("^KEGG_Enrichment_Results\\.csv$", input_dir), show_col_types = FALSE)
go_path <- find_file("^GO_Enrichment_Results\\.csv$", input_dir, required = FALSE)
go_df   <- if (!is.null(go_path)) read_csv(go_path, show_col_types = FALSE) else NULL

# --- Robust Top-Down Filtering ---
TOP_PATHWAYS     <- 5
TOP_MIRNAS_LIMIT <- 30

significant_mirnas <- dea_df %>% arrange(desc(abs(S_miRNA))) %>% head(TOP_MIRNAS_LIMIT) %>% pull(miRNA)

# Build a tripartite (miRNA -> Gene -> Term) node/edge JSON payload from any
# clusterProfiler-style enrichment table that has Description / geneID / p.adjust.
build_tripartite_json <- function(enrich_df) {
  term_expanded <- enrich_df %>% arrange(p.adjust) %>% head(TOP_PATHWAYS) %>%
    separate_rows(geneID, sep = "/") %>% select(Pathway = Description, Gene = geneID)

  targets_filtered <- val_df %>%
    filter(target_symbol %in% term_expanded$Gene, mirna %in% significant_mirnas) %>%
    select(miRNA = mirna, Gene = target_symbol)

  valid_genes   <- unique(targets_filtered$Gene)
  term_expanded <- term_expanded %>% filter(Gene %in% valid_genes)

  edges_mirna_gene   <- targets_filtered %>% rename(from = miRNA, to = Gene)
  edges_gene_pathway <- term_expanded %>% rename(from = Gene, to = Pathway)
  master_edges       <- bind_rows(edges_mirna_gene, edges_gene_pathway) %>% distinct()

  if (nrow(master_edges) == 0) return(NULL)

  all_nodes <- unique(c(master_edges$from, master_edges$to))
  node_type <- case_when(
    all_nodes %in% edges_mirna_gene$from ~ "miRNA",
    all_nodes %in% edges_gene_pathway$to ~ "Pathway",
    TRUE ~ "Gene"
  )
  node_size  <- case_when(node_type == "Pathway" ~ 40, node_type == "miRNA" ~ 30, TRUE ~ 15)
  node_color <- case_when(node_type == "Pathway" ~ "#00A087", node_type == "miRNA" ~ "#E64B35", TRUE ~ "#4DBBD5")

  nodes <- data.frame(id = all_nodes, label = all_nodes, group = node_type,
                      value = node_size, color = node_color)
  toJSON(list(nodes = nodes, edges = master_edges), pretty = FALSE)
}

# --- Generate a standalone HTML dashboard from a JSON payload ---
write_dashboard <- function(json_data, out_file, page_title, pathway_label) {
  if (is.null(json_data)) {
    log_line(sprintf("SKIP: no connected nodes for '%s' — dashboard not written.", page_title))
    return(invisible(NULL))
  }
  html_code <- c(
    '<!DOCTYPE html>',
    '<html lang="en">',
    '<head>',
    '    <meta charset="UTF-8">',
    paste0('    <title>', page_title, '</title>'),
    '    <script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>',
    '    <style>',
    '        body { font-family: "Segoe UI", Tahoma, sans-serif; margin: 0; padding: 0; background-color: #f8f9fa; }',
    '        #network-container { width: 100vw; height: 100vh; position: absolute; z-index: 1; }',
    '        #control-panel { position: absolute; top: 10px; left: 10px; z-index: 10; background: rgba(255,255,255,0.95); padding: 15px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.1); width: 300px; max-height: 95vh; display: flex; flex-direction: column; }',
    '        h2 { margin: 0 0 10px 0; font-size: 1.1rem; color: #2c3e50; text-align: center; border-bottom: 2px solid #ecf0f1; padding-bottom: 10px;}',
    '        .filter-group { margin-bottom: 15px; display: flex; flex-direction: column; min-height: 0; }',
    '        .filter-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px; font-size: 14px; }',
    '        .btn-group button { font-size: 11px; padding: 2px 6px; cursor: pointer; border: 1px solid #bdc3c7; background: #fff; border-radius: 3px; }',
    '        .btn-group button:hover { background: #ecf0f1; }',
    '        .checkbox-list { border: 1px solid #ced4da; border-radius: 4px; padding: 5px; overflow-y: auto; height: 120px; background: #fff; font-size: 13px; }',
    '        .checkbox-list label { display: block; margin-bottom: 3px; cursor: pointer; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }',
    '        .checkbox-list input { margin-right: 6px; cursor: pointer; }',
    '    </style>',
    '</head>',
    '<body>',
    '    <div id="control-panel">',
    paste0('        <h2>', page_title, '</h2>'),
    '        <div class="filter-group">',
    '            <div class="filter-header">',
    '               <span style="color:#E64B35; font-weight:bold;">miRNAs</span>',
    '               <div class="btn-group"><button onclick="selectAll(\'chk-mirna\', true)">All</button> <button onclick="selectAll(\'chk-mirna\', false)">None</button></div>',
    '            </div>',
    '            <div class="checkbox-list" id="list-mirna"></div>',
    '        </div>',
    '        <div class="filter-group">',
    '            <div class="filter-header">',
    '               <span style="color:#4DBBD5; font-weight:bold;">Genes</span>',
    '               <div class="btn-group"><button onclick="selectAll(\'chk-gene\', true)">All</button> <button onclick="selectAll(\'chk-gene\', false)">None</button></div>',
    '            </div>',
    '            <div class="checkbox-list" id="list-gene"></div>',
    '        </div>',
    '        <div class="filter-group">',
    '            <div class="filter-header">',
    paste0('               <span style="color:#00A087; font-weight:bold;">', pathway_label, '</span>'),
    '               <div class="btn-group"><button onclick="selectAll(\'chk-pathway\', true)">All</button> <button onclick="selectAll(\'chk-pathway\', false)">None</button></div>',
    '            </div>',
    '            <div class="checkbox-list" id="list-pathway"></div>',
    '        </div>',
    '    </div>',
    '    <div id="network-container"></div>',

    '    <script type="text/javascript">',
    paste0('        const networkData = ', json_data, ';'),
    '        const nodes = new vis.DataSet(networkData.nodes);',
    '        const edges = new vis.DataSet(networkData.edges);',

    '        // Map group types and adjacency for fast biological tracing',
    '        const groupMap = {}; const adj = {};',
    '        networkData.nodes.forEach(n => { groupMap[n.id] = n.group; adj[n.id] = []; });',
    '        networkData.edges.forEach(e => { adj[e.from].push(e.to); adj[e.to].push(e.from); });',

    '        // Populate Checkboxes',
    '        function populateList(groupId, containerId, chkClass) {',
    '            const container = document.getElementById(containerId);',
    '            const groupNodes = networkData.nodes.filter(n => n.group === groupId).sort((a,b) => a.id.localeCompare(b.id));',
    '            groupNodes.forEach(node => {',
    '                const div = document.createElement("div");',
    '                div.innerHTML = `<label><input type="checkbox" class="${chkClass}" value="${node.id}" checked onchange="updateFilters()"> ${node.id}</label>`;',
    '                container.appendChild(div);',
    '            });',
    '        }',
    '        populateList("miRNA", "list-mirna", "chk-mirna");',
    '        populateList("Gene", "list-gene", "chk-gene");',
    '        populateList("Pathway", "list-pathway", "chk-pathway");',

    '        // --- THE BIOLOGICAL TRACER LOGIC ---',
    '        function updateFilters() {',
    '            let explicitSelections = new Set();',
    '            document.querySelectorAll("input[type=checkbox]:checked").forEach(cb => explicitSelections.add(cb.value));',
    '            ',
    '            let visibleSet = new Set(explicitSelections);',
    '            let totalCheckboxes = document.querySelectorAll("input[type=checkbox]").length;',
    '            ',
    '            // If everything is checked, skip traversal and show all',
    '            if (explicitSelections.size !== totalCheckboxes) {',
    '                explicitSelections.forEach(nodeId => {',
    '                    let type = groupMap[nodeId];',
    '                    let neighbors = adj[nodeId] || [];',
    '                    ',
    '                    neighbors.forEach(neighbor => {',
    '                        visibleSet.add(neighbor); // Add direct connection',
    '                        let nType = groupMap[neighbor];',
    '                        ',
    '                        // If jumping from edge to center (miRNA/Pathway -> Gene), grab the opposite side too!',
    '                        if ((type === "miRNA" || type === "Pathway") && nType === "Gene") {',
    '                            (adj[neighbor] || []).forEach(nn => {',
    '                                if (groupMap[nn] !== type) visibleSet.add(nn);',
    '                            });',
    '                        }',
    '                    });',
    '                });',
    '            }',
    '            ',
    '            // Update nodes locally (fast batch update)',
    '            let nodeUpdates = [];',
    '            networkData.nodes.forEach(n => { nodeUpdates.push({id: n.id, hidden: !visibleSet.has(n.id)}); });',
    '            nodes.update(nodeUpdates);',
    '            ',
    '            // Wake up gravity to adapt to the new shape, it will auto-freeze!',
    '            network.setOptions({ physics: true });',
    '            network.stabilize(50); ',
    '        }',

    '        function selectAll(chkClass, isChecked) {',
    '            document.querySelectorAll("." + chkClass).forEach(cb => cb.checked = isChecked);',
    '            updateFilters();',
    '        }',

    '        // --- NETWORK INITIALIZATION ---',
    '        const container = document.getElementById("network-container");',
    '        const options = {',
    '            nodes: { shape: "dot", font: { size: 14, face: "Tahoma", strokeWidth: 3, strokeColor: "#ffffff" }, borderWidth: 2 },',
    '            edges: { width: 1.5, color: { color: "#bdc3c7", highlight: "#34495e" }, arrows: { to: { enabled: true, scaleFactor: 0.5 } }, smooth: false },',
    '            physics: { forceAtlas2Based: { gravitationalConstant: -60, centralGravity: 0.01, springLength: 100, springConstant: 0.08 }, solver: "forceAtlas2Based" },',
    '            interaction: { hover: true, tooltipDelay: 200 }',
    '        };',
    '        const network = new vis.Network(container, { nodes: nodes, edges: edges }, options);',
    '        ',
    '        // THE AUTO-FREEZE MAGIC: Once gravity settles the nodes, turn physics off entirely.',
    '        network.on("stabilizationIterationsDone", function () {',
    '            network.setOptions( { physics: false } );',
    '        });',
    '    </script>',
    '</body>',
    '</html>'
  )
  writeLines(html_code, out_file)
  log_line(sprintf("SUCCESS! Dashboard saved to: %s", out_file))
}

# --- KEGG dashboard ---
log_line("Building KEGG dashboard...")
kegg_json <- build_tripartite_json(kegg_df)
kegg_out  <- file.path(OUT_DIR, "Automated_Dashboard_KEGG.html")
write_dashboard(kegg_json, kegg_out,
                "UC Tripartite Network — KEGG Pathways", "KEGG Pathways")

# Keep legacy filename pointing at the KEGG dashboard for older report links.
if (!is.null(kegg_json)) {
  file.copy(kegg_out, file.path(OUT_DIR, "Automated_Dashboard.html"), overwrite = TRUE)
}

# --- GO dashboard ---
if (!is.null(go_df)) {
  log_line("Building GO dashboard...")
  go_json <- build_tripartite_json(go_df)
  write_dashboard(go_json, file.path(OUT_DIR, "Automated_Dashboard_GO.html"),
                  "UC Tripartite Network — GO Terms", "GO Terms")
} else {
  log_line("SKIP: GO_Enrichment_Results.csv not found — GO dashboard not built.")
}

log_line("Done.")
