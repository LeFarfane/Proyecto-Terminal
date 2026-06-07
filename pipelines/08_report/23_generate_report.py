#!/usr/bin/env python3
"""
Stage 17 — Final DEA HTML Report Generator  (v3)

Accepts a parent results directory that contains:
  figures_and_QC/            PCA, UMAP, heatmap, vst_normalized_matrix.csv
  DEA_Default_Replace/       DEA CSVs + volcanos + multimir/pathway/network outputs
  DEA_Strict_No_Replace/     DEA CSVs + volcanos (same downstream analyses not required)

All DEA sub-folders and comparisons are discovered automatically.
Downstream analyses (KEGG, GO, multiMiR, network) are read from the first
Default sub-folder found.

Usage:
    py -3.11 17_generate_report.py <results_dir>

    Defaults to current directory if no argument is given.

Output:
    <results_dir>/Final_Report.html
"""

import sys
import os
import re
import json
import base64
import io
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Helpers ───────────────────────────────────────────────────────────────────

def _file_b64(path):
    if not path or not os.path.exists(path):
        return None
    ext = os.path.splitext(path)[1].lstrip(".")
    with open(path, "rb") as f:
        return f"data:image/{ext};base64,{base64.b64encode(f.read()).decode()}"


def _fig_b64(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
    buf.seek(0)
    enc = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return f"data:image/png;base64,{enc}"


def _wrap(s, n=56):
    s = str(s)
    return s if len(s) <= n else s[:n] + "…"


def _gene_preview(g, n=12):
    if pd.isna(g):
        return "N/A"
    genes = str(g).split("/")
    total = len(genes)
    return ", ".join(genes[:n]) + (f" … (+{total - n} more)" if total > n else "")


def _parse_ratio(s):
    """Parse clusterProfiler GeneRatio/BgRatio string 'a/b' → float. Returns 0.0 on failure."""
    try:
        a, b = str(s).split("/")
        return int(a) / int(b)
    except Exception:
        return 0.0


def _missing(msg):
    return (f'<p style="color:#c0c0c0;font-size:.88em;padding:22px 0;'
            f'text-align:center">{msg}</p>')


def _lb_img(path, alt, style="width:100%;border-radius:6px"):
    b64 = _file_b64(path)
    if not b64:
        return _missing(f"{alt} — file not found")
    return (f'<img class="lb" src="{b64}" alt="{alt}" '
            f'title="Click to enlarge" style="{style}">')


# ── Recursive locators (outputs may live at the dataset root or nested) ───────

def find_dir(parent, name):
    """Return the shallowest directory called `name` found under `parent`, else None."""
    best = None
    best_depth = None
    for root, dirs, _ in os.walk(parent):
        if name in dirs:
            cand = os.path.join(root, name)
            depth = cand.count(os.sep)
            if best is None or depth < best_depth:
                best, best_depth = cand, depth
    return best


def find_file(parent, basename):
    """Return the shallowest file called `basename` found under `parent`, else None."""
    best = None
    best_depth = None
    for root, _, files in os.walk(parent):
        if basename in files:
            cand = os.path.join(root, basename)
            depth = cand.count(os.sep)
            if best is None or depth < best_depth:
                best, best_depth = cand, depth
    return best


# ── Discovery ─────────────────────────────────────────────────────────────────

def discover_all(parent_dir):
    """
    Walks parent_dir and returns:
      comparisons  list of dicts, one per (DEA CSV × model) combination
      qc           dict with paths/metadata for figures_and_QC
      default_dir  path to the first Default DEA sub-folder (for downstream)
    """
    comparisons = []
    default_dir = None

    # Collect all DEA_* folders recursively
    dea_folders = []
    for root, dirs, _ in os.walk(parent_dir):
        for folder_name in sorted(dirs):
            if folder_name.startswith("DEA_"):
                dea_folders.append((folder_name, os.path.join(root, folder_name)))

    for folder_name, folder_path in dea_folders:

        if "Strict" in folder_name or "No_Outlier" in folder_name:
            model = "Strict"
        elif "Default" in folder_name or "Outlier_Replacement_On" in folder_name:
            model = "Default"
        else:
            model = folder_name

        if model == "Default" and default_dir is None:
            default_dir = folder_path

        # Dispersion plot for this model
        disp_path = None
        for fn in os.listdir(folder_path):
            if fn.startswith("dispersion") and fn.endswith(".png"):
                disp_path = os.path.join(folder_path, fn)
                break

        for fname in sorted(os.listdir(folder_path)):
            if not (fname.startswith("DEA_") and fname.endswith(".csv")):
                continue

            dea_df = pd.read_csv(os.path.join(folder_path, fname))

            # Strip model prefix to get comparison name
            base = fname.replace(".csv", "")
            for prefix in ("DEA_Default_", "DEA_Strict_", "DEA_"):
                base = base.replace(prefix, "")
            comparison = base.replace("_", " ")          # "UC NonResponder"
            key = f"{base}_{model}"                      # "UC_NonResponder_Default"

            # Match volcano image
            volcano_path = os.path.join(
                folder_path,
                fname.replace("DEA_", "volcano_").replace(".csv", ".png"))
            if not os.path.exists(volcano_path):
                volcano_path = None

            sig = dea_df[dea_df["padj"] < 0.05]
            comparisons.append({
                "key":             key,
                "label":           f"{comparison} ({model})",
                "short":           comparison,
                "model":           model,
                "folder_path":     folder_path,
                "dea_df":          dea_df,
                "n_total":         len(dea_df),
                "n_up":            int((sig["log2FoldChange"] > 0).sum()),
                "n_down":          int((sig["log2FoldChange"] < 0).sum()),
                "volcano_path":    volcano_path,
                "dispersion_path": disp_path,
            })

    # Fallback: if no "Default" folder found, use first available dea folder
    if default_dir is None and dea_folders:
        default_dir = dea_folders[0][1]

    # QC folder — search recursively
    qc_dir = None
    for root, dirs, _ in os.walk(parent_dir):
        if "figures_and_QC" in dirs:
            qc_dir = os.path.join(root, "figures_and_QC")
            break
    qc = {"dir": qc_dir}
    if qc["dir"]:
        for fname, attr in (("PCA_plot.png", "pca"),
                            ("UMAP_plot.png", "umap"),
                            ("heatmap_ALL_miRNAs.png", "heatmap")):
            p = os.path.join(qc_dir, fname)
            qc[attr] = p if os.path.exists(p) else None

        vst = os.path.join(qc_dir, "vst_normalized_matrix.csv")
        if os.path.exists(vst):
            qc["vst_path"] = vst
            try:
                tmp = pd.read_csv(vst, index_col=0, nrows=1)
                qc["n_samples"] = len(tmp.columns)
            except Exception:
                qc["n_samples"] = "?"
            try:
                qc["n_mirnas"] = sum(1 for _ in open(vst, encoding="utf-8")) - 1
            except Exception:
                qc["n_mirnas"] = "?"
        else:
            qc["vst_path"] = None

    return comparisons, qc, default_dir


# ── Loaders for downstream analyses ──────────────────────────────────────────

def load_kegg(pathway_dir):
    if not pathway_dir: return None
    p = os.path.join(pathway_dir, "KEGG_Enrichment_Results.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def load_go(pathway_dir):
    if not pathway_dir: return None
    p = os.path.join(pathway_dir, "GO_Enrichment_Results.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def load_multimir(multimir_dir):
    if not multimir_dir: return None
    p = os.path.join(multimir_dir, "targets_validated_summary_by_miRNA.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def count_targets(multimir_dir):
    if not multimir_dir: return 0
    p = os.path.join(multimir_dir, "validated_target_genes_unique.csv")
    try:
        return len(pd.read_csv(p))
    except Exception:
        return 0


def load_clinical(parent_dir):
    """
    Discover clinical_outputs/ produced by clinical_roc.R, hmdd_prep.R,
    and ibd_target_overlap.R.  Returns None if the folder doesn't exist.
    """
    # Files are produced by stages 17/18/19 and may live in clinical_outputs/,
    # ibd_overlap_outputs/, or hmdd_outputs/ — locate each by name recursively.
    file_map = (
        ("roc_R_vs_NR.csv",       "roc_rvsnr"),
        ("roc_UC_vs_HC.csv",      "roc_uchc"),
        ("cv_model_R_vs_NR.csv",  "cv_model"),
        ("candidate_mirnas.csv",  "candidates"),
        ("ibd_target_overlap.csv","ibd_overlap"),
        ("hmdd_mirna_table.csv",  "hmdd_table"),
    )
    found = {key: find_file(parent_dir, fname) for fname, key in file_map}
    if not any(found.values()):
        return None

    result = {"dir": parent_dir}
    for key, p in found.items():
        try:
            result[key] = pd.read_csv(p) if p else None
        except Exception:
            result[key] = None

    # Clinical ROC figures (stage 17) live in a clinical_outputs/figures folder
    clin_dir = find_dir(parent_dir, "clinical_outputs")
    fig_dir = os.path.join(clin_dir, "figures") if clin_dir else None
    result["fig_dir"] = fig_dir if fig_dir and os.path.isdir(fig_dir) else None

    if result["fig_dir"]:
        def _pngs(prefix):
            return sorted(
                os.path.join(fig_dir, f)
                for f in os.listdir(fig_dir)
                if f.startswith(prefix) and f.endswith(".png")
            )
        result["roc_rvsnr_pngs"] = _pngs("ROC_RvsNR_")
        result["roc_uchc_pngs"]  = _pngs("ROC_UCvsHC_")
        result["violin_pngs"]    = _pngs("violin_")
        result["roc_multi"]      = os.path.join(fig_dir, "ROC_multi_RvsNR.png")
    else:
        result.update({"roc_rvsnr_pngs": [], "roc_uchc_pngs": [],
                        "violin_pngs": [], "roc_multi": None})

    # IBD enrichment figure lives in ibd_overlap_outputs/figures/ (stage 19)
    result["ibd_fig"] = find_file(parent_dir, "ibd_enrichment.png")
    return result


# ── Plotly builders ───────────────────────────────────────────────────────────

def volcano_json(dea, label):
    df = dea.copy()
    df["ny"] = -np.log10(df["padj"].clip(lower=1e-320))
    df["grp"] = "n.s."
    df.loc[(df.padj < 0.05) & (df.log2FoldChange > 0), "grp"] = "Up"
    df.loc[(df.padj < 0.05) & (df.log2FoldChange < 0), "grp"] = "Down"

    style = {"n.s.": ("#aaa", 6, .50), "Up": ("#d62728", 9, .88), "Down": ("#1f77b4", 9, .88)}
    traces = []
    for grp, (color, sz, op) in style.items():
        s = df[df.grp == grp]
        traces.append({
            "type": "scatter", "mode": "markers", "name": grp,
            "x": s["log2FoldChange"].tolist(), "y": s["ny"].tolist(),
            "text": s["miRNA"].tolist(),
            "customdata": list(zip(s["baseMean"].round(1), s["padj"].apply(lambda x: f"{x:.2e}"), s["S_miRNA"].round(2))),
            "hovertemplate": ("<b>%{text}</b><br>log₂FC: %{x:.3f}<br>"
                              "−log₁₀(padj): %{y:.2f}<br>baseMean: %{customdata[0]}<br>"
                              "padj: %{customdata[1]}<br>S_miRNA: %{customdata[2]}<extra></extra>"),
            "marker": {"color": color, "size": sz, "opacity": op, "line": {"width": 0}},
        })

    xmin = float(df["log2FoldChange"].min()) - 0.4
    xmax = float(df["log2FoldChange"].max()) + 0.4
    layout = {
        "title": {"text": f"Volcano — {label}", "font": {"size": 14}},
        "xaxis": {"title": "log₂FC", "zeroline": True, "zerolinecolor": "#ccc"},
        "yaxis": {"title": "−log₁₀(padj)"},
        "shapes": [{"type": "line", "x0": xmin, "x1": xmax,
                    "y0": float(-np.log10(0.05)), "y1": float(-np.log10(0.05)),
                    "line": {"dash": "dash", "color": "#888", "width": 1}}],
        "legend": {"orientation": "h", "y": -0.2},
        "hovermode": "closest",
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "height": 440, "margin": {"t": 45, "b": 75, "l": 60, "r": 20},
    }
    return {"data": traces, "layout": layout}


def lfc_json(dea, label, n=12):
    up   = dea.nlargest(n, "S_miRNA")
    down = dea.nsmallest(n, "S_miRNA").iloc[::-1]
    df   = pd.concat([down, up])
    colors = ["#1f77b4" if v < 0 else "#d62728" for v in df["log2FoldChange"]]
    trace = {
        "type": "bar", "orientation": "h",
        "x": df["log2FoldChange"].tolist(), "y": df["miRNA"].tolist(),
        "marker": {"color": colors},
        "customdata": list(zip(df["baseMean"].round(1), df["padj"].apply(lambda x: f"{x:.2e}"), df["S_miRNA"].round(2))),
        "hovertemplate": ("<b>%{y}</b><br>log₂FC: %{x:+.3f}<br>"
                          "baseMean: %{customdata[0]}<br>padj: %{customdata[1]}<br>"
                          "S_miRNA: %{customdata[2]}<extra></extra>"),
    }
    layout = {
        "title": {"text": f"Top {n} miRNAs — {label}", "font": {"size": 13}},
        "xaxis": {"title": "log₂FC", "zeroline": True, "zerolinecolor": "#bbb"},
        "yaxis": {"autorange": "reversed"},
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "height": 460, "margin": {"t": 45, "b": 60, "l": 150, "r": 70},
    }
    return {"data": [trace], "layout": layout}


def overview_json(comps):
    labels = [c["label"] for c in comps]
    traces = [
        {"type": "bar", "name": "Upregulated",
         "x": labels, "y": [c["n_up"] for c in comps],
         "marker": {"color": "#d62728", "opacity": .82},
         "hovertemplate": "<b>%{x}</b><br>Up: %{y}<extra></extra>"},
        {"type": "bar", "name": "Downregulated",
         "x": labels, "y": [-c["n_down"] for c in comps],
         "marker": {"color": "#1f77b4", "opacity": .82},
         "customdata": [c["n_down"] for c in comps],
         "hovertemplate": "<b>%{x}</b><br>Down: %{customdata}<extra></extra>"},
    ]
    layout = {
        "title": {"text": "DE miRNAs per Comparison (padj < 0.05)", "font": {"size": 14}},
        "barmode": "relative",
        "xaxis": {"tickangle": -20},
        "yaxis": {"title": "# miRNAs", "zeroline": True, "zerolinecolor": "#aaa"},
        "legend": {"orientation": "h", "y": -0.3},
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "height": 340, "margin": {"t": 45, "b": 90, "l": 70, "r": 20},
    }
    return json.dumps({"data": traces, "layout": layout})


_CAT_COLORS = {
    "Human Diseases":                       "#e74c3c",
    "Cellular Processes":                   "#3498db",
    "Environmental Information Processing": "#27ae60",
    "Organismal Systems":                   "#9b59b6",
    "Genetic Information Processing":       "#e67e22",
    "Metabolism":                           "#1abc9c",
}
_ONTO_COLORS = {"BP": "#27ae60", "CC": "#3498db", "MF": "#9b59b6"}
_ONTO_NAMES  = {"BP": "Biological Process", "CC": "Cellular Component", "MF": "Molecular Function"}


def kegg_json(kegg_df, n=20):
    # One bar trace per KEGG category (each its own legend entry) so clicking a
    # category in the legend toggles it and double-clicking isolates it.
    df = kegg_df.nsmallest(n, "p.adjust").sort_values("FoldEnrichment", ascending=True)
    cat_col  = "category" if "category" in df.columns else None
    gene_col = "geneID"   if "geneID"   in df.columns else None
    r_col    = "GeneRatio" if "GeneRatio" in df.columns else None

    df = df.copy()
    df["_cat"]   = df[cat_col].astype(str) if cat_col else "—"
    df["_y"]     = [_wrap(str(d)) for d in df["Description"]]
    df["_genes"] = [_gene_preview(g) for g in df[gene_col]] if gene_col else ["—"] * len(df)
    df["_ratio"] = df[r_col].tolist() if r_col else ["—"] * len(df)
    y_order      = df["_y"].tolist()

    hovertemplate = ("<b>%{y}</b><br>Category: %{customdata[4]}<br>"
                     "Fold Enrichment: %{x:.3f}<br>Gene count: %{customdata[0]}<br>"
                     "Gene ratio: %{customdata[2]}<br>FDR: %{customdata[1]}<br><br>"
                     "<i>Genes:</i> %{customdata[3]}<extra></extra>")

    traces = []
    for cat in list(dict.fromkeys(df["_cat"].tolist())):
        sub = df[df["_cat"] == cat]
        traces.append({
            "type": "bar", "orientation": "h", "showlegend": True,
            "name": cat,
            "x": sub["FoldEnrichment"].tolist(),
            "y": sub["_y"].tolist(),
            "marker": {"color": _CAT_COLORS.get(cat, "#7f8c8d"), "opacity": .88},
            "customdata": list(zip(
                sub["Count"].tolist(),
                sub["p.adjust"].apply(lambda x: f"{x:.2e}"),
                sub["_ratio"].tolist(),
                sub["_genes"].tolist(),
                sub["_cat"].tolist(),
            )),
            "hovertemplate": hovertemplate,
        })
    layout = {
        "title": {"text": f"Top {n} KEGG Pathways (FDR < 0.05)", "font": {"size": 14}},
        "barmode": "overlay",
        "xaxis": {"title": "Fold Enrichment"},
        "yaxis": {"autorange": "reversed", "tickfont": {"size": 10},
                  "categoryorder": "array", "categoryarray": y_order},
        "legend": {"title": {"text": "KEGG Category (double-click to isolate)"},
                   "x": 1.02, "y": 1, "font": {"size": 9}, "bgcolor": "rgba(255,255,255,.85)"},
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "height": 600, "margin": {"t": 50, "b": 60, "l": 310, "r": 60},
    }
    return json.dumps({"data": traces, "layout": layout})


def kegg_dotplot_json(kegg_df, n=20):
    """Plotly bubble chart replicating clusterProfiler dotplot style for KEGG.
    x = GeneRatio, y = Description, size = Count, color = p.adjust (red = significant).

    One scatter trace per KEGG category (sharing a single FDR colour scale) so that
    clicking / double-clicking a category in the legend isolates that category. A
    legendgroup + groupclick:'togglegroup' ties each data trace to its legend entry.
    The FDR colour bar sits on the right; the category legend is placed below the
    plot so the two no longer overlap.
    """
    df = kegg_df.nsmallest(n, "p.adjust").copy()
    df = df.sort_values("p.adjust", ascending=False)   # most significant ends up at top via autorange:reversed

    cat_col  = "category" if "category" in df.columns else None
    gene_col = "geneID"   if "geneID"   in df.columns else None

    df["_x"]     = [_parse_ratio(r) for r in df["GeneRatio"]]
    df["_y"]     = [_wrap(str(d)) for d in df["Description"]]
    df["_genes"] = [_gene_preview(g) for g in df[gene_col]] if gene_col else ["—"] * len(df)
    df["_cat"]   = df[cat_col].astype(str) if cat_col else "—"

    padj_all   = df["p.adjust"].tolist()
    cmin, cmax = (min(padj_all), max(padj_all)) if padj_all else (0.0, 1.0)
    max_c      = df["Count"].max() if len(df) else 1

    def _msize(c):
        return max(7, round(30 * (c / max_c) ** 0.5))

    hovertemplate = (
        "<b>%{y}</b><br>Category: %{customdata[4]}<br>"
        "GeneRatio: %{customdata[2]} (%{x:.4f})<br>"
        "Gene count: %{customdata[0]}<br>FDR: %{customdata[1]}<br><br>"
        "<i>Genes:</i> %{customdata[3]}<extra></extra>"
    )

    cats_order = list(dict.fromkeys(df["_cat"].tolist()))
    traces, first = [], True
    for cat in cats_order:
        sub    = df[df["_cat"] == cat]
        marker = {
            "size": [_msize(c) for c in sub["Count"]], "sizemode": "diameter",
            "color": sub["p.adjust"].tolist(),
            "colorscale": "RdBu", "reversescale": True,
            "cmin": cmin, "cmax": cmax,
            "line": {"width": 1.5, "color": _CAT_COLORS.get(cat, "#7f8c8d")},
            "opacity": 0.88,
        }
        if first:
            marker["showscale"] = True
            marker["colorbar"]  = {"title": "FDR", "thickness": 14, "len": 0.55,
                                   "y": 0.5, "x": 1.01, "tickformat": ".1e"}
            first = False
        else:
            marker["showscale"] = False
        traces.append({
            "type": "scatter", "mode": "markers",
            "name": cat, "legendgroup": cat, "showlegend": False,
            "x": sub["_x"].tolist(), "y": sub["_y"].tolist(),
            "marker": marker,
            "customdata": list(zip(
                sub["Count"].tolist(),
                [f"{p:.2e}" for p in sub["p.adjust"]],
                sub["GeneRatio"].tolist(),
                sub["_genes"].tolist(),
                sub["_cat"].tolist(),
            )),
            "hovertemplate": hovertemplate,
        })

    # Ghost legend traces — carry the visible category swatch and share the
    # legendgroup so a legend click toggles the whole category group.
    for cat in cats_order:
        traces.append({
            "type": "scatter", "mode": "markers", "x": [None], "y": [None],
            "name": cat, "legendgroup": cat, "showlegend": True, "hoverinfo": "none",
            "marker": {"color": _CAT_COLORS.get(cat, "#7f8c8d"), "size": 12, "symbol": "circle"},
        })

    layout = {
        "title": {"text": f"Top {n} KEGG Pathways — Dotplot (FDR < 0.05)", "font": {"size": 14}},
        "xaxis": {"title": "GeneRatio", "zeroline": False, "tickformat": ".3f"},
        "yaxis": {"autorange": "reversed", "tickfont": {"size": 10}},
        "legend": {"title": {"text": "KEGG Category (click to toggle, double-click to isolate)"},
                   "orientation": "h", "x": 0, "xanchor": "left", "y": -0.16, "yanchor": "top",
                   "font": {"size": 9}, "bgcolor": "rgba(255,255,255,.85)",
                   "groupclick": "togglegroup"},
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "height": 640, "margin": {"t": 50, "b": 130, "l": 320, "r": 90},
    }
    return json.dumps({"data": traces, "layout": layout})


def go_json(go_df, n=15):
    gene_col = "geneID"    if "geneID"    in go_df.columns else None
    r_col    = "GeneRatio" if "GeneRatio" in go_df.columns else None
    result = {}
    for onto, color in _ONTO_COLORS.items():
        sub = go_df[go_df["ONTOLOGY"] == onto].nsmallest(n, "p.adjust").copy()
        if sub.empty:
            continue
        sub = sub.sort_values("FoldEnrichment", ascending=True)
        trace = {
            "type": "bar", "orientation": "h",
            "x": sub["FoldEnrichment"].tolist(),
            "y": [_wrap(str(d)) for d in sub["Description"]],
            "marker": {"color": color, "opacity": .82},
            "customdata": list(zip(
                sub["Count"].tolist(),
                sub["p.adjust"].apply(lambda x: f"{x:.2e}"),
                sub[r_col].tolist() if r_col else ["—"] * len(sub),
                [_gene_preview(g) for g in sub[gene_col]] if gene_col else ["—"] * len(sub),
            )),
            "hovertemplate": (f"<b>%{{y}}</b><br>Ontology: {_ONTO_NAMES[onto]}<br>"
                              "Fold Enrichment: %{x:.3f}<br>Gene count: %{customdata[0]}<br>"
                              "Gene ratio: %{customdata[2]}<br>FDR: %{customdata[1]}<br><br>"
                              "<i>Genes:</i> %{customdata[3]}<extra></extra>"),
        }
        layout = {
            "title": {"text": f"Top {n} GO — {_ONTO_NAMES[onto]} (FDR < 0.05)", "font": {"size": 14}},
            "xaxis": {"title": "Fold Enrichment"},
            "yaxis": {"autorange": "reversed", "tickfont": {"size": 10}},
            "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
            "height": 550, "margin": {"t": 50, "b": 60, "l": 320, "r": 40},
        }
        result[onto] = json.dumps({"data": [trace], "layout": layout})
    return result


def go_dotplot_json(go_df, n=15):
    """Per-ontology Plotly bubble dotplot for GO. Returns {onto: json_string}."""
    gene_col = "geneID" if "geneID" in go_df.columns else None
    result = {}
    for onto, base_color in _ONTO_COLORS.items():
        sub = go_df[go_df["ONTOLOGY"] == onto].nsmallest(n, "p.adjust").copy()
        if sub.empty:
            continue
        sub = sub.sort_values("p.adjust", ascending=False)

        x_vals       = [_parse_ratio(r) for r in sub["GeneRatio"]]
        y_labels     = [_wrap(str(d)) for d in sub["Description"]]
        sizes        = sub["Count"].tolist()
        padj         = sub["p.adjust"].tolist()
        genes        = [_gene_preview(g) for g in sub[gene_col]] if gene_col else ["—"] * len(sub)
        ratios       = sub["GeneRatio"].tolist()
        max_c        = max(sizes) if sizes else 1
        marker_sizes = [max(7, round(30 * (c / max_c) ** 0.5)) for c in sizes]

        trace = {
            "type": "scatter", "mode": "markers", "showlegend": False,
            "x": x_vals, "y": y_labels,
            "marker": {
                "size": marker_sizes, "sizemode": "diameter",
                "color": padj,
                "colorscale": "RdBu", "reversescale": True,
                "colorbar": {"title": "FDR", "thickness": 14, "len": 0.6, "tickformat": ".1e"},
                "line": {"width": 1.5, "color": base_color},   # ontology-tinted border
                "opacity": 0.88,
            },
            "customdata": list(zip(sizes, [f"{p:.2e}" for p in padj], ratios, genes)),
            "hovertemplate": (
                f"<b>%{{y}}</b><br>Ontology: {_ONTO_NAMES[onto]}<br>"
                "GeneRatio: %{customdata[2]} (%{x:.4f})<br>"
                "Gene count: %{customdata[0]}<br>FDR: %{customdata[1]}<br><br>"
                "<i>Genes:</i> %{customdata[3]}<extra></extra>"
            ),
        }
        layout = {
            "title": {"text": f"Top {n} GO — {_ONTO_NAMES[onto]} — Dotplot (FDR < 0.05)", "font": {"size": 14}},
            "xaxis": {"title": "GeneRatio", "zeroline": False, "tickformat": ".3f"},
            "yaxis": {"autorange": "reversed", "tickfont": {"size": 10}},
            "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
            "height": 560, "margin": {"t": 50, "b": 70, "l": 330, "r": 140},
        }
        result[onto] = json.dumps({"data": [trace], "layout": layout})
    return result


def multimir_json(summary_df, n=20):
    df = summary_df.nlargest(n, "validated_targets").sort_values("validated_targets", ascending=True)
    trace = {
        "type": "bar", "orientation": "h",
        "x": df["validated_targets"].tolist(), "y": df["mirna"].tolist(),
        "marker": {"color": "#2ca02c", "opacity": .82},
        "text": [f"{v:,}" for v in df["validated_targets"]], "textposition": "outside",
        "hovertemplate": "<b>%{y}</b><br>Validated targets: %{x:,}<extra></extra>",
    }
    layout = {
        "title": {"text": f"Top {n} miRNAs by Validated Target Count (multiMiR)", "font": {"size": 14}},
        "xaxis": {"title": "Validated Targets"},
        "yaxis": {"autorange": "reversed"},
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "height": 520, "margin": {"t": 50, "b": 60, "l": 175, "r": 90},
    }
    return json.dumps({"data": [trace], "layout": layout})


# ── Clinical Plotly builders ──────────────────────────────────────────────────

def roc_bar_json(roc_df, label, n=20, auc_thresh=0.70):
    """Horizontal bar chart of per-miRNA AUC values."""
    if roc_df is None or roc_df.empty:
        return '{"data":[],"layout":{}}'
    df = roc_df.head(n).sort_values("AUC", ascending=True).copy()
    colors  = ["#d62728" if a >= auc_thresh else "#aec7e8" for a in df["AUC"]]
    n_pos   = df["n_pos"].tolist() if "n_pos" in df.columns else [0] * len(df)
    n_neg   = df["n_neg"].tolist() if "n_neg" in df.columns else [0] * len(df)
    trace   = {
        "type": "bar", "orientation": "h",
        "x": df["AUC"].round(3).tolist(), "y": df["miRNA"].tolist(),
        "marker": {"color": colors, "opacity": 0.85},
        "text": [f"{a:.3f}" for a in df["AUC"]], "textposition": "outside",
        "customdata": list(zip(n_pos, n_neg)),
        "hovertemplate": ("<b>%{y}</b><br>AUC: %{x:.3f}<br>"
                          "n(pos): %{customdata[0]}  n(neg): %{customdata[1]}"
                          "<extra></extra>"),
    }
    layout  = {
        "title": {"text": f"ROC-AUC per miRNA — {label}", "font": {"size": 14}},
        "xaxis": {"title": "AUC", "range": [0.3, 1.12]},
        "yaxis": {"autorange": "reversed"},
        "shapes": [
            {"type": "line", "x0": 0.5, "x1": 0.5,
             "y0": -0.5, "y1": len(df) - 0.5,
             "line": {"dash": "dash", "color": "#888", "width": 1}},
            {"type": "line", "x0": auc_thresh, "x1": auc_thresh,
             "y0": -0.5, "y1": len(df) - 0.5,
             "line": {"dash": "dot", "color": "#e67e22", "width": 1.5}},
        ],
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "height": max(380, 22 * len(df) + 120),
        "margin": {"t": 50, "b": 80, "l": 185, "r": 85},
    }
    return json.dumps({"data": [trace], "layout": layout})


def ibd_bar_json(ibd_df, n=20):
    """Horizontal bar chart of per-miRNA IBD-gene overlap counts."""
    if ibd_df is None or ibd_df.empty:
        return '{"data":[],"layout":{}}'
    df = ibd_df[ibd_df["n_ibd_overlap"] > 0].copy()
    if df.empty:
        return '{"data":[],"layout":{}}'
    df = df.head(n).sort_values("n_ibd_overlap", ascending=True)
    padj_col = "padj" if "padj" in df.columns else None
    colors   = (["#d62728" if p <= 0.05 else "#aec7e8" for p in df[padj_col]]
                if padj_col else ["#aec7e8"] * len(df))
    or_vals  = [round(float(v), 2) if not pd.isna(v) else 0
                for v in (df["OR"] if "OR" in df.columns else [0] * len(df))]
    padj_fmt = ([f"{p:.2e}" for p in df[padj_col]]
                if padj_col else ["—"] * len(df))
    trace    = {
        "type": "bar", "orientation": "h",
        "x": df["n_ibd_overlap"].tolist(), "y": df["miRNA"].tolist(),
        "marker": {"color": colors, "opacity": 0.85},
        "text": [f"OR={o:.2f}" for o in or_vals], "textposition": "outside",
        "customdata": list(zip(or_vals, padj_fmt)),
        "hovertemplate": ("<b>%{y}</b><br>IBD-overlapping targets: %{x}<br>"
                          "OR: %{customdata[0]:.2f}<br>FDR: %{customdata[1]}"
                          "<extra></extra>"),
    }
    layout   = {
        "title": {"text": "Validated Targets Overlapping IBD Gene Set", "font": {"size": 14}},
        "xaxis": {"title": "# IBD-associated validated targets"},
        "yaxis": {"autorange": "reversed"},
        "plot_bgcolor": "#fafafa", "paper_bgcolor": "white",
        "height": max(380, 22 * len(df) + 120),
        "margin": {"t": 50, "b": 60, "l": 185, "r": 105},
    }
    return json.dumps({"data": [trace], "layout": layout})


def cand_table_html(cand_df, score_by_source=None):
    """HTML table of candidate miRNAs from candidate_mirnas.csv.

    score_by_source: {source_token -> {miRNA -> S_miRNA}} lookup built from the
    Default DEA tables, used to show the miRNA score (S_miRNA) for each candidate
    from its source comparison.
    """
    if cand_df is None or cand_df.empty:
        return _missing("candidate_mirnas.csv not found — run clinical_roc.R first")
    score_by_source = score_by_source or {}
    rows = ""
    for _, r in cand_df.iterrows():
        direction = str(r.get("direction", ""))
        cls   = "badge-up"   if direction == "UP"   else "badge-down"
        badge = "▲ Up"       if direction == "UP"   else "▼ Down"
        lfc   = r.get("log2FoldChange", 0)
        padj  = r.get("padj",           1)
        src   = r.get("source", "—")
        token = str(src).split()[-1].lower() if not pd.isna(src) else ""
        score = score_by_source.get(token, {}).get(r["miRNA"])
        score_str = f"{float(score):.2f}" if score is not None and not pd.isna(score) else "—"
        rows += (f'<tr>'
                 f'<td><strong>{r["miRNA"]}</strong></td>'
                 f'<td class="num">{float(lfc):+.3f}</td>'
                 f'<td class="num">{float(padj):.2e}</td>'
                 f'<td class="num">{score_str}</td>'
                 f'<td><span class="{cls}">{badge}</span></td>'
                 f'<td>{src}</td>'
                 f'</tr>\n')
    return (
        '<table class="sortable-tbl">'
        '<thead><tr>'
        '<th onclick="sortTbl(this,0)">miRNA</th>'
        '<th class="num" onclick="sortTbl(this,1)">log₂FC</th>'
        '<th class="num" onclick="sortTbl(this,2)">padj</th>'
        '<th class="num" onclick="sortTbl(this,3)">miRNA score</th>'
        '<th>Direction</th>'
        '<th>Source comparison</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody>'
        '</table>'
    )


def ibd_table_html(ibd_df):
    """HTML table of per-miRNA IBD overlap results."""
    if ibd_df is None or ibd_df.empty:
        return _missing("ibd_target_overlap.csv not found — run ibd_target_overlap.R first")
    rows = ""
    for _, r in ibd_df[ibd_df["n_ibd_overlap"] > 0].head(20).iterrows():
        padj_v = float(r.get("padj", 1) or 1)
        sig    = " ★" if padj_v <= 0.05 else ""
        or_v   = f'{float(r["OR"]):.2f}' if "OR" in r and not pd.isna(r["OR"]) else "—"
        genes  = str(r.get("overlap_genes", "")).replace(";", ", ")
        genes  = genes[:100] + "…" if len(genes) > 100 else genes
        rows  += (f'<tr>'
                  f'<td><strong>{r["miRNA"]}</strong>{sig}</td>'
                  f'<td class="num">{int(r["n_validated"]):,}</td>'
                  f'<td class="num" style="color:#d62728;font-weight:600">{int(r["n_ibd_overlap"])}</td>'
                  f'<td class="num">{float(r.get("pct_ibd", 0)):.1f}%</td>'
                  f'<td class="num">{or_v}</td>'
                  f'<td class="num">{padj_v:.2e}</td>'
                  f'<td style="font-size:.8em;color:#555">{genes}</td>'
                  f'</tr>\n')
    return (
        '<table class="sortable-tbl">'
        '<thead><tr>'
        '<th onclick="sortTbl(this,0)">miRNA (★=FDR≤0.05)</th>'
        '<th class="num" onclick="sortTbl(this,1)">Validated targets</th>'
        '<th class="num" onclick="sortTbl(this,2)">IBD overlap</th>'
        '<th class="num" onclick="sortTbl(this,3)">% IBD</th>'
        '<th class="num" onclick="sortTbl(this,4)">OR</th>'
        '<th class="num" onclick="sortTbl(this,5)">FDR</th>'
        '<th onclick="sortTbl(this,6)">Overlapping IBD genes</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody>'
        '</table>'
    )


# ── DEA table ─────────────────────────────────────────────────────────────────

def dea_table_html(dea, n=40):
    sig = dea[dea["padj"] < 0.05].copy()
    sig = sig.reindex(sig["S_miRNA"].abs().sort_values(ascending=False).index).head(n)
    rows = ""
    for _, r in sig.iterrows():
        d   = "▲ Up" if r["log2FoldChange"] > 0 else "▼ Down"
        cls = "badge-up" if r["log2FoldChange"] > 0 else "badge-down"
        bg  = "#fff7f7" if r["log2FoldChange"] > 0 else "#f7f7ff"
        rows += (f'<tr style="background:{bg}">'
                 f'<td><strong>{r["miRNA"]}</strong></td>'
                 f'<td class="num">{r["baseMean"]:.1f}</td>'
                 f'<td class="num">{r["log2FoldChange"]:+.3f}</td>'
                 f'<td class="num">{r["padj"]:.2e}</td>'
                 f'<td class="num">{r["S_miRNA"]:.2f}</td>'
                 f'<td><span class="{cls}">{d}</span></td></tr>\n')
    return (
        '<table class="sortable-tbl">'
        '<thead><tr>'
        '<th onclick="sortTbl(this,0)">miRNA</th>'
        '<th class="num" onclick="sortTbl(this,1)">Base Mean</th>'
        '<th class="num" onclick="sortTbl(this,2)">log₂FC</th>'
        '<th class="num" onclick="sortTbl(this,3)">padj</th>'
        '<th class="num" onclick="sortTbl(this,4)">S_miRNA</th>'
        '<th>Direction</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody>'
        '</table>'
    )


# ── HTML template ─────────────────────────────────────────────────────────────

_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>miRNA DEA Report — <<<DATASET>>></title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'Segoe UI',Arial,sans-serif;background:#f4f6f9;color:#333;line-height:1.55}

.hdr{background:linear-gradient(135deg,#1a3c5e,#2e7bb5);color:#fff;padding:44px 30px;text-align:center}
.hdr h1{font-size:1.95em;letter-spacing:.4px}
.hdr .sub{margin-top:9px;opacity:.88;font-size:1.05em}
.hdr .sub2{margin-top:5px;opacity:.64;font-size:.88em}

.wrap{max-width:1200px;margin:0 auto;padding:26px 16px}

/* stat cards */
.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:22px 0}
.card{background:#fff;border-radius:10px;padding:22px 16px;text-align:center;
      box-shadow:0 2px 10px rgba(0,0,0,.07);border-top:4px solid #2e7bb5}
.card .val{font-size:2.2em;font-weight:700;color:#2e7bb5}
.card .lbl{color:#888;margin-top:5px;font-size:.86em}

/* sections */
.sec{background:#fff;border-radius:10px;padding:28px;margin:18px 0;
     box-shadow:0 2px 10px rgba(0,0,0,.07)}
.sec h2{font-size:1.3em;color:#1a3c5e;border-bottom:2px solid #e4eaf0;
        padding-bottom:10px;margin-bottom:20px}
.sec h3{font-size:1em;color:#555;margin-bottom:10px}

/* two-column grid */
.two-col{display:grid;grid-template-columns:1fr 1fr;gap:22px}

/* note box */
.note{background:#f8fafc;border-left:4px solid #2e7bb5;padding:11px 15px;
      border-radius:0 7px 7px 0;color:#555;font-size:.9em;margin-bottom:16px}
.note.warn{border-color:#e67e22;background:#fff9f3}

/* tabs — comparison selector */
.tabs{display:flex;gap:8px;margin-bottom:18px;flex-wrap:wrap}
.tab{padding:8px 18px;border:2px solid #2e7bb5;border-radius:22px;background:#fff;
     color:#2e7bb5;cursor:pointer;font-size:.9em;font-weight:600;
     transition:all .18s;white-space:nowrap}
.tab:hover{background:#eef5fd}
.tab.active{background:#2e7bb5;color:#fff}
.tab.strict{border-color:#e67e22;color:#e67e22}
.tab.strict.active{background:#e67e22;color:#fff}

/* DEA table panels */
.dea-tbl{display:none}
.dea-tbl.visible{display:block}

/* tables */
table{width:100%;border-collapse:collapse;font-size:.88em}
th{background:#2e7bb5;color:#fff;padding:9px 12px;text-align:left;
   font-weight:600;cursor:pointer;user-select:none;white-space:nowrap}
th:hover{background:#1a3c5e}
td{padding:7px 12px;border-bottom:1px solid #eef0f3}
td.num{font-variant-numeric:tabular-nums}
tr:hover td{background:#eef5fd!important}
.badge-up{background:#fde8e8;color:#c0392b;padding:2px 9px;border-radius:4px;
          font-size:.82em;font-weight:600}
.badge-down{background:#e8eeff;color:#2143a0;padding:2px 9px;border-radius:4px;
            font-size:.82em;font-weight:600}

/* QC grid */
.qc-grid{display:grid;grid-template-columns:1fr 1fr 1fr;gap:16px}
.qc-card{background:#f8fafc;border-radius:8px;padding:14px;text-align:center}
.qc-card h3{margin-bottom:10px;font-size:.95em;color:#444}

/* buttons / links */
.btn{display:inline-block;background:#2e7bb5;color:#fff;padding:10px 22px;
     border-radius:7px;text-decoration:none;margin:6px 4px;font-size:.93em;
     transition:background .2s}
.btn:hover{background:#1a3c5e}

/* overview stats table */
.ov-tbl td{vertical-align:middle}

/* lightbox */
img.lb{cursor:zoom-in;border-radius:6px;display:block}
#lb-ov{display:none;position:fixed;inset:0;background:rgba(255,255,255,.96);
        z-index:9999;align-items:center;justify-content:center;flex-direction:column}
#lb-ov.open{display:flex}
#lb-wrap{position:relative;display:flex;flex-direction:column;
         align-items:center;gap:12px;max-width:95vw}
#lb-img{max-width:90vw;max-height:82vh;object-fit:contain;transform-origin:center;
        transition:transform .14s;cursor:grab;border-radius:4px;
        box-shadow:0 4px 24px rgba(0,0,0,.18)}
#lb-img.gr{cursor:grabbing}
#lb-ctrls{display:flex;gap:10px}
.lbtn{background:#fff;border:2px solid #d0d8e4;border-radius:50%;width:38px;height:38px;
      font-size:18px;cursor:pointer;display:flex;align-items:center;
      justify-content:center;box-shadow:0 2px 6px rgba(0,0,0,.1);color:#333}
.lbtn:hover{background:#eef5fd;border-color:#2e7bb5}
#lb-x{position:absolute;top:-44px;right:0;color:#555;font-size:30px;
      cursor:pointer;line-height:1;padding:0 6px;opacity:.7}
#lb-x:hover{opacity:1;color:#c0392b}
#lb-hint{color:rgba(0,0,0,.38);font-size:.8em}

.footer{text-align:center;color:#bbb;padding:28px;font-size:.87em}
@media(max-width:700px){
  .cards{grid-template-columns:1fr 1fr}
  .two-col,.qc-grid{grid-template-columns:1fr}
}
</style>
</head>
<body>

<!-- HEADER -->
<div class="hdr">
  <h1>miRNA Differential Expression Analysis</h1>
  <p class="sub">Dataset: <strong><<<DATASET>>></strong> &ensp;&mdash;&ensp; IBD miRNA-seq</p>
  <p class="sub2">Generated <<<DATE>>></p>
</div>

<div class="wrap">

<!-- OVERVIEW STATS CARDS -->
<div class="cards">
  <div class="card">
    <div class="val"><<<N_TOTAL>>></div>
    <div class="lbl">miRNAs tested</div>
  </div>
  <div class="card">
    <div class="val"><<<N_COMPS>>></div>
    <div class="lbl">Comparisons analysed</div>
  </div>
  <div class="card">
    <div class="val" style="color:#27ae60"><<<N_SAMPLES>>></div>
    <div class="lbl">Samples (VST matrix)</div>
  </div>
  <div class="card">
    <div class="val" style="color:#27ae60"><<<N_TARGETS>>></div>
    <div class="lbl">Unique validated target genes</div>
  </div>
</div>

<!-- DE OVERVIEW CHART -->
<div class="sec">
  <h2>Overview — Differentially Expressed miRNAs</h2>
  <div id="plt-overview"></div>
  <br>
  <table class="ov-tbl">
    <thead><tr>
      <th>Comparison</th><th>Model</th>
      <th class="num">Total tested</th>
      <th class="num" style="color:#fde8e8">Upregulated</th>
      <th class="num" style="color:#e8eeff">Downregulated</th>
    </tr></thead>
    <tbody><<<STATS_ROWS>>></tbody>
  </table>
</div>

<!-- QC SECTION -->
<<<QC_SECTION>>>

<!-- DEA EXPLORER -->
<div class="sec">
  <h2>Differential Expression Explorer</h2>
  <div class="note">
    Click a tab to switch comparison. Default model uses Cook&apos;s distance
    outlier replacement; Strict model disables it (minReplicatesForReplace = Inf).
    Hover any point or bar for full statistics.
  </div>
  <div class="tabs" id="comp-tabs"><<<COMP_TABS>>></div>
  <div class="two-col">
    <div>
      <h3>Volcano Plot</h3>
      <div id="plt-volcano"></div>
    </div>
    <div>
      <h3>Top miRNAs by S_miRNA Score</h3>
      <div id="plt-lfc"></div>
    </div>
  </div>
</div>

<!-- DEA TABLES -->
<div class="sec">
  <h2>Significant miRNAs Table — Top 40 (padj &lt; 0.05)</h2>
  <div class="note">
    Sorted by |S_miRNA|. Click any column header to re-sort.
    Use the tabs above to switch comparison.
  </div>
  <div class="tabs" id="tbl-tabs"><<<TBL_TABS>>></div>
  <<<COMP_TABLES>>>
</div>

<!-- DISPERSION PLOTS -->
<div class="sec">
  <h2>DESeq2 Dispersion Estimates</h2>
  <div class="two-col">
    <div>
      <h3>Default model</h3>
      <<<DISP_DEFAULT>>>
    </div>
    <div>
      <h3>Strict model</h3>
      <<<DISP_STRICT>>>
    </div>
  </div>
</div>

<!-- KEGG -->
<div class="sec">
  <h2>KEGG Pathway Enrichment <span style="font-weight:400;font-size:.85em">(Default model — UC comparisons combined)</span></h2>
  <div class="note">Coloured by KEGG category. Hover a bar to see top genes.</div>
  <<<KEGG_CONTENT>>>
</div>

<!-- GO -->
<div class="sec">
  <h2>Gene Ontology Enrichment <span style="font-weight:400;font-size:.85em">(Default model)</span></h2>
  <div class="note">Switch ontology with the tabs. Hover a bar to see top genes.</div>
  <<<GO_CONTENT>>>
</div>

<!-- MULTIMIR -->
<div class="sec">
  <h2>Target Prediction — multiMiR Validated Targets</h2>
  <div class="note">Experimental evidence from miRTarBase, TarBase and miRecords.</div>
  <<<MULTIMIR_CONTENT>>>
</div>

<!-- NETWORK -->
<div class="sec">
  <h2>Network Visualization <span style="font-weight:400;font-size:.85em">(Default model)</span></h2>
  <div class="note">
    Interactive tripartite dashboards (miRNA &rarr; target gene &rarr; pathway / GO term).
    Use the three filter menus on the left to select nodes; connected partners are
    revealed automatically and the layout auto-freezes once it settles. Each dashboard
    opens in a new tab.
  </div>
  <<<NETWORK_LINKS>>>
  <<<TRIPARTITE_IMG>>>
</div>

<<<CLINICAL_SECTION>>>

</div><!-- /wrap -->

<div class="footer">
  miRNA-seq Pipeline &nbsp;&middot;&nbsp; IBD Project &nbsp;&middot;&nbsp; <<<DATASET>>><br>
  DESeq2 + apeglm + multiMiR + clusterProfiler + vis.js
</div>

<!-- LIGHTBOX -->
<div id="lb-ov">
  <div id="lb-wrap">
    <span id="lb-x" title="Close (Esc)">&times;</span>
    <img id="lb-img" src="" alt="">
    <div id="lb-ctrls">
      <button class="lbtn" id="lbzi">+</button>
      <button class="lbtn" id="lbzr">&#8635;</button>
      <button class="lbtn" id="lbzo">&#8722;</button>
    </div>
    <p id="lb-hint">Scroll &bull; Drag to pan &bull; Esc to close &bull; 0 to reset</p>
  </div>
</div>

<script>
/* ===== DATA ===== */
const ALL_COMPS   = <<<ALL_COMPS_JSON>>>;
const OVERVIEW    = <<<OVERVIEW_JSON>>>;
const KEGG_DATA   = <<<KEGG_JSON>>>;
const GO_DATA     = <<<GO_JSON>>>;
const KEGG_DOT    = <<<KEGG_DOT_JSON>>>;
const GO_DOT      = <<<GO_DOT_JSON>>>;
const MULTIMIR    = <<<MULTIMIR_JSON>>>;
const ROC_RVSNR   = <<<ROC_RVSNR_JSON>>>;
const ROC_UCHC    = <<<ROC_UCHC_JSON>>>;
const IBD_DATA    = <<<IBD_JSON>>>;

const CFG = {responsive:true, displaylogo:false,
             modeBarButtonsToRemove:['select2d','lasso2d','autoScale2d']};

/* ===== INIT CHARTS =====
   Wrapped so that a CDN/Plotly failure or a single bad chart can NEVER abort
   the rest of the script — tab switching and the image lightbox are defined
   below and must always work, even offline. */
const GO_ORDER = ['BP','CC','MF'];

function plot(id, fig) {
  try {
    if (typeof Plotly === 'undefined') return;
    if (fig && fig.data) Plotly.newPlot(id, fig.data, fig.layout, CFG);
  } catch (e) { console.error('Plotly init failed for', id, e); }
}

function initAllCharts() {
  plot('plt-overview', OVERVIEW);

  const compKeys = Object.keys(ALL_COMPS);
  if (compKeys.length) {
    const first = ALL_COMPS[compKeys[0]];
    plot('plt-volcano', first.volcano);
    plot('plt-lfc',     first.lfc);
  }

  plot('plt-kegg',     KEGG_DATA);
  plot('plt-kegg-dot', KEGG_DOT);
  plot('plt-multimir', MULTIMIR);
  if (ROC_RVSNR && ROC_RVSNR.data && ROC_RVSNR.data.length) plot('plt-roc-rvsnr', ROC_RVSNR);
  if (ROC_UCHC  && ROC_UCHC.data  && ROC_UCHC.data.length)  plot('plt-roc-uchc',  ROC_UCHC);
  if (IBD_DATA  && IBD_DATA.data  && IBD_DATA.data.length)   plot('plt-ibd',       IBD_DATA);

  for (const o of GO_ORDER) { if (GO_DOT[o])  { plot('plt-go-dot', GO_DOT[o]); break; } }
  for (const o of GO_ORDER) { if (GO_DATA[o]) { plot('plt-go',     GO_DATA[o]); break; } }
}

if (typeof Plotly === 'undefined') {
  /* CDN blocked or no internet — interactive charts cannot render, but the rest
     of the report (tables, static PNG figures with zoom) still works. */
  const b = document.createElement('div');
  b.style.cssText = 'background:#fff3cd;color:#664d03;border:1px solid #ffe69c;'
    + 'padding:12px 16px;margin:0 0 18px;border-radius:8px;font-size:.92em';
  b.innerHTML = '⚠️ Interactive charts could not load because the Plotly '
    + 'library (cdn.plot.ly) was unreachable. Static figures, tables, tab switching '
    + 'and image zoom still work. Reconnect to the internet and reload to see the charts.';
  const wrap = document.querySelector('.wrap') || document.body;
  wrap.insertBefore(b, wrap.firstChild);
} else {
  initAllCharts();
}

/* Shared GO state */
let _goOnto = 'BP', _goView = 'dot';

/* Ontology switcher — updates whichever chart type is currently active */
window.switchGOView = function(onto, btn) {
  _goOnto = onto;
  if (_goView === 'dot' && GO_DOT[onto])  Plotly.react('plt-go-dot', GO_DOT[onto].data,  GO_DOT[onto].layout,  CFG);
  if (_goView === 'bar' && GO_DATA[onto]) Plotly.react('plt-go',     GO_DATA[onto].data, GO_DATA[onto].layout, CFG);
  document.querySelectorAll('#go-tabs .tab').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
};

/* Chart-type toggle (Dotplot ↔ Bar chart) */
window.switchGOChartType = function(view, btn) {
  _goView = view;
  document.getElementById('go-panel-dot').style.display = view === 'dot' ? '' : 'none';
  document.getElementById('go-panel-bar').style.display = view === 'bar' ? '' : 'none';
  document.querySelectorAll('#go-view-tabs .tab').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
  requestAnimationFrame(function() {
    try {
      if (view === 'dot' && GO_DOT[_goOnto])  { Plotly.react('plt-go-dot', GO_DOT[_goOnto].data,  GO_DOT[_goOnto].layout,  CFG); Plotly.Plots.resize('plt-go-dot'); }
      if (view === 'bar' && GO_DATA[_goOnto]) { Plotly.react('plt-go',     GO_DATA[_goOnto].data, GO_DATA[_goOnto].layout, CFG); Plotly.Plots.resize('plt-go'); }
    } catch(e) {}
  });
};

/* KEGG chart-type toggle */
window.switchKeggView = function(view, btn) {
  document.getElementById('kegg-panel-dot').style.display = view === 'dot' ? '' : 'none';
  document.getElementById('kegg-panel-bar').style.display = view === 'bar' ? '' : 'none';
  document.querySelectorAll('#kegg-view-tabs .tab').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
  requestAnimationFrame(function() {
    try { Plotly.Plots.resize(view === 'dot' ? 'plt-kegg-dot' : 'plt-kegg'); } catch(e) {}
  });
};

/* ===== COMPARISON TABS ===== */
window.switchComp = function(key, btn) {
  const d = ALL_COMPS[key];
  if (!d) return;
  Plotly.react('plt-volcano', d.volcano.data, d.volcano.layout, CFG);
  Plotly.react('plt-lfc',     d.lfc.data,     d.lfc.layout,     CFG);
  document.querySelectorAll('#comp-tabs .tab').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
};
window.switchROC = function(which, btn) {
  document.getElementById('roc-panel-rvsnr').style.display = which==='rvsnr' ? '' : 'none';
  document.getElementById('roc-panel-uchc').style.display  = which==='uchc'  ? '' : 'none';
  document.querySelectorAll('#roc-tabs .tab').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
  requestAnimationFrame(function() {
    try {
      if (which==='rvsnr' && ROC_RVSNR.data && ROC_RVSNR.data.length) Plotly.Plots.resize('plt-roc-rvsnr');
      if (which==='uchc'  && ROC_UCHC.data  && ROC_UCHC.data.length)  Plotly.Plots.resize('plt-roc-uchc');
    } catch(e) {}
  });
};

window.switchTbl = function(key, btn) {
  document.querySelectorAll('.dea-tbl').forEach(t => t.classList.remove('visible'));
  const el = document.getElementById('tbl-' + key);
  if (el) el.classList.add('visible');
  document.querySelectorAll('#tbl-tabs .tab').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
};

/* ===== TABLE SORT ===== */
window.sortTbl = function(th, col) {
  const tbody = th.closest('table').querySelector('tbody');
  const rows  = Array.from(tbody.querySelectorAll('tr'));
  const asc   = th.dataset.asc !== '1';
  th.dataset.asc = asc ? '1' : '0';
  rows.sort((a,b)=>{
    const av = a.cells[col].innerText.trim();
    const bv = b.cells[col].innerText.trim();
    const an = parseFloat(av), bn = parseFloat(bv);
    const cmp = isNaN(an) ? av.localeCompare(bv) : an-bn;
    return asc ? cmp : -cmp;
  });
  rows.forEach(r => tbody.appendChild(r));
  th.closest('table').querySelectorAll('th').forEach(h => {
    h.textContent = h.textContent.replace(/[ ↑↓]$/,'');
  });
  th.textContent += asc ? ' ↑' : ' ↓';
};

/* ===== LIGHTBOX ===== */
(function(){
  const ov  = document.getElementById('lb-ov');
  const img = document.getElementById('lb-img');
  let sc=1, tx=0, ty=0, dr=false, ox=0, oy=0;
  const upd = () => { img.style.transform=`translate(${tx}px,${ty}px) scale(${sc})`; };
  const rst = () => { sc=1;tx=0;ty=0;upd(); };
  const open = src => { img.src=src; rst(); ov.classList.add('open'); document.body.style.overflow='hidden'; };
  const close = () => { ov.classList.remove('open'); document.body.style.overflow=''; };

  document.querySelectorAll('img.lb').forEach(el => el.addEventListener('click',()=>open(el.src)));
  document.getElementById('lb-x').addEventListener('click', close);
  ov.addEventListener('click', e => { if(e.target===ov) close(); });
  document.getElementById('lbzi').addEventListener('click', e => { e.stopPropagation(); sc=Math.min(sc*1.4,8); upd(); });
  document.getElementById('lbzo').addEventListener('click', e => { e.stopPropagation(); sc=Math.max(sc/1.4,.4); upd(); });
  document.getElementById('lbzr').addEventListener('click', e => { e.stopPropagation(); rst(); });
  img.addEventListener('wheel', e => { e.preventDefault(); sc=e.deltaY<0?Math.min(sc*1.12,8):Math.max(sc/1.12,.4); upd(); }, {passive:false});
  img.addEventListener('mousedown', e => { if(sc<=1) return; dr=true; ox=e.clientX-tx; oy=e.clientY-ty; img.classList.add('gr'); });
  document.addEventListener('mousemove', e => { if(!dr) return; tx=e.clientX-ox; ty=e.clientY-oy; upd(); });
  document.addEventListener('mouseup', () => { dr=false; img.classList.remove('gr'); });
  document.addEventListener('keydown', e => {
    if(!ov.classList.contains('open')) return;
    if(e.key==='Escape') close();
    if(e.key==='+'||e.key==='=') { sc=Math.min(sc*1.3,8); upd(); }
    if(e.key==='-') { sc=Math.max(sc/1.3,.4); upd(); }
    if(e.key==='0') rst();
  });
})();
</script>
</body>
</html>
"""


# ── HTML builders ─────────────────────────────────────────────────────────────

def build_qc_section(qc):
    if not qc.get("dir"):
        return ""

    n_samples = qc.get("n_samples", "?")
    n_mirnas  = qc.get("n_mirnas", "?")

    pca_img  = _lb_img(qc.get("pca"),     "PCA Plot",  "width:100%;border-radius:6px")
    umap_img = _lb_img(qc.get("umap"),    "UMAP Plot", "width:100%;border-radius:6px")
    heat_img = _lb_img(qc.get("heatmap"), "Heatmap",   "width:100%;border-radius:6px")

    vst_note = ""
    if qc.get("vst_path"):
        vst_note = (f'<div class="note" style="margin-top:16px">'
                    f'<strong>VST Normalized Matrix</strong> — '
                    f'{n_mirnas} miRNAs &times; {n_samples} samples. '
                    f'Used for PCA, UMAP and heatmap. '
                    f'File: <code>figures_and_QC/vst_normalized_matrix.csv</code></div>')

    return f"""
<div class="sec">
  <h2>Quality Control &amp; Sample Normalization</h2>
  <div class="note">
    Variance-Stabilizing Transformation (VST) applied before all QC plots.
    These plots show sample-level clustering and overall expression landscape
    <em>before</em> differential expression testing.
  </div>
  <div class="qc-grid">
    <div class="qc-card"><h3>PCA</h3>{pca_img}</div>
    <div class="qc-card"><h3>UMAP</h3>{umap_img}</div>
    <div class="qc-card"><h3>Heatmap — all DE miRNAs</h3>{heat_img}</div>
  </div>
  {vst_note}
</div>"""


def build_comp_tabs(comparisons, onclick="switchComp"):
    html = ""
    for i, c in enumerate(comparisons):
        active   = "active" if i == 0 else ""
        strict   = "strict" if c["model"] == "Strict" else ""
        html += (f'<button class="tab {active} {strict}" '
                 f'onclick="{onclick}(\'{c["key"]}\',this)">'
                 f'{c["short"]}<br>'
                 f'<small style="font-weight:400">{c["model"]}</small>'
                 f'</button>\n')
    return html


def build_comp_tables(comparisons):
    html = ""
    for i, c in enumerate(comparisons):
        visible = "visible" if i == 0 else ""
        tbl     = dea_table_html(c["dea_df"])
        html   += f'<div id="tbl-{c["key"]}" class="dea-tbl {visible}">{tbl}</div>\n'
    return html


def build_stats_rows(comparisons):
    rows = ""
    for c in comparisons:
        model_color = "#2e7bb5" if c["model"] == "Default" else "#e67e22"
        rows += (f'<tr>'
                 f'<td><strong>{c["short"]}</strong></td>'
                 f'<td><span style="color:{model_color};font-weight:600">{c["model"]}</span></td>'
                 f'<td class="num">{c["n_total"]}</td>'
                 f'<td class="num" style="color:#c0392b;font-weight:700">{c["n_up"]}</td>'
                 f'<td class="num" style="color:#2143a0;font-weight:700">{c["n_down"]}</td>'
                 f'</tr>\n')
    return rows


def build_dispersion(comparisons, model):
    for c in comparisons:
        if c["model"] == model and c.get("dispersion_path"):
            return _lb_img(c["dispersion_path"], f"Dispersion {model}",
                           "max-width:100%;border-radius:6px")
    return _missing("Dispersion plot not found")


def build_network_section(net_dir, parent_dir):
    if not net_dir:
        return _missing("Network dashboard files not found — run 16_interactive_network.R first"), ""

    # Relative path from the HTML file (which lives in parent_dir)
    net_rel = os.path.relpath(net_dir, parent_dir).replace("\\", "/")

    # Working interactive dashboards (vis-network "auto-freeze" tripartite views)
    # produced by 16_interactive_network.R — one per enrichment ontology.
    kegg_dash   = os.path.join(net_dir, "Automated_Dashboard_KEGG.html")
    go_dash     = os.path.join(net_dir, "Automated_Dashboard_GO.html")
    legacy_dash = os.path.join(net_dir, "Automated_Dashboard.html")

    links = ""
    if os.path.exists(kegg_dash):
        links += (f'<a href="{net_rel}/Automated_Dashboard_KEGG.html" class="btn">'
                  f'🕸 KEGG Network Dashboard</a>')
    elif os.path.exists(legacy_dash):
        # Fallback for results generated before the KEGG/GO split.
        links += (f'<a href="{net_rel}/Automated_Dashboard.html" class="btn">'
                  f'🕸 KEGG Network Dashboard</a>')
    if os.path.exists(go_dash):
        links += (f'<a href="{net_rel}/Automated_Dashboard_GO.html" class="btn">'
                  f'🕸 GO Network Dashboard</a>')
    if not links:
        links = _missing("Network dashboard files not found — run 16_interactive_network.R first")

    tri_path = os.path.join(net_dir, "Tripartite_Network_Plot.png")
    tri_html = ""
    if os.path.exists(tri_path):
        b64 = _file_b64(tri_path)
        if b64:
            tri_html = (f'<div style="margin-top:18px">'
                        f'<img class="lb" src="{b64}" alt="Tripartite Network" '
                        f'title="Click to enlarge" style="width:100%;border-radius:6px">'
                        f'</div>')

    return links, tri_html


def build_clinical_section(clin, score_by_source=None):
    """Build the three clinical HTML sections (ROC, IBD enrichment, HMDD/miRNet)."""
    if clin is None:
        return (
            '<div class="sec"><h2>Clinical Biomarker Analysis</h2>'
            + _missing("clinical_outputs/ folder not found — run clinical_roc.R, "
                       "hmdd_prep.R and ibd_target_overlap.R first, then regenerate this report.")
            + '</div>'
        )

    cand_df   = clin.get("candidates")
    roc_rvsnr = clin.get("roc_rvsnr")
    roc_uchc  = clin.get("roc_uchc")
    cv_model  = clin.get("cv_model")
    ibd_df    = clin.get("ibd_overlap")
    hmdd_df   = clin.get("hmdd_table")

    # ── mini stat cards ──────────────────────────────────────────────────────
    n_cand = len(cand_df) if cand_df is not None else "—"

    best_auc_str, best_mir = "—", ""
    if roc_rvsnr is not None and len(roc_rvsnr):
        best_auc_str = f'{roc_rvsnr.iloc[0]["AUC"]:.3f}'
        best_mir     = roc_rvsnr.iloc[0]["miRNA"]

    cv_auc_str = "—"
    if cv_model is not None and len(cv_model) and "AUC_CV" in cv_model.columns:
        cv_auc_str = f'{cv_model.iloc[0]["AUC_CV"]:.3f}'

    n_ibd_sig = "—"
    if ibd_df is not None and "padj" in ibd_df.columns:
        n_ibd_sig = str(int((ibd_df["padj"] <= 0.05).sum()))

    mini = (
        '<div style="display:grid;grid-template-columns:repeat(4,1fr);'
        'gap:12px;margin:16px 0 24px">'
        + "".join(
            f'<div style="background:#f8fafc;border-radius:8px;padding:16px;'
            f'text-align:center;border-top:3px solid {color}">'
            f'<div style="font-size:2em;font-weight:700;color:{color}">{val}</div>'
            f'<div style="color:#888;font-size:.83em;margin-top:5px">{lbl}</div></div>'
            for val, lbl, color in [
                (n_cand,       "Candidate miRNAs",                            "#2e7bb5"),
                (best_auc_str, f"Best AUC — R vs NR<br><em>{best_mir}</em>", "#d62728"),
                (cv_auc_str,   "CV-AUC<br><em>multi-miRNA model</em>",       "#27ae60"),
                (n_ibd_sig,    "miRNAs enriched<br><em>in IBD genes FDR≤0.05</em>", "#9b59b6"),
            ]
        )
        + "</div>"
    )

    # ── ROC image galleries ───────────────────────────────────────────────────
    def roc_gallery(pngs, n=6):
        items = "".join(
            f'<div style="text-align:center">'
            f'<p style="font-size:.78em;color:#666;margin-bottom:4px">'
            f'{os.path.basename(p).replace(".png","").replace("ROC_RvsNR_","").replace("ROC_UCvsHC_","").replace("_","-")}'
            f'</p>'
            f'{_lb_img(p, os.path.basename(p), "width:100%")}'
            f'</div>'
            for p in pngs[:n]
        )
        if not items:
            items = _missing("ROC PNGs not found — run clinical_roc.R first")
        return (
            '<div style="display:grid;grid-template-columns:repeat(3,1fr);gap:12px">'
            + items + "</div>"
        )

    rvsnr_gallery = roc_gallery(clin.get("roc_rvsnr_pngs", []))
    uchc_gallery  = roc_gallery(clin.get("roc_uchc_pngs",  []))

    # ── Multi-miRNA CV model card ─────────────────────────────────────────────
    multi_html = ""
    if cv_model is not None and len(cv_model):
        r          = cv_model.iloc[0]
        mirs_str   = str(r.get("miRNAs", "—")).replace(";", ", ")
        n_samples  = int(r.get("n", 0))
        multi_roc  = _lb_img(
            clin.get("roc_multi"), "Multi-miRNA CV ROC",
            "width:100%;max-width:420px;border-radius:6px"
        )
        multi_html = (
            '<div class="note" style="margin-top:20px">'
            f'<strong>Multi-miRNA 5-fold CV Logistic Model</strong>'
            f'&emsp;AUC<sub>CV</sub> = <strong>{cv_auc_str}</strong>'
            f'&emsp;|&emsp; n = {n_samples} samples<br>'
            f'<em>miRNAs: {mirs_str}</em></div>'
            f'<div style="text-align:center;max-width:420px;margin:14px auto">{multi_roc}</div>'
        )

    # ── Violin gallery ────────────────────────────────────────────────────────
    violin_pngs = clin.get("violin_pngs", [])
    if violin_pngs:
        vitems = "".join(
            f'<div style="text-align:center">'
            f'<p style="font-size:.78em;color:#666;margin-bottom:4px">'
            f'{os.path.basename(p).replace("violin_","").replace(".png","").replace("_","-")}'
            f'</p>{_lb_img(p, os.path.basename(p), "width:100%")}</div>'
            for p in violin_pngs[:12]
        )
        violin_html = (
            '<h3 style="margin:24px 0 10px;color:#1a3c5e">'
            'Expression by Group — HC / NR baseline / R baseline</h3>'
            '<div class="note">Violin + boxplot for each candidate miRNA. Click any image to enlarge.</div>'
            '<div style="display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin-top:12px">'
            + vitems + "</div>"
        )
    else:
        violin_html = _missing("Violin plots not found — run clinical_roc.R first")

    # ── HMDD / miRNet card ────────────────────────────────────────────────────
    dir_clin  = clin["dir"]
    n_up      = int((hmdd_df["direction"] == "UP").sum())   if hmdd_df is not None else "?"
    n_down    = int((hmdd_df["direction"] == "DOWN").sum()) if hmdd_df is not None else "?"
    hmdd_html = (
        '<p style="color:#666;margin-bottom:14px">'
        'Lists generated by <code>hmdd_prep.R</code>: '
        f'<code>clinical_outputs/hmdd_lists.txt</code> (copy-paste) &amp; '
        f'<code>clinical_outputs/mirnet_input.tsv</code> (upload file).</p>'
        f'<p style="color:#666;font-size:.9em;margin-bottom:16px">'
        f'{len(hmdd_df) if hmdd_df is not None else "?"} candidates: '
        f'{n_up} upregulated · {n_down} downregulated in UC vs HC</p>'
        '<div class="two-col">'
        '<div class="note"><strong>HMDD v4.0</strong>'
        '&ensp;<a href="https://www.cuilab.cn/hmdd" target="_blank" class="btn"'
        ' style="padding:5px 14px;font-size:.85em">Open →</a><br>'
        '<small>Batch search: paste Section 1 of <em>hmdd_lists.txt</em> (max 20 keywords)<br>'
        'Disease network: UP/DOWN lists in Section 2</small></div>'
        '<div class="note"><strong>miRNet 2.0</strong>'
        '&ensp;<a href="https://www.mirnet.ca/miRNet/upload/MirUploadView.xhtml"'
        ' target="_blank" class="btn" style="padding:5px 14px;font-size:.85em">Open →</a><br>'
        '<small>Upload <em>mirnet_input.tsv</em> → miRNA → miRBase → <em>Homo sapiens</em></small>'
        '</div></div>'
    )

    # ── IBD table ─────────────────────────────────────────────────────────────
    ibd_tbl   = ibd_table_html(ibd_df)
    ibd_fig   = _lb_img(
        clin.get("ibd_fig"), "IBD Enrichment",
        "width:100%;border-radius:6px"
    )

    # ── Assemble ──────────────────────────────────────────────────────────────
    cand_tbl = cand_table_html(cand_df, score_by_source)

    return f"""
<!-- ── CLINICAL BIOMARKER ANALYSIS ── -->
<div class="sec">
  <h2>Clinical Biomarker Analysis</h2>
  <div class="note">
    Candidate miRNAs: padj &le; 0.05, |log&#8322;FC| &ge; 1 (either DEA comparison, Default model).
    <strong>Contrast 1 (primary):</strong> Responder vs Non-Responder at baseline —
    can miRNA expression <em>before</em> treatment predict glucocorticoid response?
    &ensp;<strong>Contrast 2:</strong> UC baseline vs Healthy Controls.
  </div>

  {mini}

  <h3 style="color:#1a3c5e;margin-bottom:12px">Per-miRNA ROC / AUC</h3>
  <div class="tabs" id="roc-tabs">
    <button class="tab active" onclick="switchROC('rvsnr',this)">R vs NR &mdash; Responder Prediction</button>
    <button class="tab"        onclick="switchROC('uchc',this)">UC vs HC &mdash; Disease Detection</button>
  </div>

  <div id="roc-panel-rvsnr">
    <div class="two-col">
      <div><div id="plt-roc-rvsnr"></div></div>
      <div><h3>Top-6 ROC curves (R vs NR)</h3>{rvsnr_gallery}</div>
    </div>
  </div>

  <div id="roc-panel-uchc" style="display:none">
    <div class="two-col">
      <div><div id="plt-roc-uchc"></div></div>
      <div><h3>Top-6 ROC curves (UC vs HC)</h3>{uchc_gallery}</div>
    </div>
  </div>

  {multi_html}
  {violin_html}

  <h3 style="color:#1a3c5e;margin:24px 0 10px">Candidate miRNAs ({n_cand} total)</h3>
  <div class="note">Sorted by padj. Direction is relative to UC vs Healthy Controls.</div>
  {cand_tbl}
</div>

<!-- ── IBD TARGET ENRICHMENT ── -->
<div class="sec">
  <h2>IBD Target Gene Enrichment</h2>
  <div class="note">
    Fisher&apos;s exact test (one-tailed, BH FDR) per miRNA: are validated targets
    enriched in known IBD / UC / Crohn&apos;s disease genes (OpenTargets Platform or
    curated GWAS fallback)? &ensp; ★ = FDR &le; 0.05.
  </div>
  <div class="two-col">
    <div><div id="plt-ibd"></div></div>
    <div style="text-align:center">{ibd_fig}</div>
  </div>
  <div style="margin-top:20px">{ibd_tbl}</div>
</div>

<!-- ── EXTERNAL DATABASE SUBMISSIONS ── -->
<div class="sec">
  <h2>External Database Submissions &mdash; HMDD &amp; miRNet</h2>
  {hmdd_html}
</div>"""


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parent_dir = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
    print(f"Results dir : {parent_dir}")

    # Derive a dataset label from the folder name (e.g. "PRJNA331127", or a legacy "a_PRJNA331127" -> "PRJNA331127")
    m = re.search(r"(PRJ[A-Z]{2}\d+|GSE\d+)", os.path.basename(parent_dir), re.I)
    dataset_label = m.group(1) if m else os.path.basename(parent_dir)

    # Discover
    comparisons, qc, default_dir = discover_all(parent_dir)
    if not comparisons:
        sys.exit("ERROR: No DEA_* sub-folders with DEA_*.csv files found.")
    print(f"  Found {len(comparisons)} comparison(s): "
          + ", ".join(c['label'] for c in comparisons))

    # Locate downstream output folders anywhere under the dataset root.
    # Stages 14/15/16/19 write to the dataset root (getwd()), not into the DEA folder.
    pathway_dir  = find_dir(parent_dir, "pathway_outputs")
    multimir_dir = find_dir(parent_dir, "multimir_outputs")
    network_dir  = find_dir(parent_dir, "network_outputs")
    for label, d in (("pathway", pathway_dir), ("multimir", multimir_dir), ("network", network_dir)):
        print(f"  {label}_outputs : {os.path.relpath(d, parent_dir) if d else 'not found'}")

    kegg_df     = load_kegg(pathway_dir)
    go_df       = load_go(pathway_dir)
    multimir_df = load_multimir(multimir_dir)
    n_targets   = count_targets(multimir_dir)

    print("  Clinical analysis …")
    clin_data   = load_clinical(parent_dir)
    if clin_data:
        keys_found = [k for k, v in clin_data.items()
                      if isinstance(v, pd.DataFrame) and v is not None]
        print(f"    Found: {', '.join(keys_found)}")
    else:
        print("    clinical_outputs/ not found — clinical sections will show placeholders")

    # Build Plotly JSON for all comparisons
    print("  Building comparison charts …")
    all_comps_dict = {}
    for c in comparisons:
        all_comps_dict[c["key"]] = {
            "volcano": volcano_json(c["dea_df"], c["label"]),
            "lfc":     lfc_json(c["dea_df"], c["label"]),
        }

    # miRNA-score (S_miRNA) lookup keyed by source token, for the candidate table.
    # Token = last word of the comparison short name, e.g. "UC Responder" -> "responder",
    # matching the `source` column written by clinical_roc.R ("Responder"/"NonResponder").
    score_by_source = {}
    for c in comparisons:
        if c["model"] != "Default":
            continue
        df = c["dea_df"]
        if "S_miRNA" in df.columns:
            token = c["short"].split()[-1].lower()
            score_by_source[token] = dict(zip(df["miRNA"], df["S_miRNA"]))

    print("  Overview chart …")
    ov_js   = overview_json(comparisons)
    print("  KEGG …")
    k_js    = kegg_json(kegg_df) if kegg_df is not None else '{"data":[],"layout":{}}'
    print("  GO …")
    go_dict = go_json(go_df) if go_df is not None else {}
    go_js   = "{" + ",".join(f'"{k}":{v}' for k, v in go_dict.items()) + "}"
    print("  KEGG dotplot …")
    kd_js        = kegg_dotplot_json(kegg_df) if kegg_df is not None else '{"data":[],"layout":{}}'
    print("  GO dotplot …")
    go_dot_dict  = go_dotplot_json(go_df) if go_df is not None else {}
    go_dot_js    = "{" + ",".join(f'"{k}":{v}' for k, v in go_dot_dict.items()) + "}"
    print("  multiMiR …")
    m_js    = multimir_json(multimir_df) if multimir_df is not None else '{"data":[],"layout":{}}'
    print("  Clinical ROC / IBD charts …")
    roc_rvsnr_js = roc_bar_json(
        clin_data.get("roc_rvsnr") if clin_data else None,
        "Responder vs Non-Responder (baseline)")
    roc_uchc_js  = roc_bar_json(
        clin_data.get("roc_uchc")  if clin_data else None,
        "UC baseline vs Healthy Controls")
    ibd_js       = ibd_bar_json(
        clin_data.get("ibd_overlap") if clin_data else None)

    # HTML fragments
    qc_section   = build_qc_section(qc)
    comp_tabs    = build_comp_tabs(comparisons, "switchComp")
    tbl_tabs     = build_comp_tabs(comparisons, "switchTbl")
    comp_tables  = build_comp_tables(comparisons)
    stats_rows   = build_stats_rows(comparisons)
    disp_default = build_dispersion(comparisons, "Default")
    disp_strict  = build_dispersion(comparisons, "Strict")

    # KEGG section — dotplot primary, bar chart behind toggle, static PNG always shown
    if kegg_df is not None:
        kegg_content = (
            '<div class="tabs" id="kegg-view-tabs" style="margin-bottom:14px">'
            '<button class="tab active" onclick="switchKeggView(\'dot\',this)">Dotplot</button>'
            '<button class="tab"        onclick="switchKeggView(\'bar\',this)">Bar chart</button>'
            '</div>'
            '<div id="kegg-panel-dot"><div id="plt-kegg-dot"></div></div>'
            '<div id="kegg-panel-bar" style="display:none"><div id="plt-kegg"></div></div>'
            '<div style="margin-top:24px">'
            '<h3 style="color:#555;font-size:.95em;margin-bottom:10px">'
            'R / clusterProfiler Static Dotplot</h3>'
            '<<<KEGG_DOTPLOT>>>'
            '</div>'
        )
    else:
        kegg_content = _missing("KEGG enrichment results not found — run stage 15 first")

    # GO section — dotplot primary, bar chart behind toggle, static PNG always shown
    if go_dot_dict or go_dict:
        go_content = (
            '<div class="tabs" id="go-tabs">'
            '<button class="tab active" onclick="switchGOView(\'BP\',this)">Biological Process</button>'
            '<button class="tab"        onclick="switchGOView(\'CC\',this)">Cellular Component</button>'
            '<button class="tab"        onclick="switchGOView(\'MF\',this)">Molecular Function</button>'
            '</div>'
            '<div class="tabs" id="go-view-tabs" style="margin-top:8px">'
            '<button class="tab active" onclick="switchGOChartType(\'dot\',this)">Dotplot</button>'
            '<button class="tab"        onclick="switchGOChartType(\'bar\',this)">Bar chart</button>'
            '</div>'
            '<div id="go-panel-dot"><div id="plt-go-dot"></div></div>'
            '<div id="go-panel-bar" style="display:none"><div id="plt-go"></div></div>'
            '<div style="margin-top:24px">'
            '<h3 style="color:#555;font-size:.95em;margin-bottom:10px">'
            'R / clusterProfiler Static Dotplot</h3>'
            '<<<GO_DOTPLOT>>>'
            '</div>'
        )
    else:
        go_content = _missing("GO enrichment results not found — run stage 15 first")
    multimir_content = ('<div id="plt-multimir"></div>' if multimir_df is not None
                        else _missing("multiMiR results not found — run stage 14 first"))

    net_links, tri_img = build_network_section(network_dir, parent_dir)

    # Dotplots (embedded with lightbox)
    kegg_dotplot = _lb_img(
        os.path.join(pathway_dir, "KEGG_Dotplot.png")
        if pathway_dir else None, "KEGG Dotplot")
    go_dotplot = _lb_img(
        os.path.join(pathway_dir, "GO_Dotplot.png")
        if pathway_dir else None, "GO Dotplot")

    # Stat cards
    n_samples_str = str(qc.get("n_samples", "—"))
    n_total_str   = str(comparisons[0]["n_total"]) if comparisons else "—"

    # Assemble HTML
    html = (
        _TEMPLATE
        .replace("<<<DATASET>>>",        dataset_label)
        .replace("<<<DATE>>>",           datetime.now().strftime("%Y-%m-%d %H:%M"))
        .replace("<<<N_TOTAL>>>",        n_total_str)
        .replace("<<<N_COMPS>>>",        str(len(comparisons)))
        .replace("<<<N_SAMPLES>>>",      n_samples_str)
        .replace("<<<N_TARGETS>>>",      f"{n_targets:,}" if n_targets else "—")
        .replace("<<<STATS_ROWS>>>",     stats_rows)
        .replace("<<<QC_SECTION>>>",     qc_section)
        .replace("<<<COMP_TABS>>>",      comp_tabs)
        .replace("<<<TBL_TABS>>>",       tbl_tabs)
        .replace("<<<COMP_TABLES>>>",    comp_tables)
        .replace("<<<DISP_DEFAULT>>>",   disp_default)
        .replace("<<<DISP_STRICT>>>",    disp_strict)
        .replace("<<<KEGG_CONTENT>>>",   kegg_content)
        .replace("<<<GO_CONTENT>>>",     go_content)
        .replace("<<<KEGG_DOTPLOT>>>",   kegg_dotplot)
        .replace("<<<GO_DOTPLOT>>>",     go_dotplot)
        .replace("<<<MULTIMIR_CONTENT>>>", multimir_content)
        .replace("<<<NETWORK_LINKS>>>",  net_links)
        .replace("<<<TRIPARTITE_IMG>>>", tri_img)
        .replace("<<<ALL_COMPS_JSON>>>", json.dumps(all_comps_dict))
        .replace("<<<OVERVIEW_JSON>>>",  ov_js)
        .replace("<<<KEGG_JSON>>>",      k_js)
        .replace("<<<GO_JSON>>>",        go_js)
        .replace("<<<KEGG_DOT_JSON>>>",  kd_js)
        .replace("<<<GO_DOT_JSON>>>",    go_dot_js)
        .replace("<<<MULTIMIR_JSON>>>",  m_js)
        .replace("<<<CLINICAL_SECTION>>>", build_clinical_section(clin_data, score_by_source))
        .replace("<<<ROC_RVSNR_JSON>>>",   roc_rvsnr_js)
        .replace("<<<ROC_UCHC_JSON>>>",    roc_uchc_js)
        .replace("<<<IBD_JSON>>>",         ibd_js)
    )

    out = os.path.join(parent_dir, "Final_Report.html")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)

    size_kb = os.path.getsize(out) // 1024
    print(f"\nReport : {out}  ({size_kb} KB)")
    print(f"Open   : file:///{out.replace(os.sep, '/')}")


if __name__ == "__main__":
    main()
