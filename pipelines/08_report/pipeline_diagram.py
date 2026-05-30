#!/usr/bin/env python3
"""
Pipeline diagram generator for thesis.

Generates a high-resolution figure of the complete miRNA-seq IBD pipeline,
showing all stages, scripts, inputs, and outputs end-to-end.

Output: pipeline_diagram.png  (300 DPI, ~14 × 22 in)
        pipeline_diagram.pdf  (vector, for thesis)

Usage:
    py -3.11 pipeline_diagram.py [output_folder]
    Defaults to the same folder as this script.
"""

import sys
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

# ── Palette ───────────────────────────────────────────────────────────────────
P = {
    # Stage colours
    "s00":  "#1B4F72",   # dark navy   — search/discovery
    "s01":  "#1F618D",   # blue        — download
    "s02":  "#2471A3",   # mid-blue    — QC
    "s03":  "#117A65",   # dark teal   — quantification
    "s04":  "#148F77",   # teal        — merge/filter
    "s05":  "#6C3483",   # purple      — DEA
    "s06":  "#A04000",   # brown-red   — targets
    "s07a": "#B7770D",   # amber       — pathways
    "s07b": "#C0392B",   # red         — network/clinical
    "lit":  "#515A5A",   # dark gray   — literature
    "rep":  "#1A5276",   # dark blue   — report
    # Supporting element colours
    "db":   "#2C3E50",   # database cylinders
    "data": "#EBF5FB",   # data flow boxes
    "data_border": "#AED6F1",
    "out":  "#FEF9E7",   # output boxes
    "out_border": "#F9E79F",
    "arrow": "#566573",
    "bg":    "#FDFEFE",
    "title": "#17202A",
    "white": "#FFFFFF",
    "grid":  "#EAF2F8",
}

# ── Layout constants ──────────────────────────────────────────────────────────
FW, FH = 16, 24          # figure inches
BOX_W  = 7.2             # main stage box width
BOX_H  = 1.05            # main stage box height
CX     = 4.6             # centre-x of main pipeline
LIT_X  = 12.2            # centre-x of literature module
LIT_W  = 3.2
DATA_W = 3.8             # data-label box width
DATA_H = 0.55
FS_TITLE = 11.5          # stage title font size
FS_TOOL  = 8.0           # tool name font size
FS_IO    = 7.8           # input/output label font size
FS_STAGE = 8.5           # stage number font size

# ── Y-positions (top of each element, in figure inches) ──────────────────────
#   The figure coordinate system is inches from bottom.
#   We'll use ax.transData with xlim/ylim matching figure units.

STAGES_Y = {
    "header":  22.6,
    "db":      21.2,
    "s00":     19.4,
    "d_s00":   18.25,   # data label after s00
    "s01":     17.3,
    "d_s01":   16.15,
    "s02":     15.2,
    "d_s02":   14.05,
    "s03":     13.1,
    "d_s03":   11.95,
    "s04":     11.0,
    "d_s04":    9.85,
    "s05":      8.9,
    "d_s05":    7.75,
    "s06":      6.8,
    "d_s06":    5.65,
    "s07a":     4.7,
    "d_s07a":   3.55,
    "s07b":     2.6,
    "d_s07b":   1.5,
    "report":   0.45,
}


# ── Drawing helpers ───────────────────────────────────────────────────────────

def stage_box(ax, cx, cy, w, h, num, title, tools, scripts, color, alpha=1.0):
    """Draw a stage box with stage number badge, title, tools, and script names."""
    x = cx - w / 2
    y = cy - h / 2

    # Shadow
    shadow = FancyBboxPatch((x + 0.06, y - 0.06), w, h,
                             boxstyle="round,pad=0.08",
                             facecolor="#CACFD2", edgecolor="none",
                             linewidth=0, zorder=2, alpha=0.5)
    ax.add_patch(shadow)

    # Main box
    box = FancyBboxPatch((x, y), w, h,
                          boxstyle="round,pad=0.08",
                          facecolor=color, edgecolor="white",
                          linewidth=1.5, zorder=3, alpha=alpha)
    ax.add_patch(box)

    # Left accent strip (stage number)
    badge_w = 0.72
    badge = FancyBboxPatch((x, y), badge_w, h,
                            boxstyle="round,pad=0.08",
                            facecolor="white", edgecolor="none",
                            linewidth=0, zorder=4, alpha=0.18)
    ax.add_patch(badge)
    ax.text(x + badge_w / 2, cy, num,
            ha="center", va="center", fontsize=FS_STAGE, fontweight="bold",
            color="white", zorder=5, alpha=0.95)

    # Title
    ax.text(x + badge_w + (w - badge_w) * 0.45, cy + h * 0.18,
            title, ha="center", va="center",
            fontsize=FS_TITLE, fontweight="bold", color="white", zorder=5)

    # Tools / scripts row
    detail = "  ·  ".join(tools)
    if scripts:
        detail += f"   [{', '.join(scripts)}]"
    ax.text(x + badge_w + (w - badge_w) * 0.45, cy - h * 0.22,
            detail, ha="center", va="center",
            fontsize=FS_TOOL, color="white", alpha=0.82, zorder=5)


def data_box(ax, cx, cy, w, h, text, color=None, border=None):
    """Draw a lightweight data/file label box."""
    fc = color  or P["data"]
    ec = border or P["data_border"]
    x  = cx - w / 2
    y  = cy - h / 2
    box = FancyBboxPatch((x, y), w, h,
                          boxstyle="round,pad=0.06",
                          facecolor=fc, edgecolor=ec,
                          linewidth=1.0, zorder=3)
    ax.add_patch(box)
    ax.text(cx, cy, text, ha="center", va="center",
            fontsize=FS_IO, color=P["title"], zorder=4, style="italic")


def output_box(ax, cx, cy, w, h, text):
    data_box(ax, cx, cy, w, h, text, color=P["out"], border=P["out_border"])


def db_cylinder(ax, cx, cy, w=2.4, h=0.9, label="GEO / SRA\nNCBI Databases"):
    """Draw a simplified database 'cylinder'."""
    x = cx - w / 2
    y = cy - h / 2

    body = FancyBboxPatch((x, y), w, h,
                           boxstyle="round,pad=0.1",
                           facecolor=P["db"], edgecolor="#ABB2B9",
                           linewidth=1.2, zorder=3)
    ax.add_patch(body)

    # Ellipse top
    from matplotlib.patches import Ellipse
    ell = Ellipse((cx, cy + h / 2), width=w, height=0.3,
                  facecolor="#3D5166", edgecolor="#ABB2B9",
                  linewidth=1.0, zorder=4)
    ax.add_patch(ell)

    ax.text(cx, cy + 0.05, label, ha="center", va="center",
            fontsize=9.5, fontweight="bold", color="white", zorder=5,
            multialignment="center")


def arrow(ax, x1, y1, x2, y2, color=None, lw=1.6, head=0.18):
    c = color or P["arrow"]
    ax.annotate("",
                xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(
                    arrowstyle=f"->,head_width={head},head_length={head*0.7}",
                    color=c, lw=lw, connectionstyle="arc3,rad=0"))


def dashed_arrow(ax, x1, y1, x2, y2, color=None, lw=1.2):
    c = color or "#7F8C8D"
    ax.annotate("",
                xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(
                    arrowstyle="->,head_width=0.14,head_length=0.12",
                    color=c, lw=lw,
                    connectionstyle="arc3,rad=0",
                    linestyle="dashed"))


# ── Main drawing ──────────────────────────────────────────────────────────────

def draw():
    fig, ax = plt.subplots(figsize=(FW, FH))
    ax.set_xlim(0, FW)
    ax.set_ylim(0, FH)
    ax.set_facecolor(P["bg"])
    fig.patch.set_facecolor(P["bg"])
    ax.axis("off")

    # ── Background sections ───────────────────────────────────────────────────
    # Colour-tinted bands per pipeline section
    sections = [
        (22.0, 19.9, "#EBF5FB", "DATA ACQUISITION"),  # blue tint
        (19.9, 13.7, "#E8F8F5", "SEQUENCING PROCESSING"),  # green tint
        ( 9.5, 13.7, "#F5EEF8", "STATISTICAL ANALYSIS"),   # purple tint
        ( 1.0,  9.5, "#FEF5E7", "BIOLOGICAL INTERPRETATION"),  # orange tint
        ( 0.0,  1.0, "#EBF5FB", "REPORTING"),
    ]
    for y_top, y_bot, fc, label in sections:
        h = y_top - y_bot
        rect = mpatches.FancyBboxPatch(
            (0.18, y_bot), FW - 0.36, h,
            boxstyle="round,pad=0.08",
            facecolor=fc, edgecolor="none", zorder=0, alpha=0.5)
        ax.add_patch(rect)
        ax.text(0.42, (y_top + y_bot) / 2, label,
                ha="left", va="center", fontsize=7.5,
                color="#99A3A4", fontweight="bold", rotation=90, zorder=1,
                alpha=0.7)

    # ── Title ─────────────────────────────────────────────────────────────────
    ax.text(FW / 2, STAGES_Y["header"],
            "miRNA-seq Analysis Pipeline for IBD",
            ha="center", va="center", fontsize=17, fontweight="bold",
            color=P["title"], zorder=10)
    ax.text(FW / 2, STAGES_Y["header"] - 0.55,
            "Dataset PRJNA471862  ·  Ulcerative Colitis  ·  End-to-End Workflow",
            ha="center", va="center", fontsize=9.5,
            color="#566573", zorder=10)

    # ── Database ──────────────────────────────────────────────────────────────
    db_cylinder(ax, CX, STAGES_Y["db"], w=2.6, h=0.82)
    # PubMed on the right side for literature
    db_cylinder(ax, LIT_X, STAGES_Y["db"], w=2.5, h=0.82,
                label="PubMed / NCBI\nLiterature DB")

    # Arrow: DB → S00
    arrow(ax, CX, STAGES_Y["db"] - 0.45, CX, STAGES_Y["s00"] + BOX_H / 2 + 0.05)

    # ── Stage 00 ──────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s00"], BOX_W, BOX_H,
              "00", "Dataset Discovery & Metadata Curation",
              ["EDirect", "SRA Toolkit"],
              ["01_geo_ibd_curator.py", "02_bioproject_curator.py"],
              P["s00"])

    data_box(ax, CX, STAGES_Y["d_s00"], DATA_W + 1.0, DATA_H,
             "candidates_master.csv  ·  runinfo.csv  ·  dataset metadata")
    arrow(ax, CX, STAGES_Y["s00"] - BOX_H / 2, CX, STAGES_Y["d_s00"] + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["d_s00"] - DATA_H / 2, CX, STAGES_Y["s01"] + BOX_H / 2 + 0.05)

    # ── Stage 01 ──────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s01"], BOX_W, BOX_H,
              "01", "Raw Data Download",
              ["prefetch", "fasterq-dump", "pigz"],
              ["03_download_runs.sh"],
              P["s01"])

    data_box(ax, CX, STAGES_Y["d_s01"], DATA_W + 0.6, DATA_H,
             "Raw reads  ·  SRR*.fastq.gz  (160 samples)")
    arrow(ax, CX, STAGES_Y["s01"] - BOX_H / 2, CX, STAGES_Y["d_s01"] + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["d_s01"] - DATA_H / 2, CX, STAGES_Y["s02"] + BOX_H / 2 + 0.05)

    # ── Stage 02 ──────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s02"], BOX_W, BOX_H,
              "02", "Quality Control",
              ["FastQC", "MultiQC", "miRTrace"],
              ["04_qc_fastqc.sh", "06_qc_mirtrace.sh"],
              P["s02"])

    data_box(ax, CX, STAGES_Y["d_s02"], DATA_W + 0.8, DATA_H,
             "FastQC / MultiQC report  ·  miRTrace QC report")
    arrow(ax, CX, STAGES_Y["s02"] - BOX_H / 2, CX, STAGES_Y["d_s02"] + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["d_s02"] - DATA_H / 2, CX, STAGES_Y["s03"] + BOX_H / 2 + 0.05)

    # ── Stage 03 ──────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s03"], BOX_W, BOX_H,
              "03", "miRNA Quantification",
              ["miRge3.0", "miRBase hsa"],
              ["08_test_auto_mirge3_subfolders.sh", "09_single_end_runs_mirge3.sh"],
              P["s03"])

    data_box(ax, CX, STAGES_Y["d_s03"], DATA_W + 1.0, DATA_H,
             "miRNA.Counts.csv  ·  miRNA.RPM.csv  (per sample)")
    arrow(ax, CX, STAGES_Y["s03"] - BOX_H / 2, CX, STAGES_Y["d_s03"] + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["d_s03"] - DATA_H / 2, CX, STAGES_Y["s04"] + BOX_H / 2 + 0.05)

    # ── Stage 04 ──────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s04"], BOX_W, BOX_H,
              "04", "Count Merge & Filtering",
              ["pandas", "mean RPM ≥ 100 filter"],
              ["11_mrina_cut_off.py", "12_deseq2_metadata.py"],
              P["s04"])

    data_box(ax, CX, STAGES_Y["d_s04"], DATA_W + 1.4, DATA_H,
             "Master_Filtered_Counts_DESeq2.csv  ·  sample_metadata.csv")
    arrow(ax, CX, STAGES_Y["s04"] - BOX_H / 2, CX, STAGES_Y["d_s04"] + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["d_s04"] - DATA_H / 2, CX, STAGES_Y["s05"] + BOX_H / 2 + 0.05)

    # ── Stage 05 ──────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s05"], BOX_W, BOX_H,
              "05", "Differential Expression Analysis",
              ["DESeq2", "apeglm", "VST", "PCA / UMAP"],
              ["13_run_dea.R"],
              P["s05"])

    # Two outputs from DEA: results + QC figures
    dy = STAGES_Y["d_s05"]
    data_box(ax, CX - 1.8, dy, DATA_W, DATA_H,
             "DEA results  ·  Default & Strict models")
    data_box(ax, CX + 2.1, dy, DATA_W - 0.4, DATA_H,
             "PCA · UMAP · Heatmap · Volcano")

    arrow(ax, CX, STAGES_Y["s05"] - BOX_H / 2, CX - 1.8, dy + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["s05"] - BOX_H / 2, CX + 2.1, dy + DATA_H / 2)
    arrow(ax, CX - 1.8, dy - DATA_H / 2, CX, STAGES_Y["s06"] + BOX_H / 2 + 0.05)

    # ── Stage 06 ──────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s06"], BOX_W, BOX_H,
              "06", "Target Prediction",
              ["multiMiR", "miRTarBase", "TarBase", "miRecords"],
              ["14_multimir_targets.R"],
              P["s06"])

    data_box(ax, CX, STAGES_Y["d_s06"], DATA_W + 1.6, DATA_H,
             "targets_validated.csv  ·  targets_predicted.csv")
    arrow(ax, CX, STAGES_Y["s06"] - BOX_H / 2, CX, STAGES_Y["d_s06"] + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["d_s06"] - DATA_H / 2, CX, STAGES_Y["s07a"] + BOX_H / 2 + 0.05)

    # ── Stage 07a ─────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s07a"], BOX_W, BOX_H,
              "07a", "Pathway Enrichment",
              ["clusterProfiler", "GO", "KEGG"],
              ["15_pathway_enrich.R"],
              P["s07a"])

    data_box(ax, CX, STAGES_Y["d_s07a"], DATA_W + 1.2, DATA_H,
             "GO_Enrichment.csv  ·  KEGG_Enrichment.csv  ·  dotplots")
    arrow(ax, CX, STAGES_Y["s07a"] - BOX_H / 2, CX, STAGES_Y["d_s07a"] + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["d_s07a"] - DATA_H / 2, CX, STAGES_Y["s07b"] + BOX_H / 2 + 0.05)

    # ── Stage 07b ─────────────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["s07b"], BOX_W, BOX_H,
              "07b", "Network Construction & Clinical Evaluation",
              ["igraph", "clusterProfiler", "pROC", "Louvain"],
              ["build_network.R", "clinical_eval.R", "16_interactive_network.R"],
              P["s07b"])

    data_box(ax, CX, STAGES_Y["d_s07b"], DATA_W + 1.4, DATA_H,
             "miRNA–gene network  ·  ROC/AUC  ·  centrality scores")
    arrow(ax, CX, STAGES_Y["s07b"] - BOX_H / 2, CX, STAGES_Y["d_s07b"] + DATA_H / 2)
    arrow(ax, CX, STAGES_Y["d_s07b"] - DATA_H / 2,
          CX, STAGES_Y["report"] + BOX_H / 2 + 0.05)

    # ── Stage 08 — Report ─────────────────────────────────────────────────────
    stage_box(ax, CX, STAGES_Y["report"], BOX_W, BOX_H,
              "08", "Automated HTML Report  &  Interactive Dashboard",
              ["Plotly.js", "vis.js", "matplotlib"],
              ["17_generate_report.py", "16_interactive_network.R"],
              P["rep"])

    # ── Literature Mining module ───────────────────────────────────────────────
    lit_ys = [
        (20.8, "PubMed Search\npt_search.py · PubMed_API_0.1.py"),
        (19.1, "Paper Retrieval\n& Deduplication\nrefineria.py"),
        (17.3, "Local TF-IDF\nSearch Engine\ninteractive_cli.py"),
        (15.6, "Filtered corpus\n(papers.csv / .jsonl)"),
    ]

    for i, (ly, label) in enumerate(lit_ys):
        is_data = i == 3
        fc = P["lit"] if not is_data else P["data"]
        ec = "white"   if not is_data else P["data_border"]
        tc = "white"   if not is_data else P["title"]
        box = FancyBboxPatch((LIT_X - LIT_W / 2, ly - 0.48), LIT_W, 0.95,
                              boxstyle="round,pad=0.08",
                              facecolor=fc, edgecolor=ec, linewidth=1.2, zorder=3)
        ax.add_patch(box)
        ax.text(LIT_X, ly, label, ha="center", va="center",
                fontsize=7.5, color=tc, zorder=4, multialignment="center",
                style="italic" if is_data else "normal")
        if i < len(lit_ys) - 1:
            arrow(ax, LIT_X, ly - 0.5, LIT_X, lit_ys[i + 1][0] + 0.5,
                  color="#808B96", lw=1.2, head=0.13)

    # Arrow from PubMed DB down to first lit box
    arrow(ax, LIT_X, STAGES_Y["db"] - 0.45,
          LIT_X, lit_ys[0][0] + 0.5, color="#808B96", lw=1.3, head=0.13)

    # Dashed arrow from lit corpus to stage 05/06 (enriches biological context)
    dashed_arrow(ax,
                 LIT_X - LIT_W / 2, lit_ys[3][0],
                 CX + BOX_W / 2 + 0.1, STAGES_Y["s06"],
                 color="#AAB7B8", lw=1.1)
    ax.text((LIT_X - LIT_W / 2 + CX + BOX_W / 2) / 2, (lit_ys[3][0] + STAGES_Y["s06"]) / 2 + 0.2,
            "context", fontsize=7, color="#AAB7B8", ha="center", style="italic")

    # Literature module label
    ax.text(LIT_X, 21.6, "Literature Mining\nModule",
            ha="center", va="center", fontsize=9, fontweight="bold",
            color=P["lit"], multialignment="center")
    ax.annotate("", xy=(LIT_X, 21.3), xytext=(LIT_X, 21.25),
                arrowprops=dict(arrowstyle="-", color=P["lit"], lw=1))

    # Bracket around literature module
    bx = LIT_X - LIT_W / 2 - 0.15
    ax.plot([bx, bx], [15.1, 21.8], color="#BDC3C7", lw=1.2, zorder=1)
    ax.plot([bx, bx + 0.15], [15.1, 15.1], color="#BDC3C7", lw=1.2, zorder=1)
    ax.plot([bx, bx + 0.15], [21.8, 21.8], color="#BDC3C7", lw=1.2, zorder=1)

    # ── Legend ────────────────────────────────────────────────────────────────
    leg_items = [
        mpatches.Patch(facecolor=P["s00"], label="Data Acquisition (00–01)"),
        mpatches.Patch(facecolor=P["s02"], label="Quality Control (02)"),
        mpatches.Patch(facecolor=P["s03"], label="Quantification & Processing (03–04)"),
        mpatches.Patch(facecolor=P["s05"], label="Statistical Analysis (05)"),
        mpatches.Patch(facecolor=P["s06"], label="Biological Interpretation (06–07)"),
        mpatches.Patch(facecolor=P["lit"], label="Literature Mining (parallel)"),
        mpatches.Patch(facecolor=P["data"], edgecolor=P["data_border"],
                        label="Intermediate data / files"),
        Line2D([0], [0], color=P["arrow"], linewidth=1.6,
               label="Data flow"),
        Line2D([0], [0], color="#AAB7B8", linewidth=1.2, linestyle="dashed",
               label="Supplementary link"),
    ]
    ax.legend(handles=leg_items, loc="lower left",
              bbox_to_anchor=(0.01, 0.01),
              fontsize=7.5, ncol=2, framealpha=0.92,
              edgecolor="#BDC3C7", facecolor="white",
              title="Legend", title_fontsize=8)

    # ── Outer border ──────────────────────────────────────────────────────────
    border = FancyBboxPatch((0.1, 0.1), FW - 0.2, FH - 0.2,
                             boxstyle="round,pad=0.05",
                             facecolor="none", edgecolor="#BDC3C7",
                             linewidth=1.5, zorder=0)
    ax.add_patch(border)

    return fig


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    out_dir = os.path.abspath(sys.argv[1] if len(sys.argv) > 1
                              else os.path.dirname(__file__))
    os.makedirs(out_dir, exist_ok=True)

    print("Drawing pipeline diagram …")
    fig = draw()

    png_path = os.path.join(out_dir, "pipeline_diagram.png")
    pdf_path = os.path.join(out_dir, "pipeline_diagram.pdf")

    fig.savefig(png_path, dpi=300, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print(f"PNG (300 DPI) : {png_path}")

    fig.savefig(pdf_path, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print(f"PDF (vector)  : {pdf_path}")

    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
