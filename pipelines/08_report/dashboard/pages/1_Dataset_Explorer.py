"""Per-dataset explorer: QC, differential expression, enrichment, IBD overlap."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import streamlit.components.v1 as components

import data as D

st.set_page_config(page_title="Dataset Explorer", page_icon="📊", layout="wide")
st.title("📊 Dataset Explorer")

analyzed = D.analyzed()
if not analyzed:
    st.warning("No analyzed datasets found in the manifest.")
    st.stop()

ds_id = st.selectbox("Dataset", [d["id"] for d in analyzed])
d = D.get_dataset(ds_id)
odir = D.resolve(d["output_dir"])

# ── Header ────────────────────────────────────────────────────────────────────
st.markdown(f"### {d.get('title') or ds_id}")
hc = st.columns(4)
hc[0].metric("Samples", d.get("n_samples_analyzed") or "—")
hc[1].metric("Tissue", d.get("tissue") or "—")
hc[2].metric("Strategy", d.get("library_strategy") or "—")
hc[3].metric("Comparisons", len(d.get("comparisons") or []))
if d.get("url"):
    st.markdown(f"[NCBI record]({d['url']})" + (f"  ·  PMID {d['pmid']}" if d.get("pmid") else ""))

tab_qc, tab_dea, tab_enrich, tab_ibd, tab_reports = st.tabs(
    ["Samples & QC", "Differential expression", "Enrichment", "IBD overlap", "Reports"]
)

# ── Samples & QC ──────────────────────────────────────────────────────────────
with tab_qc:
    meta = D.sample_metadata(d)
    if meta is not None and "disease" in meta.columns:
        st.markdown("**Sample composition**")
        counts = meta["disease"].value_counts().reset_index()
        counts.columns = ["disease", "n"]
        fig = px.bar(counts, x="disease", y="n", color="disease", text="n")
        fig.update_layout(showlegend=False, height=320, margin=dict(t=10, b=10))
        st.plotly_chart(fig, use_container_width=True)
        with st.expander("Sample metadata table"):
            st.dataframe(meta, use_container_width=True, hide_index=True)
    else:
        st.info("No sample metadata CSV registered for this dataset.")

    figs = D.qc_figures(odir)
    if figs:
        st.markdown("**QC figures**")
        fcols = st.columns(len(figs))
        for col, (label, path) in zip(fcols, figs.items()):
            col.image(str(path), caption=label, use_container_width=True)
    else:
        st.info("No QC figures found.")

# ── Differential expression ───────────────────────────────────────────────────
with tab_dea:
    dea = D.discover_dea(odir)
    if not dea:
        st.info("No DEA result CSVs found.")
    else:
        cc = st.columns([1, 1, 1, 1])
        mode = cc[0].selectbox("Replacement mode", list(dea.keys()),
                               format_func=lambda m: D._MODE_LABELS.get(m, m))
        comps = list(dea[mode].keys())
        comparison = cc[1].selectbox("Comparison", comps)
        padj_cut = cc[2].slider("padj ≤", 0.0, 0.25, 0.05, 0.005)
        lfc_cut = cc[3].slider("|log2FC| ≥", 0.0, 3.0, 1.0, 0.25)

        df = D.load_dea(dea[mode][comparison])
        df = df[df["padj"].notna()].copy()
        df["sig"] = (df["padj"] <= padj_cut) & (df["log2FoldChange"].abs() >= lfc_cut)
        df["dir"] = np.where(~df["sig"], "ns",
                             np.where(df["log2FoldChange"] >= 0, "UP", "DOWN"))
        df["neglog10padj"] = -np.log10(df["padj"].clip(lower=df["padj"][df["padj"] > 0].min() or 1e-300))

        n_up = int((df["dir"] == "UP").sum())
        n_down = int((df["dir"] == "DOWN").sum())
        m1, m2, m3 = st.columns(3)
        m1.metric("Significant", n_up + n_down)
        m2.metric("Up", n_up)
        m3.metric("Down", n_down)

        fig = px.scatter(
            df, x="log2FoldChange", y="neglog10padj", color="dir",
            color_discrete_map={"UP": "#d62728", "DOWN": "#1f77b4", "ns": "#cccccc"},
            hover_name="miRNA", opacity=0.75,
            labels={"neglog10padj": "-log10(padj)", "dir": ""},
        )
        fig.add_vline(x=lfc_cut, line_dash="dot", line_color="grey")
        fig.add_vline(x=-lfc_cut, line_dash="dot", line_color="grey")
        fig.add_hline(y=-np.log10(padj_cut), line_dash="dot", line_color="grey")
        fig.update_layout(height=480, margin=dict(t=10, b=10),
                          title=f"Volcano — {comparison} ({mode})")
        st.plotly_chart(fig, use_container_width=True)

        sig_tbl = (df[df["sig"]]
                   .sort_values("padj")
                   [["miRNA", "log2FoldChange", "padj", "baseMean", "dir"]]
                   .reset_index(drop=True))
        st.markdown(f"**Significant miRNAs ({len(sig_tbl)})**")
        st.dataframe(sig_tbl, use_container_width=True, hide_index=True)
        st.download_button(
            "⬇️ Download significant table (CSV)",
            sig_tbl.to_csv(index=False).encode(),
            file_name=f"{ds_id}_{mode}_{comparison}_DE.csv",
            mime="text/csv",
        )

# ── Enrichment ────────────────────────────────────────────────────────────────
with tab_enrich:
    for kind in ("KEGG", "GO"):
        st.markdown(f"#### {kind}")
        dot = D.pathway_dotplot(odir, kind)
        tbl = D.pathway_table(odir, kind)
        ecols = st.columns([1, 1])
        if dot:
            ecols[0].image(str(dot), caption=f"{kind} dotplot", use_container_width=True)
        else:
            ecols[0].info(f"No {kind} dotplot.")
        if tbl is not None:
            keep = [c for c in ["ONTOLOGY", "category", "ID", "Description",
                                "FoldEnrichment", "p.adjust", "Count"] if c in tbl.columns]
            view = tbl[keep] if keep else tbl
            ecols[1].dataframe(view.head(50), use_container_width=True, hide_index=True)
        else:
            ecols[1].info(f"No {kind} enrichment table.")
        st.divider()

# ── IBD overlap ───────────────────────────────────────────────────────────────
with tab_ibd:
    ibd = D.ibd_overlap_table(odir)
    if ibd is None:
        st.info("No IBD target-overlap table found.")
    else:
        show = [c for c in ["miRNA", "n_validated", "n_ibd_overlap", "pct_ibd",
                            "OR", "p_fisher", "padj"] if c in ibd.columns]
        sort_col = "padj" if "padj" in ibd.columns else show[-1]
        st.caption("miRNAs whose validated targets are enriched for known IBD genes "
                   "(Fisher's exact test).")
        st.dataframe(ibd[show].sort_values(sort_col).reset_index(drop=True),
                     use_container_width=True, hide_index=True)

# ── Reports ───────────────────────────────────────────────────────────────────
with tab_reports:
    report = D.resolve(d.get("final_report"))
    if report and report.exists():
        html = report.read_text(encoding="utf-8", errors="ignore")
        st.download_button("⬇️ Download Final_Report.html", html.encode(),
                           file_name=f"{ds_id}_Final_Report.html", mime="text/html")
        with st.expander("Preview Final_Report.html inline", expanded=False):
            components.html(html, height=900, scrolling=True)
    else:
        st.info("No Final_Report.html for this dataset.")

    nets = D.network_html(odir)
    if nets:
        st.markdown("**Interactive network dashboards**")
        for kind, path in nets.items():
            st.download_button(f"⬇️ {kind} network dashboard",
                               path.read_bytes(),
                               file_name=f"{ds_id}_network_{kind}.html",
                               mime="text/html", key=f"net_{kind}")
