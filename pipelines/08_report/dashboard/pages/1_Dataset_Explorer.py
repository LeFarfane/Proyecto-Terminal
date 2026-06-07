"""Per-dataset explorer: QC, differential expression, enrichment, IBD overlap."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

import data as D

st.set_page_config(page_title="Dataset Explorer", page_icon="📊", layout="wide")
st.markdown(D.LOADING_OVERLAY, unsafe_allow_html=True)
st.title("📊 Dataset Explorer")

analyzed = D.analyzed()
if not analyzed:
    st.warning("No analyzed datasets found in the manifest.")
    st.stop()

ds_id = st.selectbox("Dataset", [d["id"] for d in analyzed], format_func=D.ds_label)
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

# The Clinical tab only exists for datasets with treatment-response (responder vs
# non-responder) data, so it's inserted conditionally before Reports.
clin = D.clinical_results(odir)
tab_labels = ["Samples & QC", "Differential expression", "Enrichment", "IBD overlap"]
if clin:
    tab_labels.append("Clinical (treatment response)")
tab_labels.append("Reports")
tabs = st.tabs(tab_labels)
tab_qc, tab_dea, tab_enrich, tab_ibd = tabs[0], tabs[1], tabs[2], tabs[3]
tab_clin = tabs[4] if clin else None
tab_reports = tabs[-1]

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
            D.df_with_copy(meta, use_container_width=True, hide_index=True)
            st.download_button("⬇️ Download sample metadata (CSV)",
                               meta.to_csv(index=False).encode(),
                               file_name=f"{ds_id}_sample_metadata.csv", mime="text/csv",
                               key=f"dl_meta_{ds_id}")
    else:
        st.info("No sample metadata CSV registered for this dataset.")

    figs = D.qc_figures(odir)
    if figs:
        st.markdown("**QC figures**")
        fcols = st.columns(len(figs))
        for col, (label, path) in zip(fcols, figs.items()):
            with col:
                D.zoomable_image(path, caption=label)
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
        # Step 0.01 so 0.58 (≈1.5-fold, the project's biologically meaningful
        # minimum used on the cards / cross-dataset page) is selectable.
        lfc_cut = cc[3].slider("|log2FC| ≥", 0.0, 3.0, 0.58, 0.01)

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
        # Roomy top margin so the title clears the metric row above it.
        fig.update_layout(height=480, margin=dict(t=60, b=10),
                          title=dict(text=f"Volcano — {comparison} ({mode})", y=0.97))
        st.plotly_chart(fig, use_container_width=True)

        sig_tbl = (df[df["sig"]]
                   .sort_values("padj")
                   [["miRNA", "log2FoldChange", "padj", "baseMean", "dir"]]
                   .reset_index(drop=True))
        st.markdown(f"**Significant miRNAs ({len(sig_tbl)})**")
        D.df_with_copy(sig_tbl, use_container_width=True, hide_index=True,
                     column_config=D.sci_format("padj"))
        st.download_button(
            "⬇️ Download significant table (CSV)",
            sig_tbl.to_csv(index=False).encode(),
            file_name=f"{ds_id}_{mode}_{comparison}_DE.csv",
            mime="text/csv",
            key=f"dl_dea_{ds_id}_{mode}_{comparison}",
        )

# ── Enrichment ────────────────────────────────────────────────────────────────
with tab_enrich:
    for kind in ("KEGG", "GO"):
        st.markdown(f"#### {kind}")
        dot = D.pathway_dotplot(odir, kind)
        tbl = D.pathway_table(odir, kind)
        ecols = st.columns([1, 1])
        if dot:
            with ecols[0]:
                D.zoomable_image(dot, caption=f"{kind} dotplot")
        else:
            ecols[0].info(f"No {kind} dotplot.")
        if tbl is not None:
            keep = [c for c in ["ONTOLOGY", "category", "ID", "Description",
                                "FoldEnrichment", "p.adjust", "Count"] if c in tbl.columns]
            view = tbl[keep] if keep else tbl
            ecols[1].dataframe(view.head(50), use_container_width=True, hide_index=True)
            ecols[1].download_button(f"⬇️ Download {kind} table (CSV)",
                                     view.to_csv(index=False).encode(),
                                     file_name=f"{ds_id}_{kind}_enrichment.csv",
                                     mime="text/csv", key=f"dl_enrich_{ds_id}_{kind}")
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
                   "(Fisher's exact test). `padj` here is the Fisher test adjusted "
                   "p-value — unrelated to the DESeq2 `padj` in the DEA tab.")
        ibd_tbl = ibd[show].sort_values(sort_col).reset_index(drop=True)
        D.df_with_copy(ibd_tbl, use_container_width=True, hide_index=True,
                     column_config=D.sci_format("p_fisher", "padj"))
        st.download_button("⬇️ Download IBD overlap table (CSV)",
                           ibd_tbl.to_csv(index=False).encode(),
                           file_name=f"{ds_id}_ibd_target_overlap.csv", mime="text/csv",
                           key=f"dl_ibd_{ds_id}")

    # ── Validated vs Predicted vs Baseline ───────────────────────────────────
    cmp = D.ibd_overlap_comparison(odir)
    if cmp is None:
        st.divider()
        st.caption("ℹ️ *Validated vs Predicted vs Baseline* comparison not generated "
                   "for this dataset yet. Run `06_targets_enrich/15_multimir_targets_"
                   "baseline.R` then re-run `07_networks_clin/19_ibd_target_overlap.R`.")
    elif {"evidence", "is_significant", "OR", "pct_ibd"}.issubset(cmp.columns):
        st.divider()
        st.markdown("#### Validated vs Predicted vs Baseline")
        st.caption(
            "Same Fisher's-exact IBD-overlap test, three layers. **Validated** & "
            "**Predicted** are the DE-significant miRNAs tested on their validated vs "
            "computationally-predicted (non-validated) targets. **Baseline** is *every* "
            "miRNA tested in the DEA (no cutoff) — the reference for what IBD overlap "
            "looks like for an arbitrary miRNA. Each evidence type uses its own pooled "
            "target universe, so OR>1 means 'more IBD-enriched than that targetome's "
            "average'.")

        cmp = cmp.copy()
        cmp["is_significant"] = cmp["is_significant"].astype(bool)
        sig = cmp[cmp["is_significant"]]
        metric = st.radio("Compare by", ["OR", "pct_ibd"], horizontal=True,
                          format_func=lambda m: {"OR": "Enrichment (odds ratio)",
                                                 "pct_ibd": "% targets that are IBD genes"}[m],
                          key=f"ibdcmp_metric_{ds_id}")

        # Per-miRNA grouped bars: validated vs predicted for the significant miRNAs,
        # with the baseline median of each evidence drawn as a reference line.
        if sig.empty:
            st.info("No significant miRNAs are present in the comparison table.")
        else:
            base_med = (cmp.groupby("evidence")[metric].median())
            fig = px.bar(sig.sort_values("miRNA"), x="miRNA", y=metric, color="evidence",
                         barmode="group",
                         color_discrete_map={"validated": "#1f77b4", "predicted": "#ff7f0e"},
                         labels={metric: ("Odds ratio" if metric == "OR" else "% IBD targets"),
                                 "evidence": "Evidence"})
            for ev, col in (("validated", "#1f77b4"), ("predicted", "#ff7f0e")):
                if ev in base_med.index:
                    fig.add_hline(y=float(base_med[ev]), line_dash="dot", line_color=col,
                                  annotation_text=f"{ev} baseline median",
                                  annotation_font_color=col, annotation_position="top left")
            if metric == "OR":
                fig.add_hline(y=1.0, line_color="#bbb")
            fig.update_layout(height=420, margin=dict(t=30, b=10))
            st.plotly_chart(fig, use_container_width=True)
            st.caption("**How to read:** bars are the significant miRNAs (blue = validated "
                       "targets, orange = predicted). Dotted lines = the baseline (all-miRNA) "
                       "median for each evidence type — bars above their line are more "
                       "IBD-enriched than a typical miRNA in this dataset.")

            # Side-by-side significant table (validated vs predicted per miRNA).
            wide = sig.pivot_table(index="miRNA", columns="evidence",
                                   values=["pct_ibd", "OR", "padj"], aggfunc="first")
            wide.columns = [f"{a}_{b}" for a, b in wide.columns]
            st.markdown("**Significant miRNAs — validated vs predicted**")
            D.df_with_copy(wide.reset_index(), use_container_width=True, hide_index=True,
                         column_config=D.sci_format("padj_validated", "padj_predicted"))

        st.download_button("⬇️ Download comparison table (CSV)",
                           cmp.to_csv(index=False).encode(),
                           file_name=f"{ds_id}_ibd_overlap_comparison.csv", mime="text/csv",
                           key=f"dl_ibdcmp_{ds_id}")

# ── Clinical (treatment response) — only for responder vs non-responder data ──
if clin:
    with tab_clin:
        st.markdown("#### Treatment-response biomarkers")
        st.caption(
            "Only for datasets with responder vs non-responder samples at baseline "
            "(here: glucocorticoid response in UC). Per-miRNA ROC/AUC from "
            "`07_networks_clin/17_treatment_response_roc.R` on the VST matrix, for two "
            "contrasts: **R vs NR** (can baseline expression predict response?) and "
            "**UC vs HC** (disease vs control). Candidate miRNAs are the union of the "
            "Responder & NonResponder DE-significant sets (padj ≤ 0.05, |log2FC| ≥ 0.58).")

        roc_rn = clin.get("roc_rvsnr")
        roc_uh = clin.get("roc_uchc")
        cv = clin.get("cv_model")
        cand = clin.get("candidates")

        mc = st.columns(4)
        mc[0].metric("Top AUC · R vs NR",
                     f"{roc_rn['AUC'].max():.3f}" if roc_rn is not None and not roc_rn.empty else "—")
        mc[1].metric("Top AUC · UC vs HC",
                     f"{roc_uh['AUC'].max():.3f}" if roc_uh is not None and not roc_uh.empty else "—")
        mc[2].metric("Multi-miRNA CV-AUC",
                     f"{cv['AUC_CV'].iloc[0]:.3f}" if cv is not None and not cv.empty else "—")
        mc[3].metric("Candidate miRNAs", len(cand) if cand is not None else "—")

        rc = st.columns(2)
        with rc[0]:
            st.markdown("**R vs NR — per-miRNA AUC**")
            if roc_rn is not None and not roc_rn.empty:
                t = roc_rn.sort_values("AUC", ascending=False).reset_index(drop=True)
                D.df_with_copy(t, use_container_width=True, hide_index=True,
                             column_config={"AUC": st.column_config.NumberColumn(format="%.3f")})
                st.download_button("⬇️ Download R-vs-NR ROC (CSV)",
                                   t.to_csv(index=False).encode(),
                                   file_name=f"{ds_id}_roc_R_vs_NR.csv", mime="text/csv",
                                   key=f"dl_roc_rn_{ds_id}")
            else:
                st.info("No R-vs-NR ROC table.")
        with rc[1]:
            st.markdown("**UC vs HC — per-miRNA AUC**")
            if roc_uh is not None and not roc_uh.empty:
                t = roc_uh.sort_values("AUC", ascending=False).reset_index(drop=True)
                D.df_with_copy(t, use_container_width=True, hide_index=True,
                             column_config={"AUC": st.column_config.NumberColumn(format="%.3f")})
                st.download_button("⬇️ Download UC-vs-HC ROC (CSV)",
                                   t.to_csv(index=False).encode(),
                                   file_name=f"{ds_id}_roc_UC_vs_HC.csv", mime="text/csv",
                                   key=f"dl_roc_uh_{ds_id}")
            else:
                st.info("No UC-vs-HC ROC table.")

        if cv is not None and not cv.empty:
            st.caption(f"**Multi-miRNA CV model (R vs NR):** {cv['model'].iloc[0]} · "
                       f"AUC_CV = {cv['AUC_CV'].iloc[0]:.3f} · n = {int(cv['n'].iloc[0])} · "
                       f"miRNAs: {cv['miRNAs'].iloc[0]}")

        if cand is not None and not cand.empty:
            with st.expander(f"Candidate miRNAs ({len(cand)})"):
                D.df_with_copy(cand, use_container_width=True, hide_index=True,
                             column_config=D.sci_format("padj"))
                st.download_button("⬇️ Download candidates (CSV)",
                                   cand.to_csv(index=False).encode(),
                                   file_name=f"{ds_id}_candidate_mirnas.csv", mime="text/csv",
                                   key=f"dl_cand_{ds_id}")

        st.divider()
        multi = clin.get("roc_multi")
        if multi:
            st.markdown("**Multi-miRNA CV ROC**")
            mcol, _ = st.columns([1, 2])
            with mcol:
                D.zoomable_image(multi)

        def _gallery(title, pngs):
            if not pngs:
                return
            with st.expander(f"{title} ({len(pngs)})"):
                cols = st.columns(3)
                for i, p in enumerate(pngs):
                    with cols[i % 3]:
                        D.zoomable_image(p, caption=p.stem)

        _gallery("ROC curves — R vs NR", clin.get("roc_rvsnr_pngs", []))
        _gallery("ROC curves — UC vs HC", clin.get("roc_uchc_pngs", []))
        _gallery("Violin plots (HC / NR / R)", clin.get("violin_pngs", []))

# ── Reports ───────────────────────────────────────────────────────────────────
with tab_reports:
    # Each report: open it in the OS browser (left) or download it (right).
    reports = []  # (label, Path, download_filename)
    final_report = D.resolve(d.get("final_report"))
    if final_report and final_report.exists():
        reports.append(("📄 Final Report", final_report, f"{ds_id}_Final_Report.html"))
    for label, qc_path in D.find_qc_reports(odir).items():
        icon = "🔬" if label == "MultiQC" else "🧬"
        reports.append((f"{icon} {label}", qc_path, f"{ds_id}_{label}.html"))

    if not reports:
        st.info("No reports found for this dataset.")
    else:
        st.caption("Open a report in your browser, or download a copy.")
        for i, (label, path, fname) in enumerate(reports):
            nc, oc, dc = st.columns([2, 1, 1])
            nc.markdown(f"**{label}**")
            if oc.button("🌐 Open", key=f"open_rep_{i}", use_container_width=True):
                D.open_in_os(path)
            dc.download_button("⬇️ Download", path.read_bytes(), file_name=fname,
                               mime="text/html", key=f"dl_rep_{i}",
                               use_container_width=True)

    nets = D.network_html(odir)
    if nets:
        st.markdown("**Interactive network dashboards**")
        for kind, path in nets.items():
            nc, oc, dc = st.columns([2, 1, 1])
            nc.markdown(f"**🕸️ {kind} network**")
            if oc.button("🌐 Open", key=f"open_net_{kind}", use_container_width=True):
                D.open_in_os(path)
            dc.download_button("⬇️ Download", path.read_bytes(),
                               file_name=f"{ds_id}_network_{kind}.html",
                               mime="text/html", key=f"dl_net_{kind}",
                               use_container_width=True)
