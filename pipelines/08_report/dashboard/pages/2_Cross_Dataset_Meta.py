"""Cross-dataset meta-analysis: overlap (Venn/UpSet), direction matrix,
formal Stouffer meta-analysis, effect-size concordance, and shared pathways."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from upsetplot import from_contents, UpSet
from scipy.stats import spearmanr
import streamlit as st

import data as D

st.set_page_config(page_title="Cross-Dataset Meta-Analysis", page_icon="🔬", layout="wide")
st.title("🔬 Cross-Dataset Meta-Analysis")
st.caption(
    "Find miRNA patterns that replicate across studies. Note: PRJNA717018 is "
    "peripheral **blood** while the others are **mucosa/intestinal** — overlap "
    "across tissues means something different from within-tissue overlap."
)

analyzed = D.analyzed()
if len(analyzed) < 2:
    st.warning("Need at least two analyzed datasets for a meta-analysis.")
    st.stop()


def comparison_selectors(mode, key):
    """Render one comparison selectbox per analyzed dataset; return {id: comparison}."""
    sel = {}
    cols = st.columns(min(len(analyzed), 3))
    for i, d in enumerate(analyzed):
        comps = D.list_comparisons(d, mode)
        with cols[i % len(cols)]:
            sel[d["id"]] = st.selectbox(d["id"], comps, key=f"{key}_{d['id']}") if comps else None
            if not comps:
                st.caption("no comparisons")
    return sel


# ── Shared threshold controls (Overlap + Direction tabs) ──────────────────────
cc = st.columns(3)
mode = cc[0].selectbox("Replacement mode", ["Default", "Strict"],
                       format_func=lambda m: D._MODE_LABELS.get(m, m))
padj_cut = cc[1].slider("padj ≤", 0.0, 0.25, 0.05, 0.005)
lfc_cut = cc[2].slider("|log2FC| ≥", 0.0, 3.0, 1.0, 0.25)

tab_overlap, tab_matrix, tab_meta, tab_concord, tab_path = st.tabs(
    ["Overlap (Venn / UpSet)", "Direction matrix", "Meta-analysis",
     "Concordance", "Shared pathways"]
)

# ══ OVERLAP: Venn + UpSet ═════════════════════════════════════════════════════
with tab_overlap:
    sets = D.de_sets(padj_cut, lfc_cut, mode=mode)
    sets = {k: v for k, v in sets.items() if v}
    if len(sets) < 2:
        st.info("Fewer than two datasets have DE miRNAs at these thresholds. Loosen the cutoffs.")
    else:
        st.markdown("#### Venn diagram")
        st.caption("Best for ≤3 sets — pick which datasets to compare.")
        chosen = st.multiselect("Datasets in the Venn (2–3)", list(sets.keys()),
                                default=list(sets.keys())[:3])
        if len(chosen) in (2, 3):
            fig, ax = plt.subplots(figsize=(6, 5))
            labels = [s.replace("_", " ") for s in chosen]
            data_sets = [sets[c] for c in chosen]
            if len(chosen) == 2:
                venn2(data_sets, set_labels=labels, ax=ax)
            else:
                venn3(data_sets, set_labels=labels, ax=ax)
            st.pyplot(fig)
            plt.close(fig)
            # shared list
            shared = set.intersection(*data_sets)
            st.caption(f"Shared by all {len(chosen)} selected: "
                       + (", ".join(sorted(shared)) if shared else "none"))
        else:
            st.info("Select 2 or 3 datasets for the Venn (use UpSet below for more).")

        st.divider()
        st.markdown("#### UpSet plot")
        st.caption("Scales to any number of datasets — each bar is one intersection.")
        try:
            contents = from_contents(sets)
            fig = plt.figure(figsize=(9, 4.5))
            UpSet(contents, show_counts=True, sort_by="cardinality").plot(fig=fig)
            st.pyplot(fig)
            plt.close(fig)
        except Exception as e:
            st.warning(f"UpSet could not be drawn: {e}")

# ══ DIRECTION MATRIX (signed heatmap + recurrence + consensus) ════════════════
with tab_matrix:
    de_long = D.collect_de_mirnas(padj_cut, lfc_cut, mode=mode)
    if de_long.empty:
        st.info("No DE miRNAs at these thresholds.")
    else:
        mat, recurrence = D.consensus_matrix(de_long)
        n_datasets = mat.shape[1]
        recurrent = recurrence[recurrence >= 2]

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Datasets", n_datasets)
        m2.metric("Unique DE miRNAs", int(mat.shape[0]))
        m3.metric("Recurrent (≥2)", int((recurrence >= 2).sum()))
        m4.metric("Core (all)", int((recurrence == n_datasets).sum()))

        dist = recurrence.value_counts().sort_index().reset_index()
        dist.columns = ["n_datasets", "n_miRNAs"]
        fig = px.bar(dist, x="n_datasets", y="n_miRNAs", text="n_miRNAs",
                     labels={"n_datasets": "# datasets a miRNA is DE in", "n_miRNAs": "# miRNAs"})
        fig.update_layout(height=280, margin=dict(t=10, b=10), showlegend=False)
        fig.update_xaxes(dtick=1)
        st.plotly_chart(fig, use_container_width=True)

        if recurrent.empty:
            st.info("No miRNA is DE in two or more datasets at these thresholds.")
        else:
            cap = st.slider("Max miRNAs in heatmap", 10, max(20, len(recurrent)),
                            min(60, len(recurrent)), 5)
            sub = mat.loc[recurrent.index[:cap]]
            fig = px.imshow(sub, color_continuous_scale=["#1f77b4", "#dddddd", "#d62728"],
                            zmin=-1, zmax=1, aspect="auto",
                            labels=dict(color="direction", x="dataset", y="miRNA"))
            fig.update_layout(height=max(360, 18 * len(sub)), margin=dict(t=10, b=10),
                              coloraxis_colorbar=dict(tickvals=[-1, 0, 1],
                                                      ticktext=["DOWN", "mixed", "UP"]))
            st.plotly_chart(fig, use_container_width=True)
            st.caption("Red = up, blue = down, grey = mixed across comparisons, white = not DE.")

# ══ META-ANALYSIS (Stouffer) ══════════════════════════════════════════════════
with tab_meta:
    st.markdown("#### Sample-size-weighted Stouffer meta-analysis")
    st.caption("Combines per-dataset effect sizes into one ranked consensus signature "
               "with a BH-adjusted combined p-value. Pick one comparison per dataset.")
    msel = comparison_selectors(mode, key="meta")
    mc = st.columns(2)
    min_ds = mc[0].slider("Require miRNA in ≥ N datasets", 2, len(analyzed), 2, key="meta_minds")
    sig_cut = mc[1].slider("Highlight combined padj ≤", 0.0, 0.25, 0.05, 0.005, key="meta_sig")

    meta = D.meta_stouffer(msel, mode=mode, min_datasets=min_ds)
    if meta.empty:
        st.info("No miRNA met the criteria. Lower the dataset requirement or check comparisons.")
    else:
        n_sig = int((meta["combined_padj"] <= sig_cut).sum())
        a, b, c = st.columns(3)
        a.metric("miRNAs tested", len(meta))
        b.metric(f"Consensus hits (padj≤{sig_cut})", n_sig)
        c.metric("Max datasets", int(meta["n_datasets"].max()))

        plot = meta.copy()
        plot["neglog10padj"] = -np.log10(plot["combined_padj"].clip(lower=1e-300))
        plot["sig"] = np.where(plot["combined_padj"] <= sig_cut, plot["direction"], "ns")
        fig = px.scatter(plot, x="combined_Z", y="neglog10padj", color="sig",
                         color_discrete_map={"UP": "#d62728", "DOWN": "#1f77b4", "ns": "#cccccc"},
                         size="n_datasets", hover_name="miRNA",
                         labels={"combined_Z": "Stouffer Z (signed)", "neglog10padj": "-log10(combined padj)", "sig": ""})
        fig.add_hline(y=-np.log10(sig_cut), line_dash="dot", line_color="grey")
        fig.update_layout(height=460, margin=dict(t=10, b=10))
        st.plotly_chart(fig, use_container_width=True)

        show = meta[["miRNA", "n_datasets", "direction", "combined_Z",
                     "mean_log2FC", "combined_p", "combined_padj", "datasets"]]
        st.dataframe(show, use_container_width=True, hide_index=True)
        st.download_button("⬇️ Download consensus signature (CSV)",
                           show.to_csv(index=False).encode(),
                           file_name="meta_analysis_consensus.csv", mime="text/csv")

# ══ CONCORDANCE (pairwise log2FC) ═════════════════════════════════════════════
with tab_concord:
    st.markdown("#### Effect-size concordance")
    st.caption("Do shared miRNAs agree in direction and magnitude across datasets? "
               "Points on the diagonal = reproducible signal.")
    csel = comparison_selectors(mode, key="conc")
    wide = D.effect_size_wide(csel, mode=mode)
    valid_cols = [c for c in wide.columns if csel.get(c)]
    wide = wide[valid_cols] if valid_cols else wide

    if wide.shape[1] < 2:
        st.info("Need at least two datasets with a selected comparison.")
    else:
        if wide.shape[1] > 2:
            st.markdown("**Pairwise concordance (Spearman ρ on shared miRNAs)**")
            cols = list(wide.columns)
            corr = pd.DataFrame(np.eye(len(cols)), index=cols, columns=cols)
            for i in range(len(cols)):
                for j in range(i + 1, len(cols)):
                    pair = wide[[cols[i], cols[j]]].dropna()
                    r = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1]).statistic if len(pair) > 2 else np.nan
                    corr.iloc[i, j] = corr.iloc[j, i] = r
            fig = px.imshow(corr, text_auto=".2f", color_continuous_scale="RdBu_r",
                            zmin=-1, zmax=1, aspect="auto")
            fig.update_layout(height=320, margin=dict(t=10, b=10))
            st.plotly_chart(fig, use_container_width=True)

        pcols = st.columns(2)
        dx = pcols[0].selectbox("X axis", list(wide.columns), index=0, key="conc_x")
        dy = pcols[1].selectbox("Y axis", list(wide.columns),
                                index=1 if wide.shape[1] > 1 else 0, key="conc_y")
        if dx != dy:
            pair = wide[[dx, dy]].dropna()
            pair = pair.reset_index().rename(columns={"index": "miRNA"})
            if len(pair) > 2:
                rho = spearmanr(pair[dx], pair[dy]).statistic
                lim = float(np.nanmax(np.abs(pair[[dx, dy]].values))) * 1.05
                fig = px.scatter(pair, x=dx, y=dy, hover_name="miRNA", opacity=0.7,
                                 labels={dx: f"log2FC — {dx}", dy: f"log2FC — {dy}"})
                fig.add_shape(type="line", x0=-lim, y0=-lim, x1=lim, y1=lim,
                              line=dict(dash="dot", color="grey"))
                fig.add_hline(y=0, line_color="#eee"); fig.add_vline(x=0, line_color="#eee")
                fig.update_layout(height=460, margin=dict(t=30, b=10),
                                  title=f"{len(pair)} shared miRNAs · Spearman ρ = {rho:.2f}")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("Too few shared miRNAs to compare.")

# ══ SHARED PATHWAYS ═══════════════════════════════════════════════════════════
with tab_path:
    st.caption("Even when miRNAs differ, do datasets converge on the same biology? "
               "Enriched terms (p.adjust ≤ 0.05) shared across datasets.")
    pc1, pc2 = st.columns(2)
    for col, kind in ((pc1, "KEGG"), (pc2, "GO")):
        with col:
            st.markdown(f"**{kind}**")
            pc = D.pathway_consensus(kind, padj=0.05, top=40)
            if pc.empty:
                st.info(f"No shared {kind} terms found.")
            else:
                st.dataframe(pc, use_container_width=True, hide_index=True)
