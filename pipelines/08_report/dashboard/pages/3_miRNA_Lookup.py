"""miRNA Lookup — per-miRNA profile across datasets (DE + IBD relevance)."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pandas as pd
import plotly.express as px
import streamlit as st

import data as D

st.set_page_config(page_title="miRNA Lookup", page_icon="🔎", layout="wide")
st.markdown(D.LOADING_OVERLAY, unsafe_allow_html=True)
st.title("🔎 miRNA Lookup")
st.caption("Pick a miRNA to see everything farmed across datasets — where it is "
           "differentially expressed and whether its targets are IBD-relevant. The "
           "list is every miRNA DE in ≥1 dataset (padj ≤ 0.05, |log2FC| ≥ 0.58).")

universe = list(D.de_universe())
if not universe:
    st.warning("No DE miRNAs found across the analyzed datasets.")
    st.stop()

mirna = st.selectbox("miRNA", universe, index=0,
                     help="Type to filter. List = miRNAs DE in ≥1 dataset.")

dea = D.mirna_dea_profile(mirna)
ibd = D.mirna_ibd_profile(mirna)
n_total = len(D.analyzed())

# ── Header / verdict ──────────────────────────────────────────────────────────
sig = dea[dea["significant"]] if not dea.empty else dea
n_sig_ds = sig["dataset"].nunique() if not sig.empty else 0
dirs = set(sig["direction"]) if not sig.empty else set()
if dirs == {"UP"}:
    verdict = "consistently **UP**"
elif dirs == {"DOWN"}:
    verdict = "consistently **DOWN**"
elif dirs:
    verdict = "**mixed** direction"
else:
    verdict = "not significant in any dataset"

st.markdown(f"## {mirna}")
st.markdown(f"DE in **{n_sig_ds}/{n_total}** datasets — {verdict}  ·  "
            f"[miRBase ↗](https://www.mirbase.org/results/?query={mirna})")
st.divider()

# ── Section 1: DE across datasets ─────────────────────────────────────────────
st.markdown("### Differential expression across datasets")
if dea.empty:
    st.info("This miRNA was not found in any DEA table.")
else:
    modes = sorted(dea["mode"].unique())
    hmode = st.radio("Replacement mode (heatmap)", modes, horizontal=True,
                     format_func=lambda m: D._MODE_LABELS.get(m, m))
    sub = dea[dea["mode"] == hmode]
    pivot = (sub.pivot_table(index="dataset", columns="comparison",
                             values="log2FoldChange", aggfunc="first")
             .rename(index=D.ds_label))
    fig = px.imshow(pivot, color_continuous_scale="RdBu_r",
                    color_continuous_midpoint=0, aspect="auto", text_auto=".2f",
                    labels=dict(color="log2FC", x="comparison", y="dataset"))
    fig.update_layout(height=max(260, 70 * len(pivot) + 80), margin=dict(t=10, b=10))
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Red = up in the comparison, blue = down, blank = not present. "
               "See the table for padj / significance.")

    disp = dea.copy()
    disp["Dataset"] = disp["dataset"].map(D.ds_label)
    disp = (disp[["Dataset", "mode", "comparison", "log2FoldChange", "padj",
                  "pvalue", "baseMean", "direction", "significant"]]
            .sort_values(["significant", "padj"], ascending=[False, True])
            .reset_index(drop=True))
    D.df_with_copy(disp, use_container_width=True, hide_index=True,
                 column_config=D.sci_format("padj", "pvalue"))
    st.download_button("⬇️ Download DE profile (CSV)",
                       dea.to_csv(index=False).encode(),
                       file_name=f"{mirna}_DE_profile.csv", mime="text/csv",
                       key="dl_mirna_dea")

st.divider()

# ── Section 2: IBD relevance + target genes ───────────────────────────────────
st.markdown("### IBD relevance & target genes")
st.caption("From each dataset's validated-target IBD enrichment (Fisher's exact "
           "test). `OR`/`padj` here describe target-set enrichment, not DE.")
if ibd.empty:
    st.info("No IBD target-overlap entry for this miRNA in any dataset.")
else:
    keep = [c for c in ["dataset", "tissue", "OR", "p_fisher", "padj",
                        "n_validated", "n_ibd_overlap", "pct_ibd"] if c in ibd.columns]
    show = ibd[keep].copy()
    show["dataset"] = show["dataset"].map(D.ds_label)
    show = show.rename(columns={"dataset": "Dataset"})
    D.df_with_copy(show, use_container_width=True, hide_index=True,
                 column_config=D.sci_format("p_fisher", "padj"))

    genes = set()
    if "overlap_genes" in ibd.columns:
        for g in ibd["overlap_genes"].dropna():
            genes.update(x for x in str(g).split(";") if x)
    if genes:
        st.markdown(f"**IBD genes targeted (union across datasets): {len(genes)}**")
        for _, row in ibd.iterrows():
            members = sorted(x for x in str(row.get("overlap_genes") or "").split(";") if x)
            with st.expander(f"{D.ds_label(row['dataset'])} — {row.get('n_ibd_overlap')} IBD genes"):
                st.write(", ".join(members) or "—")
        st.download_button("⬇️ Download IBD union gene list (CSV)",
                           pd.DataFrame(sorted(genes), columns=["ibd_gene"])
                           .to_csv(index=False).encode(),
                           file_name=f"{mirna}_IBD_genes.csv", mime="text/csv",
                           key="dl_mirna_ibd_genes")
    st.download_button("⬇️ Download IBD overlap profile (CSV)",
                       ibd.to_csv(index=False).encode(),
                       file_name=f"{mirna}_IBD_profile.csv", mime="text/csv",
                       key="dl_mirna_ibd")
