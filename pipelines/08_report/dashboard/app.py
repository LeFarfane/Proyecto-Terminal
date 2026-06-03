"""
IBD miRNA-seq dashboard — Home / Overview.

Run from the repo root with:
    streamlit run pipelines/08_report/dashboard/app.py
(or use the run_dashboard.ps1 launcher).
"""
import streamlit as st

import data as D

st.set_page_config(
    page_title="IBD miRNA-seq Dashboard",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 IBD miRNA-seq — Cross-Dataset Dashboard")
st.caption(
    "Driven by `datasets.yaml`. Use the pages in the sidebar to explore a single "
    "dataset or compare differentially expressed miRNAs across all analyzed datasets."
)

mani = D.load_manifest()
proj = mani.get("project", {})
all_ds = D.datasets()
analyzed = D.analyzed()

# ── Headline metrics ──────────────────────────────────────────────────────────
n_analyzed = len(analyzed)
n_inprog = len(D.datasets(status="in_progress"))
n_susp = len(D.datasets(status="suspended"))
n_samples = sum((d.get("n_samples_analyzed") or 0) for d in analyzed)

c1, c2, c3, c4 = st.columns(4)
c1.metric("Datasets", len(all_ds))
c2.metric("Analyzed", n_analyzed)
c3.metric("In progress", n_inprog)
c4.metric("Samples analyzed", n_samples)

st.divider()

# ── Registry table ────────────────────────────────────────────────────────────
st.subheader("Dataset registry")
st.dataframe(
    D.manifest_overview_df(),
    use_container_width=True,
    hide_index=True,
)

st.divider()

# ── Per-dataset cards ─────────────────────────────────────────────────────────
st.subheader("At a glance")
cols = st.columns(3)
for i, d in enumerate(all_ds):
    with cols[i % 3]:
        with st.container(border=True):
            st.markdown(f"**{d['id']}**  ·  {D.STATUS_LABEL.get(d.get('status'), '')}")
            title = d.get("title") or "_(metadata not yet curated)_"
            st.caption(title)
            meta = []
            if d.get("tissue"):
                meta.append(f"🧫 {d['tissue']}")
            if d.get("n_samples_analyzed"):
                meta.append(f"👥 {d['n_samples_analyzed']} samples")
            if d.get("comparisons"):
                meta.append(f"⚖️ {', '.join(d['comparisons'])}")
            if meta:
                st.markdown("  \n".join(meta))
            if d.get("url"):
                st.markdown(f"[NCBI record]({d['url']})")

st.divider()
st.caption(
    f"Project: {proj.get('name', '—')}  ·  organism: {proj.get('organism', '—')}  ·  "
    f"manifest: `datasets.yaml`"
)
