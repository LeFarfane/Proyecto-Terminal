"""
IBD miRNA-seq dashboard — Home / Overview.

Run from the repo root with:
    streamlit run pipelines/08_report/dashboard/Dashboard.py
(or use the 24_launch_dashboard.sh launcher).
"""
import streamlit as st

import data as D

# Review-layer tag colours (human curation, parsed from outputs/docs/*.md).
# Unknown tags fall back to slate-blue.
_TAG_COLORS = {
    "to_review": "#8a6d00",                              # amber
    "priority": "#8a1f1f", "high_priority": "#8a1f1f",   # red
    "validated": "#1f5a2e", "reviewed": "#1f5a2e", "done": "#1f5a2e",  # green
    "ignore": "#444", "exclude": "#444", "skip": "#444", # grey
}


def _tag_chips(tags: list[str]) -> str:
    """Render a list of tags as small coloured HTML chips."""
    spans = []
    for t in tags:
        c = _TAG_COLORS.get(t.lower(), "#34506b")
        spans.append(
            f"<span style='background:{c};color:#fff;border-radius:6px;"
            f"padding:1px 7px;margin:0 4px 2px 0;font-size:0.72rem;"
            f"display:inline-block;white-space:nowrap;'>#{t}</span>"
        )
    return "<div style='margin:2px 0 4px;'>" + "".join(spans) + "</div>"


st.set_page_config(
    page_title="IBD miRNA-seq Dashboard",
    page_icon="🧬",
    layout="wide",
)

# Blocking semi-transparent overlay while the app is loading/rerunning.
st.markdown(D.LOADING_OVERLAY, unsafe_allow_html=True)

st.title("🧬 IBD miRNA-seq — Cross-Dataset Dashboard")
st.caption(
    "Driven by `datasets.yaml`. Use the pages in the sidebar to explore a single "
    "dataset or compare differentially expressed miRNAs across all analyzed datasets."
)

# ── Manifest refresh control ──────────────────────────────────────────────────
# Make primary buttons blue (Streamlit's default primary is red). Several
# selectors cover testid changes across Streamlit versions.
st.markdown(
    "<style>"
    'button[kind="primary"],'
    '[data-testid="stBaseButton-primary"],'
    '[data-testid="baseButton-primary"]'
    "{background-color:#1f77b4 !important;border-color:#1f77b4 !important;color:#fff !important;}"
    'button[kind="primary"]:hover,'
    '[data-testid="stBaseButton-primary"]:hover,'
    '[data-testid="baseButton-primary"]:hover'
    "{background-color:#155a8a !important;border-color:#155a8a !important;}"
    "</style>",
    unsafe_allow_html=True,
)

with st.sidebar:
    st.markdown("### Manifest")
    add_new = st.checkbox(
        "Add new datasets", value=False,
        help="Also append skeleton entries for docs with no matching dataset "
             "(05_sync_yaml.py --add-new).",
    )
    if st.button("🔄 Update", type="primary", use_container_width=True,
                 help="Re-sync datasets.yaml from outputs/docs/*.md "
                      "(runs 05_sync_yaml.py in py_env)."):
        with st.spinner("Syncing datasets.yaml from docs…"):
            ok, out = D.run_sync(add_new=add_new)
        if ok:
            D.clear_caches()          # so the render below reads fresh values
            st.success("Manifest updated.")
        else:
            st.error("Sync failed.")
        with st.expander("Sync log", expanded=not ok):
            st.code(out or "(no output)")

    # Review-layer filter — tags come from the md docs, independent of status.
    _tag_universe = D.all_tags()
    sel_tags = st.multiselect(
        "Filter by tag", _tag_universe, default=[],
        help="Show only datasets whose doc carries any of these tags "
             "(parsed from outputs/docs/*.md).",
    ) if _tag_universe else []

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
D.df_with_copy(
    D.manifest_overview_df(),
    use_container_width=True,
    hide_index=True,
)

st.divider()

# ── Per-dataset cards ─────────────────────────────────────────────────────────
st.subheader("At a glance")

# Apply the sidebar tag filter (a dataset matches if it carries ANY selected tag).
cards_ds = all_ds
if sel_tags:
    _sel = {t.lower() for t in sel_tags}
    cards_ds = [d for d in all_ds if {t.lower() for t in d.get("tags", [])} & _sel]
    st.caption(f"Filtered to {len(cards_ds)} dataset(s) tagged {', '.join('#'+t for t in sel_tags)}.")
if not cards_ds:
    st.info("No datasets match the selected tags.")

# Render row by row (fresh st.columns per group) so every row aligns to a clean
# baseline and the vertical gaps between cards stay uniform. Cards keep their own
# natural height — only the spacing is regularized.
NCOLS = 3
for row_start in range(0, len(cards_ds), NCOLS):
    row = cards_ds[row_start:row_start + NCOLS]
    cols = st.columns(NCOLS, gap="large")
    for col, d in zip(cols, row):
        with col:
            with st.container(border=True):
                # Header: ID · status [· PMID link]
                header = f"**{d['id']}**  ·  {D.STATUS_LABEL.get(d.get('status'), '')}"
                if d.get("pmid"):
                    header += f"  ·  [PMID {d['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{d['pmid']}/)"
                st.markdown(header)

                # Status is disk-derived; flag when the manifest disagrees.
                if d.get("status_yaml") and d["status_yaml"] != d["status"]:
                    st.caption(
                        f"⚠️ manifest says _{D.STATUS_LABEL.get(d['status_yaml'], d['status_yaml'])}_"
                        " — card reflects disk"
                    )

                # Review-layer tags (from the md doc)
                if d.get("tags"):
                    st.markdown(_tag_chips(d["tags"]), unsafe_allow_html=True)

                # Study title
                st.caption(d.get("title") or "_(metadata not yet curated)_")

                # Tissue + sequencing strategy
                tech = []
                if d.get("tissue"):
                    tech.append(f"🧫 {d['tissue']}")
                if d.get("library_strategy"):
                    tech.append(f"📊 {d['library_strategy']}")
                if tech:
                    st.caption("  ·  ".join(tech))

                # Samples analysed vs total runs
                n_samp = d.get("n_samples_analyzed")
                n_runs = d.get("n_runs_total")
                if n_samp:
                    sample_str = f"👥 {n_samp} samples"
                    if n_runs and n_runs != n_samp:
                        sample_str += f" (of {n_runs} runs)"
                    st.caption(sample_str)

                # Comparisons
                if d.get("comparisons"):
                    st.caption(f"⚖️ {', '.join(d['comparisons'])}")

                # DE summary (analyzed datasets only)
                if d.get("status") == "analyzed":
                    stats = D.card_dea_stats(d["id"], lfc=0.58)
                    if stats and stats["n_sig"] > 0:
                        st.caption(f"🔬 {stats['n_sig']} DE miRNAs (padj<0.05, |LFC|≥0.58)")

                # File-open buttons (subprocess → explorer.exe, works in WSL)
                file_items = []
                final_report = D.resolve(d.get("final_report"))
                if final_report and final_report.exists():
                    file_items.append(("📄 Report", final_report))
                odir = D.resolve(d.get("output_dir"))
                if odir and odir.exists():
                    for label, qc_path in D.find_qc_reports(odir).items():
                        icon = "🔬" if label == "MultiQC" else "🧬"
                        file_items.append((f"{icon} {label}", qc_path))
                if file_items:
                    btn_cols = st.columns(len(file_items))
                    for j, (label, fpath) in enumerate(file_items):
                        with btn_cols[j]:
                            if st.button(label, key=f"open_{d['id']}_{j}",
                                         use_container_width=True):
                                D.open_in_os(fpath)

                # NCBI record as a secondary reference
                if d.get("url"):
                    st.caption(f"[↗ NCBI record]({d['url']})")

st.divider()
st.caption(
    f"Project: {proj.get('name', '—')}  ·  organism: {proj.get('organism', '—')}  ·  "
    f"manifest: `datasets.yaml`"
)
