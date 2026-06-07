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
import upset_compat  # noqa: F401  patches upsetplot 0.9.0 for modern pandas/matplotlib

st.set_page_config(page_title="Cross-Dataset Meta-Analysis", page_icon="🔬", layout="wide")
st.markdown(D.LOADING_OVERLAY, unsafe_allow_html=True)
st.title("🔬 Cross-Dataset Meta-Analysis")
st.caption(
    "Find miRNA patterns that replicate across studies. **Which miRNAs each tab "
    "uses:** *Overlap* and *Direction matrix* use each dataset's **significant DE "
    "miRNAs** (DESeq2, at the thresholds below); *Meta-analysis* and *Concordance* "
    "use the **full per-miRNA statistics** (no significance pre-filter); "
    "*IBD-target overlap* uses the **validated-target IBD-enrichment** results "
    "(its own controls). Note: PRJNA717018 is peripheral **blood** while the others "
    "are **mucosa/intestinal** — cross-tissue overlap means something different "
    "from within-tissue overlap."
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
            sel[d["id"]] = st.selectbox(D.ds_label(d["id"]), comps,
                                        key=f"{key}_{d['id']}") if comps else None
            if not comps:
                st.caption("no comparisons")
    return sel


# ── Shared threshold controls (Overlap + Direction tabs) ──────────────────────
cc = st.columns(3)
mode = cc[0].selectbox("Replacement mode", ["Default", "Strict"],
                       format_func=lambda m: D._MODE_LABELS.get(m, m))
padj_cut = cc[1].slider("padj ≤", 0.0, 0.25, 0.05, 0.005)
# Floor at 0.58 (≈1.5-fold) — the biologically meaningful minimum used on the cards.
lfc_cut = cc[2].slider("|log2FC| ≥", 0.58, 3.0, 0.58, 0.1)

(tab_overlap, tab_matrix, tab_meta, tab_concord, tab_path, tab_ibd,
 tab_cmp) = st.tabs(
    ["Overlap (Venn / UpSet)", "Direction matrix", "Meta-analysis",
     "Concordance", "Shared pathways", "IBD-target overlap",
     "Validated vs Predicted vs Baseline"]
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
                                default=list(sets.keys())[:3], format_func=D.ds_label)
        if len(chosen) in (2, 3):
            fig, ax = plt.subplots(figsize=(5, 4.5))
            labels = [D.ds_label(s) for s in chosen]
            data_sets = [sets[c] for c in chosen]
            v = (venn2(data_sets, set_labels=labels, ax=ax) if len(chosen) == 2
                 else venn3(data_sets, set_labels=labels, ax=ax))
            # Shrink the dataset name labels; the subset counts keep their size.
            for t in (v.set_labels or []):
                if t:
                    t.set_fontsize(8)
            # Keep the native figure size so it doesn't stretch to fill the page.
            D.zoomable_pyplot(fig)
            plt.close(fig)
            st.caption("**How to read:** each circle is one dataset's DE miRNAs; "
                       "overlaps are miRNAs DE in several datasets. The numbers count "
                       "how many miRNAs fall in each region (listed in the table below).")

            # Which miRNAs fall in each Venn region.
            named = {D.ds_label(c): sets[c] for c in chosen}
            regions = D.overlap_regions(named)
            reg_df = pd.DataFrame(
                [{"Region": r, "Count": len(members), "miRNAs": ", ".join(members)}
                 for r, members in regions.items()]
            )
            st.markdown("**miRNAs by region**")
            D.df_with_copy(reg_df, use_container_width=True, hide_index=True)
            tidy = pd.DataFrame(
                [(r, m) for r, members in regions.items() for m in members],
                columns=["Region", "miRNA"],
            )
            st.download_button("⬇️ Download Venn regions (CSV)",
                               tidy.to_csv(index=False).encode(),
                               file_name="venn_regions.csv", mime="text/csv",
                               key="dl_venn_regions")
        else:
            st.info("Select 2 or 3 datasets for the Venn (use UpSet below for more).")

        st.divider()
        st.markdown("#### UpSet plot")
        st.caption("Scales to any number of datasets — each bar is one intersection.")
        named_up = {D.ds_label(k): v for k, v in sets.items()}
        try:
            contents = from_contents(named_up)
            fig = plt.figure(figsize=(4, 3.2))
            UpSet(contents, show_counts=True, sort_by="cardinality").plot(fig=fig)
            D.zoomable_pyplot(fig)
            plt.close(fig)
            st.caption("**How to read:** each column is one exact dataset combination "
                       "(the dots show which datasets); bar height = how many miRNAs are "
                       "DE in exactly those datasets and no others. The Venn's scalable "
                       "sibling for >3 datasets.")
        except Exception as e:
            st.warning(f"UpSet could not be drawn: {e}")

        # Which miRNAs make up each UpSet bar (exclusive intersections).
        inter = D.exclusive_intersections(named_up)
        up_df = pd.DataFrame(
            [{"Datasets": lbl, "Count": len(members), "miRNAs": ", ".join(members)}
             for lbl, members in inter.items()]
        )
        st.markdown("**miRNAs by intersection**")
        D.df_with_copy(up_df, use_container_width=True, hide_index=True)
        up_tidy = pd.DataFrame(
            [(lbl, m) for lbl, members in inter.items() for m in members],
            columns=["Datasets", "miRNA"],
        )
        st.download_button("⬇️ Download UpSet intersections (CSV)",
                           up_tidy.to_csv(index=False).encode(),
                           file_name="upset_intersections.csv", mime="text/csv",
                           key="dl_upset")

# ══ DIRECTION MATRIX (signed heatmap + recurrence + consensus) ════════════════
with tab_matrix:
    de_long = D.collect_de_mirnas(padj_cut, lfc_cut, mode=mode)
    if de_long.empty:
        st.info("No DE miRNAs at these thresholds.")
    else:
        mat, recurrence = D.consensus_matrix(de_long)
        n_datasets = mat.shape[1]

        # Controls mirror the Meta-analysis tab: require a miRNA in ≥ N datasets
        # and highlight those that also pass a combined (Stouffer) padj — so you
        # see the per-dataset DIRECTION *and* the formal cross-dataset significance.
        ctrl = st.columns(2)
        min_ds = ctrl[0].slider("Require miRNA in ≥ N datasets", 2, max(2, n_datasets),
                                2, key="matrix_minds")
        sig_cut = ctrl[1].slider("Highlight combined padj ≤", 0.0, 0.25, 0.05, 0.005,
                                 key="matrix_sig")
        st.caption("**Combined padj** is the Stouffer meta-analysis FDR — pick one "
                   "comparison per dataset to define it. miRNAs passing it are marked "
                   "with ★ in the heatmap and listed first.")
        msel = comparison_selectors(mode, key="matrix")
        meta = D.meta_stouffer(msel, mode=mode, min_datasets=2)
        padj_map = (meta.set_index("miRNA")["combined_padj"].to_dict()
                    if not meta.empty else {})
        cdir_map = (meta.set_index("miRNA")["direction"].to_dict()
                    if not meta.empty else {})

        recurrent = recurrence[recurrence >= min_ds]

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Datasets", n_datasets)
        m2.metric(f"In ≥{min_ds} datasets", int((recurrence >= min_ds).sum()))
        m3.metric("Core (all)", int((recurrence == n_datasets).sum()))
        n_core_sig = sum(1 for m in recurrent.index
                         if padj_map.get(m, 1.0) <= sig_cut)
        m4.metric(f"≥{min_ds} & padj≤{sig_cut:g}", n_core_sig)

        dist = recurrence.value_counts().sort_index().reset_index()
        dist.columns = ["n_datasets", "n_miRNAs"]
        fig = px.bar(dist, x="n_datasets", y="n_miRNAs", text="n_miRNAs",
                     labels={"n_datasets": "# datasets a miRNA is DE in", "n_miRNAs": "# miRNAs"})
        fig.update_layout(height=280, margin=dict(t=10, b=10), showlegend=False)
        fig.update_xaxes(dtick=1)
        st.plotly_chart(fig, use_container_width=True)
        st.caption("**How to read:** how many miRNAs are DE in exactly 1, 2, 3… "
                   "datasets. Bars further right = more reproducible miRNAs (shared by "
                   "more studies); the leftmost bar is dataset-specific noise/biology.")

        if recurrent.empty:
            st.info(f"No miRNA is DE in ≥{min_ds} datasets at these thresholds.")
        else:
            only_sig = st.checkbox(
                f"Show only miRNAs with combined padj ≤ {sig_cut:g}",
                value=False, key="matrix_onlysig",
                disabled=not padj_map)

            idx = list(recurrent.index)
            if only_sig and padj_map:
                idx = [m for m in idx if padj_map.get(m, 1.0) <= sig_cut]
            if not idx:
                st.info("No miRNA meets the combined-padj cutoff. Loosen it, untick the "
                        "filter, or check your comparison selections.")
            else:
                # Order: most significant first (missing padj sinks to the bottom),
                # tie-broken by how many datasets the miRNA recurs in.
                idx = sorted(idx, key=lambda m: (padj_map.get(m, np.inf),
                                                 -int(recurrence[m])))
                cap = st.slider("Max miRNAs in heatmap", 10, max(20, len(idx)),
                                min(60, len(idx)), 5, key="matrix_cap")
                idx = idx[:cap]

                def _ylab(m):
                    p = padj_map.get(m)
                    return f"★ {m}" if (p is not None and p <= sig_cut) else m

                sub = mat.loc[idx].rename(columns=D.ds_label)
                sub.index = [_ylab(m) for m in idx]
                fig = px.imshow(sub, color_continuous_scale=["#1f77b4", "#dddddd", "#d62728"],
                                zmin=-1, zmax=1, aspect="auto",
                                labels=dict(color="direction", x="dataset", y="miRNA"))
                fig.update_layout(height=max(360, 18 * len(sub)), margin=dict(t=10, b=10),
                                  coloraxis_colorbar=dict(tickvals=[-1, 0, 1],
                                                          ticktext=["DOWN", "mixed", "UP"]))
                st.plotly_chart(fig, use_container_width=True)
                st.caption("**How to read:** rows = miRNAs in ≥N datasets (★ = passes the "
                           "combined-padj cutoff), columns = datasets. Red = up, blue = "
                           "down, grey = mixed across that dataset's comparisons, white = "
                           "not DE there. A row that's all red or all blue **and** starred "
                           "is a reproducible, same-direction, statistically-backed signal "
                           "— your strongest cross-tissue biomarker candidates.")

                # Companion table: per-dataset direction + combined stats, downloadable.
                dirs = (mat.loc[idx].rename(columns=D.ds_label)
                        .replace({1: "UP", -1: "DOWN", 0: "mixed"}))
                dirs = dirs.where(dirs.notna(), "—")
                out_tbl = pd.DataFrame({
                    "miRNA": idx,
                    "n_datasets": [int(recurrence[m]) for m in idx],
                    "combined_padj": [padj_map.get(m) for m in idx],
                    "consensus_dir": [cdir_map.get(m, "—") for m in idx],
                })
                out_tbl = pd.concat([out_tbl, dirs.reset_index(drop=True)], axis=1)
                st.markdown("**Direction & combined significance**")
                D.df_with_copy(out_tbl, use_container_width=True, hide_index=True,
                             column_config=D.sci_format("combined_padj"))
                st.download_button("⬇️ Download direction matrix (CSV)",
                                   out_tbl.to_csv(index=False).encode(),
                                   file_name="direction_matrix.csv", mime="text/csv",
                                   key="dl_direction_matrix")

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
        st.caption("**How to read:** each point is a miRNA. X = combined Stouffer Z — "
                   "right means consensus **UP** across studies, left means consensus "
                   "**DOWN**. Y = combined significance (higher = stronger; dashed line "
                   "is your padj cutoff). Dot size = how many datasets back it. The "
                   "top-left and top-right points are your reproducible signature.")

        show = meta[["miRNA", "n_datasets", "direction", "combined_Z",
                     "mean_log2FC", "combined_p", "combined_padj", "datasets"]]
        D.df_with_copy(show, use_container_width=True, hide_index=True,
                     column_config=D.sci_format("combined_p", "combined_padj"))
        st.download_button("⬇️ Download consensus signature (CSV)",
                           show.to_csv(index=False).encode(),
                           file_name="meta_analysis_consensus.csv", mime="text/csv",
                           key="dl_meta_consensus")

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
            fig = px.imshow(corr.rename(index=D.ds_label, columns=D.ds_label),
                            text_auto=".2f", color_continuous_scale="RdBu_r",
                            zmin=-1, zmax=1, aspect="auto")
            fig.update_layout(height=320, margin=dict(t=10, b=10))
            st.plotly_chart(fig, use_container_width=True)
            st.caption("**How to read:** Spearman ρ of log2FC between each pair of "
                       "datasets on their shared miRNAs. +1 (red) = they agree on "
                       "direction & magnitude; 0 (white) = unrelated; −1 (blue) = "
                       "opposite. Cross-tissue pairs (blood vs mucosa) usually agree least.")

        pcols = st.columns(2)
        dx = pcols[0].selectbox("X axis", list(wide.columns), index=0, key="conc_x",
                                format_func=D.ds_label)
        dy = pcols[1].selectbox("Y axis", list(wide.columns),
                                index=1 if wide.shape[1] > 1 else 0, key="conc_y",
                                format_func=D.ds_label)
        if dx != dy:
            pair = wide[[dx, dy]].dropna()
            pair = pair.reset_index().rename(columns={"index": "miRNA"})
            if len(pair) > 2:
                rho = spearmanr(pair[dx], pair[dy]).statistic
                lim = float(np.nanmax(np.abs(pair[[dx, dy]].values))) * 1.05
                fig = px.scatter(pair, x=dx, y=dy, hover_name="miRNA", opacity=0.7,
                                 labels={dx: f"log2FC — {D.ds_label(dx)}",
                                         dy: f"log2FC — {D.ds_label(dy)}"})
                fig.add_shape(type="line", x0=-lim, y0=-lim, x1=lim, y1=lim,
                              line=dict(dash="dot", color="grey"))
                fig.add_hline(y=0, line_color="#eee"); fig.add_vline(x=0, line_color="#eee")
                fig.update_layout(height=460, margin=dict(t=30, b=10),
                                  title=f"{len(pair)} shared miRNAs · Spearman ρ = {rho:.2f}")
                st.plotly_chart(fig, use_container_width=True)
                st.caption("**How to read:** each point is a miRNA shared by the two "
                           "datasets; axes are its log2FC in each. Points hugging the "
                           "dotted diagonal = same direction & size in both "
                           "(reproducible); a shapeless cloud = no agreement; the "
                           "opposite diagonal = opposite responses. ρ in the title "
                           "summarizes the overall agreement.")
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
                D.df_with_copy(pc, use_container_width=True, hide_index=True)
                st.download_button(f"⬇️ Download shared {kind} terms (CSV)",
                                   pc.to_csv(index=False).encode(),
                                   file_name=f"shared_{kind}_terms.csv", mime="text/csv",
                                   key=f"dl_shared_{kind}")

# ══ IBD-TARGET OVERLAP (validated-target enrichment, per dataset) ══════════════
with tab_ibd:
    st.markdown("#### IBD-target-overlap miRNAs across datasets")
    st.caption("Which miRNAs have **validated targets enriched for known IBD genes**, "
               "and do they recur across datasets? Built from each dataset's "
               "`ibd_target_overlap.csv` (Fisher's exact test) — independent of the "
               "DE thresholds at the top of the page.")

    ic = st.columns(3)
    metric = ic[0].selectbox(
        "Significance metric", ["padj", "p_fisher"],
        format_func=lambda m: {"padj": "BH-adjusted p (padj)",
                               "p_fisher": "raw Fisher p"}[m],
        key="ibd_metric")
    cutoff = ic[1].slider(f"{metric} ≤", 0.0, 1.0, 0.05, 0.005, key="ibd_cut")
    enriched_only = ic[2].checkbox("Enriched only (OR > 1)", value=True, key="ibd_or")

    ibd_long = D.ibd_overlap_long(metric=metric, cutoff=cutoff,
                                  enriched_only=enriched_only)
    if ibd_long.empty:
        st.info("No IBD-overlap miRNAs pass these filters — loosen the cutoff or "
                "untick 'Enriched only'.")
    else:
        sets = {ds: set(g["miRNA"]) for ds, g in ibd_long.groupby("dataset")}
        all_mir = sorted(set().union(*sets.values()))
        rec = {m: sum(m in s for s in sets.values()) for m in all_mir}

        m1, m2, m3 = st.columns(3)
        m1.metric("Datasets with IBD data", len(sets))
        m2.metric("Unique IBD miRNAs", len(all_mir))
        m3.metric("Recurrent (≥2)", sum(v >= 2 for v in rec.values()))

        # Tell the user which datasets have IBD data but no miRNA passing the filter,
        # so a missing dataset in the Venn/heatmap isn't mistaken for a bug.
        with_file = [d["id"] for d in analyzed
                     if D.ibd_overlap_table(D.resolve(d["output_dir"])) is not None]
        dropped = [x for x in with_file if x not in sets]
        if dropped:
            st.caption(
                f"⚠️ No miRNA passes the current filter for: **{', '.join(dropped)}** "
                f"(IBD data present, but nothing meets {metric} ≤ {cutoff}"
                + (" with OR > 1" if enriched_only else "")
                + "). They're excluded from the Venn and heatmap — loosen the cutoff "
                "or switch to raw Fisher p to include them.")

        # Venn of the IBD-overlap miRNA sets (pick 2–3 datasets).
        if len(sets) >= 2:
            st.markdown("**Venn of IBD miRNA sets**")
            ibd_chosen = st.multiselect("Datasets in the Venn (2–3)", list(sets.keys()),
                                        default=list(sets.keys())[:3],
                                        format_func=D.ds_label, key="ibd_venn_sets")
            if len(ibd_chosen) in (2, 3):
                fig_v, ax_v = plt.subplots(figsize=(5, 4.5))
                vlabels = [D.ds_label(s) for s in ibd_chosen]
                vdata = [sets[c] for c in ibd_chosen]
                vobj = (venn2(vdata, set_labels=vlabels, ax=ax_v) if len(ibd_chosen) == 2
                        else venn3(vdata, set_labels=vlabels, ax=ax_v))
                for t in (vobj.set_labels or []):
                    if t:
                        t.set_fontsize(8)
                D.zoomable_pyplot(fig_v)
                plt.close(fig_v)
                st.caption("**How to read:** like the Overlap Venn, but on the "
                           "IBD-target-enriched miRNAs — circles are datasets, overlaps "
                           "are miRNAs shared at the current metric/cutoff/OR filter.")

                named_v = {D.ds_label(c): sets[c] for c in ibd_chosen}
                regions_v = D.overlap_regions(named_v)
                regv_df = pd.DataFrame(
                    [{"Region": r, "Count": len(mem), "miRNAs": ", ".join(mem)}
                     for r, mem in regions_v.items()]
                )
                D.df_with_copy(regv_df, use_container_width=True, hide_index=True)
                tidy_v = pd.DataFrame(
                    [(r, m) for r, mem in regions_v.items() for m in mem],
                    columns=["Region", "miRNA"])
                st.download_button("⬇️ Download IBD Venn regions (CSV)",
                                   tidy_v.to_csv(index=False).encode(),
                                   file_name="ibd_venn_regions.csv", mime="text/csv",
                                   key="dl_ibd_venn")
            else:
                st.info("Select 2 or 3 datasets for the Venn.")
            st.divider()

        # Odds-ratio matrix (miRNA × dataset), miRNAs ordered by recurrence then OR.
        # This IBD layer is a per-dataset Fisher test (no Stouffer), so the
        # cross-dataset "consensus" analog is recurrence: a ★ marks miRNAs that
        # pass the current metric/cutoff/OR filter in ≥ N datasets.
        max_rec = max(rec.values()) if rec else 1
        hc = st.columns(2)
        min_ds_ibd = hc[0].slider("Require IBD miRNA in ≥ N datasets", 1,
                                  max(1, max_rec), 1, key="ibd_minds")
        only_rec = hc[1].checkbox(f"Show only miRNAs in ≥{min_ds_ibd} datasets",
                                  value=min_ds_ibd > 1, key="ibd_onlyrec")

        ormat = ibd_long.pivot_table(index="miRNA", columns="dataset", values="OR",
                                     aggfunc="first")
        order = sorted(ormat.index, key=lambda m: (-rec[m], -ormat.loc[m].mean()))
        if only_rec:
            order = [m for m in order if rec[m] >= min_ds_ibd]

        if not order:
            st.info(f"No IBD miRNA recurs in ≥{min_ds_ibd} datasets at the current "
                    "filter. Lower the requirement or loosen the metric/cutoff.")
        else:
            if len(order) > 12:
                cap = st.slider("Max miRNAs in heatmap", 5, len(order),
                                min(40, len(order)), 1, key="ibd_cap")
            else:
                cap = len(order)
            shown = order[:cap]
            ormat = ormat.loc[shown]
            # ★ marks miRNAs recurrent in ≥ N datasets (the consensus highlight).
            ylabels = [f"★ {m}" if rec[m] >= min_ds_ibd else m for m in shown]
            ormat.index = ylabels
            ormat = ormat.rename(columns=D.ds_label)
            fig = px.imshow(ormat, color_continuous_scale="RdBu_r",
                            color_continuous_midpoint=1, aspect="auto",
                            labels=dict(color="odds ratio", x="dataset", y="miRNA"))
            fig.update_layout(height=max(360, 20 * len(ormat)), margin=dict(t=10, b=10))
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"Red = enriched (OR>1), blue = depleted (OR<1), white ≈ 1, "
                       f"blank = miRNA absent from that dataset's IBD table. "
                       f"★ = IBD-enriched in ≥{min_ds_ibd} datasets (recurrent "
                       "consensus candidate).")

        # Recurrence summary.
        summary = (ibd_long.groupby("miRNA")
                   .agg(n_datasets=("dataset", "nunique"),
                        datasets=("dataset", lambda s: ", ".join(sorted(set(s)))),
                        mean_OR=("OR", "mean"),
                        best=(metric, "min"))
                   .reset_index()
                   .sort_values(["n_datasets", "mean_OR"], ascending=[False, False])
                   .rename(columns={"best": f"min_{metric}"})
                   .reset_index(drop=True))
        st.markdown("**Recurrence summary**")
        D.df_with_copy(summary, use_container_width=True, hide_index=True,
                     column_config=D.sci_format(f"min_{metric}"))
        st.download_button("⬇️ Download IBD cross-dataset table (CSV)",
                           ibd_long.to_csv(index=False).encode(),
                           file_name="ibd_overlap_cross_dataset.csv", mime="text/csv",
                           key="dl_ibd_cross")

# ══ VALIDATED vs PREDICTED vs BASELINE (cross-dataset) ════════════════════════
with tab_cmp:
    st.markdown("#### IBD-target enrichment by evidence layer, across datasets")
    st.caption(
        "Three layers of the *same* Fisher IBD-overlap test, pooled across datasets. "
        "**Validated** / **Predicted** = the DE-significant miRNAs tested on their "
        "validated vs computationally-predicted targets; **Baseline** = *every* "
        "DEA-tested miRNA (no cutoff) — the reference for an arbitrary miRNA. Built "
        "from each dataset's `ibd_target_overlap_comparison.csv` "
        "(15_multimir_targets_baseline.R → 19_ibd_target_overlap.R).")

    cmp_long = D.ibd_comparison_long()
    if cmp_long.empty:
        st.info("No comparison files found. Run 06_targets_enrich/15_multimir_targets_"
                "baseline.R then 07_networks_clin/19_ibd_target_overlap.R in each dataset.")
    elif not {"evidence", "is_significant", "OR", "pct_ibd", "padj"}.issubset(cmp_long.columns):
        st.warning("Comparison files are missing expected columns — re-run script 19.")
    else:
        present = sorted(cmp_long["dataset"].unique())
        st.caption(f"Datasets with comparison data: **{len(present)}** "
                   f"({', '.join(D.ds_label(x) for x in present)}).")

        cc = st.columns(2)
        metric = cc[0].radio("Metric", ["OR", "pct_ibd"], horizontal=True,
                             format_func=lambda m: {"OR": "Odds ratio",
                                                    "pct_ibd": "% IBD targets"}[m],
                             key="cmp_metric")
        sig_cut = cc[1].slider("Significant: padj ≤", 0.0, 0.25, 0.05, 0.005,
                               key="cmp_sig")
        mlabel = "Median odds ratio" if metric == "OR" else "Median % IBD targets"

        # 1) Four paired bars per dataset: each evidence's significant miRNAs next to
        # ITS OWN all-miRNA baseline (same universe), so the comparison stays honest —
        # the predicted universe sits higher on its own, which a single shared baseline
        # bar would hide. Red family = validated, blue family = predicted; dark =
        # significant, light = baseline.
        L_VS, L_VB = "Validated · significant", "Validated · baseline"
        L_PS, L_PB = "Predicted · significant", "Predicted · baseline"
        rows = []
        for ds, g in cmp_long.groupby("dataset"):
            v = g[g["evidence"] == "validated"]; p = g[g["evidence"] == "predicted"]
            lab = D.ds_label(ds)
            rows += [
                {"dataset": lab, "layer": L_VS, metric: v.loc[v["is_significant"], metric].median()},
                {"dataset": lab, "layer": L_VB, metric: v[metric].median()},
                {"dataset": lab, "layer": L_PS, metric: p.loc[p["is_significant"], metric].median()},
                {"dataset": lab, "layer": L_PB, metric: p[metric].median()},
            ]
        bar = pd.DataFrame(rows)
        cmap = {L_VS: "#d62728", L_VB: "#f2b6b6", L_PS: "#1f77b4", L_PB: "#a9cce3"}
        fig = px.bar(bar, x="dataset", y=metric, color="layer", barmode="group",
                     color_discrete_map=cmap,
                     category_orders={"layer": [L_VS, L_VB, L_PS, L_PB]},
                     labels={metric: mlabel, "layer": ""})
        if metric == "OR":
            fig.add_hline(y=1.0, line_color="#ccc")
        fig.update_layout(height=400, margin=dict(t=30, b=10))
        st.plotly_chart(fig, use_container_width=True)
        st.caption("**How to read:** four bars per dataset — each evidence's "
                   "**significant** miRNAs (dark) beside *its own* **baseline** of all "
                   "miRNAs (light). Compare dark-vs-light *within the same color* "
                   "(same universe): that gap is how much your DE miRNAs beat a typical "
                   "miRNA. The predicted (blue) pair sitting higher than the validated "
                   "(red) pair is a **universe-level** shift, not a stronger DE signal — "
                   "which is why we keep each baseline next to its own significant bar.")

        # 2) Recurrence heatmap: significant miRNAs (in ≥N datasets) × dataset.
        st.divider()
        st.markdown("##### Recurrent significant miRNAs across datasets")
        h1, _ = st.columns([1, 2])
        evidence = h1.radio("Evidence layer (heatmap)", ["validated", "predicted"],
                            horizontal=True, key="cmp_ev")
        sub = cmp_long[cmp_long["evidence"] == evidence].copy()
        sigsub = sub[sub["is_significant"]]
        rec = sigsub.groupby("miRNA")["dataset"].nunique()
        if rec.empty:
            st.info("No significant miRNAs in the comparison for this evidence layer.")
        else:
            max_rec = int(rec.max())
            min_ds = st.slider("Require significant miRNA in ≥ N datasets", 1,
                               max(1, max_rec), 2 if max_rec >= 2 else 1, key="cmp_minds")
            keep = set(rec[rec >= min_ds].index)
            if not keep:
                st.info(f"No significant miRNA recurs in ≥{min_ds} datasets for "
                        f"{evidence} targets. Lower the requirement.")
            else:
                # Show the metric across ALL datasets where each kept miRNA exists.
                hm = sub[sub["miRNA"].isin(keep)]
                mat = hm.pivot_table(index="miRNA", columns="dataset",
                                     values=metric, aggfunc="first")
                order = sorted(mat.index, key=lambda m: (-int(rec.get(m, 0)),
                                                         -mat.loc[m].mean(skipna=True)))
                padj_hits = set(hm.loc[hm["padj"] <= sig_cut, "miRNA"])
                ylabels = [f"★ {m}" if m in padj_hits else m for m in order]
                mat = mat.loc[order].rename(columns=D.ds_label)
                mat.index = ylabels
                if metric == "OR":
                    fig2 = px.imshow(mat, color_continuous_scale="RdBu_r",
                                     color_continuous_midpoint=1, aspect="auto",
                                     labels=dict(color="odds ratio", x="dataset", y="miRNA"))
                else:
                    fig2 = px.imshow(mat, color_continuous_scale="Reds", aspect="auto",
                                     labels=dict(color="% IBD targets", x="dataset", y="miRNA"))
                fig2.update_layout(height=max(300, 26 * len(mat)), margin=dict(t=10, b=10))
                st.plotly_chart(fig2, use_container_width=True)
                st.caption(f"Significant miRNAs recurrent in ≥{min_ds} datasets, shown "
                           f"across every dataset where detected ({evidence} targets). "
                           f"★ = padj ≤ {sig_cut:g} in at least one dataset. Blank = miRNA "
                           "absent from that dataset's comparison. These — especially the "
                           "starred, same-color rows — are your cross-tissue candidates.")

        st.download_button("⬇️ Download pooled comparison (CSV)",
                           cmp_long.to_csv(index=False).encode(),
                           file_name="ibd_comparison_cross_dataset.csv", mime="text/csv",
                           key="dl_cmp_cross")
