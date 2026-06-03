"""
Data access layer for the IBD miRNA-seq dashboard.

Everything is driven by the repo-root `datasets.yaml` manifest and the on-disk
result layout produced by the pipeline. Nothing here is hard-coded per dataset,
so adding a new dataset to the manifest makes it appear in the app automatically.

Result layout assumed per analyzed dataset (output_dir):
    DEA_results/
        <DEA_*_folder>/           DEA_<mode>_<comparison>.csv + volcano/dispersion PNGs
        figures_and_QC/           PCA_plot.png, UMAP_plot.png, heatmap_ALL_miRNAs.png,
                                  vst_normalized_matrix.csv
    pathway_outputs/              GO/KEGG_Enrichment_Results.csv + *_Dotplot.png
    multimir_outputs/             validated/predicted target CSVs
    ibd_overlap_outputs/          ibd_target_overlap.csv, ibd_gene_list.csv, ...
    hmdd_outputs/                 hmdd_mirna_table.csv
    network_outputs/              Automated_Dashboard_{GO,KEGG}.html
    Final_Report.html
"""
from __future__ import annotations

import re
from pathlib import Path
from functools import lru_cache

import numpy as np
import pandas as pd
import yaml
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

# Repo root = .../pipelines/08_report/dashboard/data.py -> parents[3]
REPO_ROOT = Path(__file__).resolve().parents[3]
MANIFEST_PATH = REPO_ROOT / "datasets.yaml"

STATUS_ORDER = {"analyzed": 0, "in_progress": 1, "suspended": 2}
STATUS_LABEL = {
    "analyzed": "✅ Analyzed",
    "in_progress": "🟡 In progress",
    "suspended": "⏸️ Suspended",
}


# ── Path helpers ──────────────────────────────────────────────────────────────

def resolve(path_str: str | None) -> Path | None:
    """Resolve a manifest path (relative to repo root) to an absolute Path."""
    if not path_str:
        return None
    p = Path(path_str)
    if not p.is_absolute():
        p = REPO_ROOT / p
    return p


# ── Manifest ──────────────────────────────────────────────────────────────────

@lru_cache(maxsize=1)
def load_manifest() -> dict:
    with open(MANIFEST_PATH, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def datasets(status: str | None = None) -> list[dict]:
    ds = load_manifest().get("datasets", [])
    if status:
        ds = [d for d in ds if d.get("status") == status]
    return sorted(ds, key=lambda d: (STATUS_ORDER.get(d.get("status"), 9), d["id"]))


def analyzed() -> list[dict]:
    return datasets(status="analyzed")


def get_dataset(ds_id: str) -> dict | None:
    return next((d for d in load_manifest()["datasets"] if d["id"] == ds_id), None)


def manifest_overview_df() -> pd.DataFrame:
    rows = []
    for d in datasets():
        rows.append({
            "Dataset": d["id"],
            "Status": STATUS_LABEL.get(d.get("status"), d.get("status")),
            "Tissue": d.get("tissue") or "—",
            "Comparisons": ", ".join(d.get("comparisons") or []) or "—",
            "Samples": d.get("n_samples_analyzed") or "—",
            "Runs": d.get("n_runs_total") or "—",
            "Strategy": d.get("library_strategy") or "—",
            "Report": "Yes" if resolve(d.get("final_report")) and resolve(d["final_report"]).exists() else "—",
        })
    return pd.DataFrame(rows)


# ── DEA discovery & loading ───────────────────────────────────────────────────

_MODE_LABELS = {
    "Default": "Default (outlier replacement on)",
    "Strict": "Strict (no replacement)",
}
_DEA_RE = re.compile(r"^DEA_([A-Za-z]+)_(.+)\.csv$")


def discover_dea(output_dir: Path) -> dict[str, dict[str, Path]]:
    """Return {mode: {comparison: csv_path}} discovered under DEA_results/."""
    root = output_dir / "DEA_results"
    out: dict[str, dict[str, Path]] = {}
    if not root.exists():
        return out
    for csv in root.rglob("DEA_*.csv"):
        m = _DEA_RE.match(csv.name)
        if not m:
            continue
        mode, comparison = m.group(1), m.group(2)
        out.setdefault(mode, {})[comparison] = csv
    return out


def load_dea(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # Normalize the miRNA id column name
    first = df.columns[0]
    if first.lower() in {"mirna", "mir", "unnamed: 0", ""}:
        df = df.rename(columns={first: "miRNA"})
    for col in ("log2FoldChange", "padj", "pvalue", "baseMean"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def dea_significant(df: pd.DataFrame, padj: float, lfc: float) -> pd.DataFrame:
    sig = df[(df["padj"].notna()) & (df["padj"] <= padj) & (df["log2FoldChange"].abs() >= lfc)]
    return sig.copy()


def add_direction(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["direction"] = np.where(df["log2FoldChange"] >= 0, "UP", "DOWN")
    return df


# ── Figures & tabular result files ────────────────────────────────────────────

def qc_figures(output_dir: Path) -> dict[str, Path]:
    base = output_dir / "DEA_results" / "figures_and_QC"
    names = {
        "PCA": "PCA_plot.png",
        "UMAP": "UMAP_plot.png",
        "Heatmap (all miRNAs)": "heatmap_ALL_miRNAs.png",
    }
    return {label: base / fn for label, fn in names.items() if (base / fn).exists()}


def pathway_table(output_dir: Path, kind: str) -> pd.DataFrame | None:
    """kind in {'GO', 'KEGG'}."""
    f = output_dir / "pathway_outputs" / f"{kind}_Enrichment_Results.csv"
    if not f.exists():
        return None
    return pd.read_csv(f)


def pathway_dotplot(output_dir: Path, kind: str) -> Path | None:
    f = output_dir / "pathway_outputs" / f"{kind}_Dotplot.png"
    return f if f.exists() else None


def ibd_overlap_table(output_dir: Path) -> pd.DataFrame | None:
    f = output_dir / "ibd_overlap_outputs" / "ibd_target_overlap.csv"
    if not f.exists():
        return None
    return pd.read_csv(f)


def sample_metadata(d: dict) -> pd.DataFrame | None:
    f = resolve(d.get("metadata_csv"))
    if not f or not f.exists():
        return None
    return pd.read_csv(f)


def network_html(output_dir: Path) -> dict[str, Path]:
    base = output_dir / "network_outputs"
    out = {}
    for kind in ("GO", "KEGG"):
        f = base / f"Automated_Dashboard_{kind}.html"
        if f.exists():
            out[kind] = f
    return out


# ── Cross-dataset meta-analysis ───────────────────────────────────────────────

def collect_de_mirnas(padj: float, lfc: float, mode: str = "Default") -> pd.DataFrame:
    """
    Long-format DE-miRNA table across all analyzed datasets.

    A miRNA is collected if it is significant in *any* comparison of the dataset
    (for the chosen replacement `mode`). Columns:
        miRNA, dataset, comparison, log2FoldChange, padj, direction
    """
    rows = []
    for d in analyzed():
        odir = resolve(d["output_dir"])
        dea = discover_dea(odir)
        comp_map = dea.get(mode) or (next(iter(dea.values())) if dea else {})
        for comparison, csv in comp_map.items():
            try:
                df = load_dea(csv)
            except Exception:
                continue
            sig = dea_significant(df, padj, lfc)
            if sig.empty:
                continue
            sig = add_direction(sig)
            for _, r in sig.iterrows():
                rows.append({
                    "miRNA": r["miRNA"],
                    "dataset": d["id"],
                    "comparison": comparison,
                    "log2FoldChange": r["log2FoldChange"],
                    "padj": r["padj"],
                    "direction": r["direction"],
                })
    return pd.DataFrame(rows)


def consensus_matrix(de_long: pd.DataFrame):
    """
    Build a miRNA × dataset matrix of signed presence.

    Cell value: +1 (UP in that dataset), -1 (DOWN), 0 (mixed across comparisons),
    NaN (not DE). Returns (matrix_df, recurrence_series).
    """
    if de_long.empty:
        return pd.DataFrame(), pd.Series(dtype=int)

    # Per (miRNA, dataset): consensus direction across comparisons
    def _dir(group):
        dirs = set(group["direction"])
        if dirs == {"UP"}:
            return 1
        if dirs == {"DOWN"}:
            return -1
        return 0  # mixed

    signed = (de_long.groupby(["miRNA", "dataset"])
              .apply(_dir, include_groups=False)
              .rename("val").reset_index())
    mat = signed.pivot(index="miRNA", columns="dataset", values="val")
    recurrence = mat.notna().sum(axis=1).sort_values(ascending=False)
    mat = mat.loc[recurrence.index]
    return mat, recurrence


def de_sets(padj: float, lfc: float, mode: str = "Default") -> dict[str, set]:
    """{dataset_id: set(DE miRNAs)} — a miRNA is DE if significant in ANY of the
    dataset's comparisons. Drives the Venn and UpSet views."""
    de_long = collect_de_mirnas(padj, lfc, mode=mode)
    if de_long.empty:
        return {}
    return {ds: set(g["miRNA"]) for ds, g in de_long.groupby("dataset")}


def list_comparisons(d: dict, mode: str = "Default") -> list[str]:
    dea = discover_dea(resolve(d["output_dir"]))
    comp_map = dea.get(mode) or next(iter(dea.values()), {})
    return list(comp_map.keys())


def load_comparison(d: dict, comparison: str, mode: str = "Default") -> pd.DataFrame | None:
    dea = discover_dea(resolve(d["output_dir"]))
    comp_map = dea.get(mode) or next(iter(dea.values()), {})
    csv = comp_map.get(comparison)
    return load_dea(csv) if csv is not None else None


def meta_stouffer(selections: dict[str, str], mode: str = "Default",
                  min_datasets: int = 2) -> pd.DataFrame:
    """
    Sample-size-weighted Stouffer meta-analysis across datasets.

    `selections` maps dataset_id -> the single comparison to use for that dataset.
    For each miRNA, per-dataset two-sided p-values are converted to signed z-scores
    (sign from log2FC), combined as Z = Σ wᵢzᵢ / √Σwᵢ², with wᵢ = √n_samples.
    Returns one ranked row per miRNA present in ≥ min_datasets, with a BH-adjusted
    combined p-value. This is the formal "consensus signature".
    """
    rows = []
    for ds_id, comparison in selections.items():
        if not comparison:
            continue
        d = get_dataset(ds_id)
        df = load_comparison(d, comparison, mode)
        if df is None:
            continue
        sub = df[["miRNA", "log2FoldChange", "pvalue"]].dropna(subset=["pvalue"])
        if sub.empty:
            continue
        w = float(d.get("n_samples_analyzed") or 1) ** 0.5
        p = sub["pvalue"].clip(1e-300, 1 - 1e-16)
        z = np.sign(sub["log2FoldChange"]) * norm.isf(p / 2)
        rows.append(pd.DataFrame({"miRNA": sub["miRNA"].values, "dataset": ds_id,
                                  "z": z.values, "w": w,
                                  "lfc": sub["log2FoldChange"].values}))
    if not rows:
        return pd.DataFrame()
    long = pd.concat(rows, ignore_index=True)

    def _combine(x):
        return pd.Series({
            "n_datasets": x["dataset"].nunique(),
            "combined_Z": (x["w"] * x["z"]).sum() / np.sqrt((x["w"] ** 2).sum()),
            "mean_log2FC": x["lfc"].mean(),
            "datasets": ", ".join(sorted(x["dataset"].unique())),
        })

    agg = long.groupby("miRNA").apply(_combine, include_groups=False).reset_index()
    agg = agg[agg["n_datasets"] >= min_datasets].copy()
    if agg.empty:
        return agg
    agg["combined_p"] = 2 * norm.sf(agg["combined_Z"].abs())
    agg["combined_padj"] = multipletests(agg["combined_p"], method="fdr_bh")[1]
    agg["direction"] = np.where(agg["combined_Z"] >= 0, "UP", "DOWN")
    agg["n_datasets"] = agg["n_datasets"].astype(int)
    return agg.sort_values("combined_p").reset_index(drop=True)


def effect_size_wide(selections: dict[str, str], mode: str = "Default") -> pd.DataFrame:
    """Wide log2FC matrix (index miRNA, columns dataset) for concordance plots."""
    series = {}
    for ds_id, comparison in selections.items():
        if not comparison:
            continue
        df = load_comparison(get_dataset(ds_id), comparison, mode)
        if df is None:
            continue
        series[ds_id] = df.dropna(subset=["log2FoldChange"]).set_index("miRNA")["log2FoldChange"]
    return pd.DataFrame(series) if series else pd.DataFrame()


def pathway_consensus(kind: str, padj: float = 0.05, top: int | None = None) -> pd.DataFrame:
    """
    Count how many analyzed datasets share each enriched term (p.adjust <= padj).
    kind in {'GO','KEGG'}.
    """
    counts: dict[str, dict] = {}
    for d in analyzed():
        tbl = pathway_table(resolve(d["output_dir"]), kind)
        if tbl is None or "Description" not in tbl.columns:
            continue
        padj_col = "p.adjust" if "p.adjust" in tbl.columns else "pvalue"
        sig = tbl[pd.to_numeric(tbl[padj_col], errors="coerce") <= padj]
        for desc in sig["Description"].dropna().unique():
            entry = counts.setdefault(desc, {"Term": desc, "n_datasets": 0, "datasets": []})
            entry["n_datasets"] += 1
            entry["datasets"].append(d["id"])
    if not counts:
        return pd.DataFrame()
    out = pd.DataFrame(counts.values())
    out["datasets"] = out["datasets"].apply(lambda xs: ", ".join(sorted(xs)))
    out = out.sort_values(["n_datasets", "Term"], ascending=[False, True]).reset_index(drop=True)
    return out.head(top) if top else out
