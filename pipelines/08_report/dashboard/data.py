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

import base64
import os
import re
import sys
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
SYNC_SCRIPT = REPO_ROOT / "pipelines" / "00_search" / "05_sync_yaml.py"
DOCS_DIR = REPO_ROOT / "outputs" / "docs"

# Leading YAML frontmatter block delimited by `---` fences. Tags are read ONLY
# from here; inline #hashtags elsewhere in the doc are ignored.
_FRONTMATTER_RE = re.compile(r"^\s*---\s*\n(.*?)\n---\s*(?:\n|$)", re.DOTALL)

# Full-screen blocking overlay shown while the app is loading. It is hidden by
# default and revealed purely by CSS whenever Streamlit's running indicator
# (`stStatusWidget`) is in the DOM, so it needs no per-action wiring. Render it
# once near the top of every page: st.markdown(D.LOADING_OVERLAY, unsafe_allow_html=True)
LOADING_OVERLAY = """
<style>
#stLoadingOverlay {
    position: fixed; inset: 0; z-index: 2147483646;
    display: none; align-items: center; justify-content: center;
    background: rgba(14, 17, 23, 0.45);
    -webkit-backdrop-filter: blur(1.5px); backdrop-filter: blur(1.5px);
    cursor: wait;
}
/* Reveal only while a script run is in progress (running indicator present). */
[data-testid="stApp"]:has([data-testid="stStatusWidget"]) #stLoadingOverlay {
    display: flex;
}
#stLoadingOverlay .lo-dots { display: flex; gap: 14px; }
#stLoadingOverlay .lo-dots span {
    width: 16px; height: 16px; border-radius: 50%;
    background: #1f77b4;
    animation: lo-bounce 1.2s ease-in-out infinite;
}
#stLoadingOverlay .lo-dots span:nth-child(2) { animation-delay: 0.16s; }
#stLoadingOverlay .lo-dots span:nth-child(3) { animation-delay: 0.32s; }
@keyframes lo-bounce {
    0%, 80%, 100% { transform: translateY(0);     opacity: 0.55; }
    40%           { transform: translateY(-18px);  opacity: 1; }
}
</style>
<div id="stLoadingOverlay"><div class="lo-dots"><span></span><span></span><span></span></div></div>
"""

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

def find_doc(d: dict) -> Path | None:
    """Locate the outputs/docs/*.md for a dataset (id → accession → bioproject)."""
    for key in ("id", "accession", "bioproject"):
        stem = d.get(key)
        if stem:
            p = DOCS_DIR / f"{stem}.md"
            if p.exists():
                return p
    return None


def doc_tags(d: dict) -> list[str]:
    """Tags from the YAML frontmatter `tags:` list at the top of the dataset's doc.

    Only a leading `---`-fenced frontmatter block is read, e.g.::

        ---
        tags:
          - finalized
          - network_pending
        ---

    A scalar (`tags: finalized` or comma/space-separated) is also accepted; any
    leading '#' is stripped. Returns a deduped, order-preserving list.
    """
    doc = find_doc(d)
    if not doc:
        return []
    try:
        text = doc.read_text(encoding="utf-8")
    except OSError:
        return []
    m = _FRONTMATTER_RE.match(text)
    if not m:
        return []
    try:
        meta = yaml.safe_load(m.group(1)) or {}
    except yaml.YAMLError:
        return []
    raw = meta.get("tags") if isinstance(meta, dict) else None
    if raw is None:
        return []
    if isinstance(raw, str):
        raw = raw.replace(",", " ").split()
    elif not isinstance(raw, (list, tuple)):
        raw = [raw]
    seen, tags = set(), []
    for t in raw:
        t = str(t).lstrip("#").strip()
        if t and t.lower() not in seen:
            seen.add(t.lower())
            tags.append(t)
    return tags


def all_tags() -> list[str]:
    """Sorted union of tags across every dataset (drives the sidebar filter)."""
    uniq: set[str] = set()
    for d in load_manifest().get("datasets", []):
        uniq.update(d.get("tags", []))
    return sorted(uniq, key=str.lower)


def derive_status(d: dict) -> str:
    """Compute a dataset's pipeline status from on-disk artifacts.

    The filesystem is the source of truth — a dataset is ``analyzed`` once its
    Final_Report.html exists, otherwise ``in_progress``. ``suspended`` is a human
    "parked intentionally" decision that can't be read off disk, so an explicit
    suspended in the manifest is preserved as an override.
    """
    if d.get("status") == "suspended":
        return "suspended"
    report = resolve(d.get("final_report"))
    if report is None:
        odir = resolve(d.get("output_dir"))
        report = (odir / "Final_Report.html") if odir else None
    return "analyzed" if (report and report.exists()) else "in_progress"


@lru_cache(maxsize=1)
def load_manifest() -> dict:
    with open(MANIFEST_PATH, encoding="utf-8") as fh:
        m = yaml.safe_load(fh)
    # Status is derived from disk, not trusted from the file. Keep what the file
    # said as `status_yaml` so the UI can flag drift between the two.
    for d in m.get("datasets", []):
        d["status_yaml"] = d.get("status")
        d["status"] = derive_status(d)       # pipeline progress, from disk
        d["tags"] = doc_tags(d)              # human/review layer, from md docs
    return m


def clear_caches() -> None:
    """Drop every lru_cache that reads from disk so the next render is fresh.

    Streamlit reruns the script on each interaction but keeps this module
    loaded, so the caches would otherwise survive a manifest edit.
    """
    load_manifest.cache_clear()
    card_dea_stats.cache_clear()
    find_qc_reports.cache_clear()
    de_universe.cache_clear()


def _find_conda() -> str | None:
    """Locate the conda binary the way 21_launch_dashboard.sh does."""
    import shutil
    home = Path.home()
    for c in (home / "miniconda3" / "bin" / "conda",
              home / "anaconda3" / "bin" / "conda"):
        if c.exists():
            return str(c)
    return shutil.which("conda")


def run_sync(add_new: bool = False, timeout: int = 120) -> tuple[bool, str]:
    """Run 05_sync_yaml.py --confirm to refill datasets.yaml from outputs/docs.

    The script needs ruamel.yaml (lives in py_env, not the dashboard's env), so
    it is invoked via `conda run -n py_env`. Falls back to the current
    interpreter if conda can't be found. Returns (ok, combined stdout+stderr).
    """
    import subprocess
    # Pass absolute --yaml/--docs so the script's cwd-relative defaults
    # (../../datasets.yaml) can't resolve to the wrong directory.
    args = [
        str(SYNC_SCRIPT),
        "--yaml", str(MANIFEST_PATH),
        "--docs", str(REPO_ROOT / "outputs" / "docs"),
        "--confirm",
    ]
    if add_new:
        args.append("--add-new")

    conda = _find_conda()
    if conda:
        cmd = [conda, "run", "--no-capture-output", "-n", "py_env", "python", *args]
    else:
        cmd = [sys.executable, *args]

    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           cwd=str(REPO_ROOT), timeout=timeout)
        return p.returncode == 0, (p.stdout + p.stderr).strip()
    except Exception as e:  # noqa: BLE001 — surface any launch failure to the UI
        return False, f"Could not run the sync script: {e}"


def datasets(status: str | None = None) -> list[dict]:
    ds = load_manifest().get("datasets", [])
    if status:
        ds = [d for d in ds if d.get("status") == status]
    return sorted(ds, key=lambda d: (STATUS_ORDER.get(d.get("status"), 9), d["id"]))


def analyzed() -> list[dict]:
    return datasets(status="analyzed")


def get_dataset(ds_id: str) -> dict | None:
    return next((d for d in load_manifest()["datasets"] if d["id"] == ds_id), None)


def ds_label(ds_id: str) -> str:
    """Display label tagging a dataset id with its tissue, e.g. 'GSE272890 [colon tissue]'.

    Used as the `format_func` for selectboxes/multiselects and for figure labels so
    it's easy to remember which dataset is which. Falls back to the bare id."""
    d = get_dataset(ds_id)
    tissue = d.get("tissue") if d else None
    return f"{ds_id} [{tissue}]" if tissue else ds_id


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


def sci_format(*cols: str, fmt: str = "%.2e") -> dict:
    """Streamlit ``column_config`` that renders the given columns in scientific
    notation (e.g. ``4.00e-15``) while keeping them numerically sortable — so tiny
    p-values / padj stop showing as ``0.000000000000004``. Streamlit is imported
    lazily so the data layer stays usable outside the dashboard.
    """
    import streamlit as st
    return {c: st.column_config.NumberColumn(format=fmt) for c in cols}


def _df_to_markdown(df: pd.DataFrame) -> str:
    """GitHub/Obsidian-flavoured Markdown table (pipes escaped, newlines flattened)."""
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |",
             "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in df.iterrows():
        cells = ["" if pd.isna(v) else str(v).replace("|", "\\|").replace("\n", " ")
                 for v in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def df_with_copy(df: pd.DataFrame, *, copy_label: str = "📋 Copy table (Markdown / CSV / TSV)",
                 max_copy_rows: int = 2000, **kwargs) -> None:
    """st.dataframe + a collapsed expander that exposes the table as copyable text.

    Streamlit's grid only copies as TSV (great for Word/Excel, not for Obsidian),
    so this adds a **Markdown** block (pastes as a real table in Obsidian) and a
    **TSV** block (pastes as a table in Word/Excel). Both use st.code, whose
    built-in top-right icon copies the whole block in one click.

    Extra **kwargs (use_container_width, hide_index, column_config, …) pass straight
    through to st.dataframe, so call sites swap st.dataframe → D.df_with_copy as-is.
    """
    import streamlit as st
    st.dataframe(df, **kwargs)
    with st.expander(copy_label):
        export, truncated = df, False
        if len(df) > max_copy_rows:
            export, truncated = df.head(max_copy_rows), True
        st.caption("Each box has a copy icon at its top-right — one click copies the "
                   "whole block (the grid above can't, and Ctrl+C there triggers a rerun).")
        st.caption("**Markdown** (Obsidian) — pastes as a real table.")
        st.code(_df_to_markdown(export), language="markdown")
        st.caption("**CSV** — comma-separated, for spreadsheets / scripts.")
        st.code(export.to_csv(index=False), language="text")
        st.caption("**TSV** (Word / Excel) — paste and it auto-formats as a table.")
        st.code(export.to_csv(sep="\t", index=False), language="text")
        if truncated:
            st.caption(f"⚠️ Showing the first {max_copy_rows:,} of {len(df):,} rows — "
                       "use the CSV download for the full table.")


def ibd_overlap_table(output_dir: Path) -> pd.DataFrame | None:
    f = output_dir / "ibd_overlap_outputs" / "ibd_target_overlap.csv"
    if not f.exists():
        return None
    return pd.read_csv(f)


def ibd_overlap_comparison(output_dir: Path) -> pd.DataFrame | None:
    """Validated vs Predicted vs Baseline IBD-overlap comparison for one dataset.

    Produced by 07_networks_clin/19_ibd_target_overlap.R after the all-miRNA
    baseline targets (06_targets_enrich/15_multimir_targets_baseline.R) exist.
    One row per (miRNA, evidence) with a shared per-evidence universe:
      evidence ∈ {validated, predicted}; is_significant flags the DE candidates.
    The three thesis layers are slices of this single table:
      validated     = is_significant & evidence == 'validated'
      non-validated = is_significant & evidence == 'predicted'
      baseline      = all miRNAs (the full per-evidence distribution)
    Returns None if the comparison hasn't been generated for this dataset.
    """
    f = output_dir / "ibd_overlap_outputs" / "ibd_target_overlap_comparison.csv"
    if not f.exists():
        return None
    return pd.read_csv(f)


def _shallowest_dir(parent: Path, name: str) -> Path | None:
    """Shallowest directory called `name` found under `parent`, else None.

    Mirrors the report generator's `find_dir` so the dashboard locates artifacts
    that live deep in the raw pipeline tree (e.g. clinical outputs nested under
    fastq/muestras_clasificadas/DEA_results_round_2/) the same way."""
    best, best_depth = None, None
    for root, dirs, _ in os.walk(parent):
        if name in dirs:
            cand = Path(root) / name
            depth = len(cand.parts)
            if best is None or depth < best_depth:
                best, best_depth = cand, depth
    return best


@lru_cache(maxsize=32)
def clinical_results(output_dir: Path) -> dict | None:
    """Treatment-response biomarker results for datasets with responder vs
    non-responder samples at baseline (07_networks_clin/17_treatment_response_roc.R).

    Most datasets are simple case/control and never run this step — for those this
    returns None (and the dashboard hides the Clinical tab). When present, returns
    the ROC tables and figure paths discovered in the dataset's clinical_outputs/
    folder. Cached because it walks the (large) raw output tree to find it.
    """
    clin_dir = _shallowest_dir(output_dir, "clinical_outputs")
    if clin_dir is None:
        return None

    def _csv(name: str) -> pd.DataFrame | None:
        f = clin_dir / name
        try:
            return pd.read_csv(f) if f.exists() else None
        except Exception:
            return None

    res = {
        "dir": clin_dir,
        "roc_rvsnr": _csv("roc_R_vs_NR.csv"),
        "roc_uchc": _csv("roc_UC_vs_HC.csv"),
        "cv_model": _csv("cv_model_R_vs_NR.csv"),
        "candidates": _csv("candidate_mirnas.csv"),
    }
    # The folder also holds IBD/HMDD files for every dataset, so only treat it as a
    # clinical layer when the ROC analysis itself actually produced something.
    if all(res[k] is None for k in ("roc_rvsnr", "roc_uchc", "candidates")):
        return None

    fig_dir = clin_dir / "figures"

    def _pngs(prefix: str) -> list[Path]:
        return sorted(fig_dir.glob(f"{prefix}*.png")) if fig_dir.is_dir() else []

    res["roc_rvsnr_pngs"] = _pngs("ROC_RvsNR_")
    res["roc_uchc_pngs"] = _pngs("ROC_UCvsHC_")
    res["violin_pngs"] = _pngs("violin_")
    multi = fig_dir / "ROC_multi_RvsNR.png"
    res["roc_multi"] = multi if multi.exists() else None
    return res


def ibd_comparison_long() -> pd.DataFrame:
    """Pool every analyzed dataset's ibd_target_overlap_comparison.csv into one
    long table with a leading 'dataset' column. One row per (dataset, miRNA,
    evidence). Drives the cross-dataset Validated/Predicted/Baseline view.
    Returns an empty frame if no dataset has the comparison generated yet.
    """
    frames = []
    for d in analyzed():
        tbl = ibd_overlap_comparison(resolve(d["output_dir"]))
        if tbl is None or tbl.empty:
            continue
        sub = tbl.copy()
        sub.insert(0, "dataset", d["id"])
        frames.append(sub)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    if "is_significant" in out.columns:
        out["is_significant"] = out["is_significant"].astype(bool)
    return out


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


# ── Card helpers ──────────────────────────────────────────────────────────────

def to_file_url(path: Path) -> str:
    """Absolute path → file:// URL, translating WSL /mnt/<drive>/ prefixes."""
    s = str(path)
    if sys.platform != "win32" and s.startswith("/mnt/") and len(s) > 6 and s[5] != "/":
        drive = s[5].upper()
        rest = s[6:].lstrip("/")
        return f"file:///{drive}:/{rest}"
    return path.as_uri()


@lru_cache(maxsize=256)
def _img_data_uri(path: Path) -> str:
    """base64 data-URI for a local PNG so it can be embedded in inline HTML
    (cached — galleries re-render on every interaction)."""
    return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()


def _lightbox_html(uri: str, caption: str | None, uid: str) -> str:
    """Inline CSS-only lightbox: click the image to enlarge it to full screen in
    place, click again (anywhere) to close. Uses the :target pseudo-class driven by
    a #hash anchor — the same mechanism as Streamlit's own header anchors — so it
    works in Chrome (which blocks opening data: URIs in a new tab) with no JS."""
    cap = (f'<div style="font-size:.8rem;color:#9aa0a6;text-align:center;'
           f'margin-top:4px">{caption}</div>') if caption else ""
    return f"""
<style>
#fig_{uid}{{position:relative;margin:0}}
#fig_{uid} .zoomi{{width:100%;border-radius:6px;display:block;cursor:zoom-in}}
#fig_{uid} .zopen{{position:absolute;top:0;left:0;right:0;bottom:0}}
#fig_{uid} .zclose{{display:none}}
#{uid}:target{{position:fixed;top:0;left:0;right:0;bottom:0;margin:auto;
  width:auto;height:auto;max-width:96vw;max-height:96vh;z-index:99990;
  cursor:zoom-out;box-shadow:0 0 0 100vmax rgba(0,0,0,.92)}}
#{uid}:target ~ .zclose{{display:block;position:fixed;top:0;left:0;right:0;
  bottom:0;z-index:99991;cursor:zoom-out}}
</style>
<div id="fig_{uid}">
  <img class="zoomi" id="{uid}" src="{uri}"/>
  <a class="zopen" href="#{uid}" title="Click to zoom"></a>
  <a class="zclose" href="#_" title="Close"></a>
  {cap}
</div>
"""


def zoomable_image(path: Path, caption: str | None = None) -> None:
    """Render an image that enlarges to full screen on click (CSS lightbox) with an
    optional caption — for the dense ROC / volcano / QC / Venn figures."""
    import hashlib
    import streamlit as st
    uri = _img_data_uri(Path(path))
    uid = "z" + hashlib.md5(str(path).encode()).hexdigest()[:10]
    st.markdown(_lightbox_html(uri, caption, uid), unsafe_allow_html=True)


def zoomable_pyplot(fig, caption: str | None = None) -> None:
    """Click-to-zoom wrapper for an in-memory matplotlib figure (Venn / UpSet):
    renders it to a PNG, then shows it in the same CSS lightbox as zoomable_image."""
    import hashlib
    import io
    import streamlit as st
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    data = buf.getvalue()
    uri = "data:image/png;base64," + base64.b64encode(data).decode()
    uid = "z" + hashlib.md5(data[:4096]).hexdigest()[:10]
    st.markdown(_lightbox_html(uri, caption, uid), unsafe_allow_html=True)


def open_in_os(path: Path) -> None:
    """Open a file in its default OS application.

    In WSL the browser blocks file:// links, so we call explorer.exe directly —
    it is always reachable from WSL2 and respects Windows file associations.
    """
    import subprocess
    s = str(path)
    if sys.platform != "win32" and s.startswith("/mnt/") and len(s) > 6 and s[5] != "/":
        drive = s[5].upper()
        rest = s[6:].lstrip("/").replace("/", "\\")
        subprocess.Popen(["explorer.exe", f"{drive}:\\{rest}"])
    elif sys.platform == "win32":
        subprocess.Popen(["explorer.exe", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


@lru_cache(maxsize=32)
def find_qc_reports(output_dir: Path) -> dict[str, Path]:
    """Locate MultiQC and miRTrace HTML reports under output_dir (cached)."""
    out: dict[str, Path] = {}
    hits = list(output_dir.glob("**/multiqc_report.html"))
    if hits:
        out["MultiQC"] = hits[0]
    hits = list(output_dir.glob("**/mirtrace-report.html"))
    if hits:
        out["miRTrace"] = hits[0]
    return out


@lru_cache(maxsize=32)
def card_dea_stats(ds_id: str, padj: float = 0.05, lfc: float = 0.58) -> dict | None:
    """Quick DE summary for a card: total unique significant miRNAs (cached)."""
    d = get_dataset(ds_id)
    if not d or not d.get("output_dir"):
        return None
    odir = resolve(d["output_dir"])
    dea = discover_dea(odir)
    if not dea:
        return None
    mode = "Default" if "Default" in dea else next(iter(dea))
    total_sig: set[str] = set()
    per_comp: dict[str, dict] = {}
    for comp, csv in dea[mode].items():
        try:
            df = dea_significant(load_dea(csv), padj, lfc)
            df = add_direction(df)
            per_comp[comp] = {
                "up": int((df["direction"] == "UP").sum()),
                "down": int((df["direction"] == "DOWN").sum()),
            }
            total_sig.update(df["miRNA"].tolist())
        except Exception:
            continue
    return {"n_sig": len(total_sig), "per_comparison": per_comp}


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


def ibd_overlap_long(metric: str = "padj", cutoff: float = 0.05,
                     enriched_only: bool = True) -> pd.DataFrame:
    """IBD target-overlap rows pooled across all analyzed datasets.

    One row per (dataset, miRNA) from each dataset's ibd_target_overlap.csv.
    `metric` ∈ {'padj','p_fisher'} is filtered at <= cutoff; `enriched_only`
    additionally keeps OR > 1. A leading 'dataset' column is added. This is the
    validated-target IBD layer — independent of the DE significance thresholds.
    """
    frames = []
    for d in analyzed():
        tbl = ibd_overlap_table(resolve(d["output_dir"]))
        if tbl is None or "miRNA" not in tbl.columns:
            continue
        sub = tbl.copy()
        if metric in sub.columns:
            sub = sub[pd.to_numeric(sub[metric], errors="coerce") <= cutoff]
        if enriched_only and "OR" in sub.columns:
            sub = sub[pd.to_numeric(sub["OR"], errors="coerce") > 1]
        if sub.empty:
            continue
        sub.insert(0, "dataset", d["id"])
        frames.append(sub)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def overlap_regions(named: "dict[str, set]") -> "dict[str, list]":
    """Exclusive Venn regions for 2 or 3 named sets, mirroring the diagram.

    Returns an ordered {region_label: sorted member list}. For 3 sets this yields
    the 3 'only', 3 pairwise-exclusive, and 1 triple-intersection regions.
    """
    names = list(named)
    out: "dict[str, list]" = {}
    if len(names) == 2:
        a, b = names
        A, B = named[a], named[b]
        out[f"{a} only"] = sorted(A - B)
        out[f"{b} only"] = sorted(B - A)
        out[f"{a} ∩ {b}"] = sorted(A & B)
    elif len(names) == 3:
        a, b, c = names
        A, B, C = named[a], named[b], named[c]
        out[f"{a} only"] = sorted(A - B - C)
        out[f"{b} only"] = sorted(B - A - C)
        out[f"{c} only"] = sorted(C - A - B)
        out[f"{a} ∩ {b}"] = sorted((A & B) - C)
        out[f"{a} ∩ {c}"] = sorted((A & C) - B)
        out[f"{b} ∩ {c}"] = sorted((B & C) - A)
        out[f"{a} ∩ {b} ∩ {c}"] = sorted(A & B & C)
    return out


def exclusive_intersections(named: "dict[str, set]") -> "dict[str, list]":
    """Exclusive (UpSet-style) intersections for any number of named sets.

    For every non-empty combination of datasets, returns the members that are DE
    in exactly those datasets and none of the others — i.e. one entry per UpSet
    bar. Ordered by descending size to mirror UpSet's `sort_by='cardinality'`.
    """
    from itertools import combinations
    names = list(named)
    rows: list[tuple[str, list]] = []
    for r in range(len(names), 0, -1):
        for combo in combinations(names, r):
            inside = set.intersection(*(named[n] for n in combo))
            others = [named[n] for n in names if n not in combo]
            outside = set().union(*others) if others else set()
            members = sorted(inside - outside)
            if members:
                label = f"{combo[0]} only" if r == 1 else " ∩ ".join(combo)
                rows.append((label, members))
    rows.sort(key=lambda kv: len(kv[1]), reverse=True)
    return dict(rows)


@lru_cache(maxsize=1)
def de_universe() -> tuple:
    """Sorted miRNAs DE (padj≤0.05, |log2FC|≥0.58) in ≥1 dataset/comparison, across
    either replacement mode. Drives the miRNA Lookup dropdown."""
    seen: set = set()
    for mode in ("Default", "Strict"):
        df = collect_de_mirnas(0.05, 0.58, mode=mode)
        if not df.empty:
            seen.update(df["miRNA"].astype(str))
    return tuple(sorted(seen, key=str.lower))


def _norm_mirna(name) -> str:
    return str(name).strip().lower()


def mirna_dea_profile(mirna: str) -> pd.DataFrame:
    """Every (dataset, mode, comparison) DESeq2 result for one miRNA.

    Columns: dataset, tissue, mode, comparison, log2FoldChange, padj, pvalue,
    baseMean, direction, significant — one row per place the miRNA appears."""
    key = _norm_mirna(mirna)
    rows = []
    for d in analyzed():
        dea = discover_dea(resolve(d["output_dir"]))
        for mode, comp_map in dea.items():
            for comparison, csv in comp_map.items():
                try:
                    df = load_dea(csv)
                except Exception:
                    continue
                hit = df[df["miRNA"].astype(str).str.lower() == key]
                if hit.empty:
                    continue
                r = hit.iloc[0]
                lfc, padj = r.get("log2FoldChange"), r.get("padj")
                rows.append({
                    "dataset": d["id"], "tissue": d.get("tissue"),
                    "mode": mode, "comparison": comparison,
                    "log2FoldChange": lfc, "padj": padj,
                    "pvalue": r.get("pvalue"), "baseMean": r.get("baseMean"),
                    "direction": (("UP" if lfc >= 0 else "DOWN")
                                  if pd.notna(lfc) else None),
                    "significant": bool(pd.notna(padj) and padj <= 0.05
                                        and pd.notna(lfc) and abs(lfc) >= 0.58),
                })
    return pd.DataFrame(rows)


def mirna_ibd_profile(mirna: str) -> pd.DataFrame:
    """Per-dataset IBD target-overlap row(s) for one miRNA (with a tissue column)."""
    key = _norm_mirna(mirna)
    frames = []
    for d in analyzed():
        tbl = ibd_overlap_table(resolve(d["output_dir"]))
        if tbl is None or "miRNA" not in tbl.columns:
            continue
        hit = tbl[tbl["miRNA"].astype(str).str.lower() == key].copy()
        if hit.empty:
            continue
        hit.insert(0, "tissue", d.get("tissue"))
        hit.insert(0, "dataset", d["id"])
        frames.append(hit)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


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
