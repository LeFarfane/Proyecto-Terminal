#!/usr/bin/env python3
"""
03_geo_id_curator.py  —  The Targeted GSE Edition

Fetches metadata and generates an outputs/docs/{GSE_ID}.md file for one or
more explicit GEO accessions (GSE...).  Mirrors the output format of
01_geo_ibd_curator.py, but skips the broad query step and goes straight to
the accessions you name.

Usage:
    python 03_geo_id_curator.py --id GSE114603
    python 03_geo_id_curator.py --input my_gse_ids.txt
    python 03_geo_id_curator.py --id GSE114603 --out ../../outputs
"""

import argparse
import csv
import io
import os
import re
import subprocess
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


KEEP_COLS = [
    "Run", "avgLength", "size_MB", "LibraryStrategy", "LibrarySelection",
    "LibrarySource", "LibraryLayout", "Platform", "Model", "BioProject",
    "BioSample", "ScientificName",
]


# ── Shell helpers ──────────────────────────────────────────────────────────────

def run_local_cmd(cmd: List[str], input_text: Optional[str] = None) -> str:
    print(f"\n🚀 [CMD]: {' '.join(cmd)}")
    try:
        p = subprocess.run(
            cmd,
            input=input_text.encode("utf-8") if input_text is not None else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if p.returncode != 0:
            print(f"❌ [STDERR]:\n{p.stderr.decode('utf-8', errors='replace')}")
        out = p.stdout.decode("utf-8", errors="replace")
        preview = out.strip()
        if len(preview) > 500:
            preview = preview[:500] + "\n... [OUTPUT TRUNCATED]"
        print(f"📥 [STDOUT PREVIEW]:\n{preview if preview else '<EMPTY>'}")
        return out
    except FileNotFoundError:
        raise RuntimeError(f"Command not found: {cmd[0]}. Ensure EDirect is in your PATH.")
    except Exception as e:
        print(f"💥 [EXCEPTION]: {e}")
        return ""


def run_network_cmd(cmd: List[str], delay: float, input_text: Optional[str] = None) -> str:
    time.sleep(delay)
    return run_local_cmd(cmd, input_text)


# ── GEO fetch ─────────────────────────────────────────────────────────────────

def fetch_gse_entry(gse_id: str, delay: float) -> Optional[Dict[str, str]]:
    """Query GEO by exact accession and return {GSE_ID, Title, Summary, PMID}."""
    print(f"\n🔍 Fetching GEO entry for {gse_id}...")
    search_out = run_network_cmd(
        ["esearch", "-db", "gds", "-query", f"{gse_id}[ACCN]"], delay
    )
    docsum_xml = run_network_cmd(
        ["efetch", "-db", "gds", "-format", "docsum"], delay, input_text=search_out
    )

    try:
        root = ET.fromstring(docsum_xml)
    except ET.ParseError as e:
        print(f"⚠️  XML parse error: {e}")
        return None

    for docsum in root.findall(".//DocumentSummary"):
        acc = docsum.findtext("Accession", default="").strip()
        if acc.upper() != gse_id.upper():
            continue
        pmid = ""
        pmid_node = docsum.find(".//PubMedIds/int")
        if pmid_node is not None and pmid_node.text:
            pmid = pmid_node.text.strip()
        return {
            "GSE_ID": acc,
            "Title": docsum.findtext("title", default="").strip(),
            "Summary": docsum.findtext("summary", default="").strip(),
            "PMID": pmid,
        }

    print(f"⚠️  No GEO entry found for {gse_id}.")
    return None


# ── SRA metadata (identical to 01_geo_ibd_curator) ───────────────────────────

def get_rich_sra_metadata(gse_id: str, delay: float) -> Dict[str, Any]:
    meta: Dict[str, Any] = {
        "available": False, "srr_count": 0, "total_mb": 0.0, "total_gb": 0.0,
        "bioprojects": set(), "strategies": set(), "layouts": set(), "platforms": set(),
        "diseases": set(), "tissues": set(),
        "merged_csv_data": [], "merged_csv_headers": [],
    }

    try:
        # STEP 1: GSE → BioProject
        print("\n--- STEP 1: GSE → BioProject ---")
        search_bp = run_network_cmd(["esearch", "-db", "bioproject", "-query", gse_id], delay)
        docsum_bp = run_network_cmd(["efetch", "-format", "docsum"], delay, input_text=search_bp)
        xtract_bp = run_local_cmd(
            ["xtract", "-pattern", "DocumentSummary", "-element", "Project_Acc"],
            input_text=docsum_bp,
        )
        bp_accs = [a.strip() for a in xtract_bp.splitlines() if a.strip().startswith("PRJ")]

        if not bp_accs:
            print(f"⚠️  Direct BioProject search failed for {gse_id}. Trying elink fallback...")
            link_bp = run_network_cmd(
                ["elink", "-db", "gds", "-id", gse_id, "-target", "bioproject"], delay
            )
            docsum_bp = run_network_cmd(["efetch", "-format", "docsum"], delay, input_text=link_bp)
            xtract_bp = run_local_cmd(
                ["xtract", "-pattern", "DocumentSummary", "-element", "Project_Acc"],
                input_text=docsum_bp,
            )
            bp_accs = [a.strip() for a in xtract_bp.splitlines() if a.strip().startswith("PRJ")]

        print(f"🎯 BioProjects found: {bp_accs}")
        if not bp_accs:
            print(f"🛑 No BioProjects found for {gse_id}.")
            return meta

        meta["bioprojects"].update(bp_accs)

        # STEP 2: BioProject → SRA RunInfo
        print("\n--- STEP 2: BioProject → SRA RunInfo ---")
        runs: List[Dict] = []
        biosample_ids: set = set()

        for bp_acc in bp_accs:
            search_bp_sra = run_network_cmd(
                ["esearch", "-db", "bioproject", "-query", bp_acc], delay
            )
            link_sra = run_network_cmd(["elink", "-target", "sra"], delay, input_text=search_bp_sra)
            runinfo_raw = run_network_cmd(["efetch", "-format", "runinfo"], delay, input_text=link_sra)

            if not runinfo_raw.strip():
                continue
            for row in csv.DictReader(io.StringIO(runinfo_raw)):
                if not row.get("Run", "").startswith("SRR"):
                    continue
                runs.append({k: row.get(k, "") for k in KEEP_COLS if k in row})
                if row.get("BioSample", "").startswith("SAMN"):
                    biosample_ids.add(row["BioSample"])
                try:
                    meta["total_mb"] += float(row.get("size_MB") or 0)
                except ValueError:
                    pass
                meta["strategies"].add(row.get("LibraryStrategy", ""))
                meta["layouts"].add(row.get("LibraryLayout", ""))

        print(f"🎯 Extracted {len(runs)} runs and {len(biosample_ids)} BioSample IDs.")
        if not runs:
            print(f"🛑 No SRA runs found for {gse_id}.")
            return meta

        meta["available"] = True
        meta["srr_count"] = len(runs)
        meta["total_gb"] = meta["total_mb"] / 1024

        # STEP 3: BioSample attributes
        print("\n--- STEP 3: BioSample Attributes ---")
        biosample_attrs: Dict[str, Dict[str, str]] = {}
        all_attr_keys: set = set()

        for i in range(0, len(list(biosample_ids)), 200):
            chunk = list(biosample_ids)[i:i + 200]
            print(f"📦 Fetching BioSample batch {i + 1}–{i + len(chunk)}...")
            bs_xml = run_network_cmd(
                ["efetch", "-db", "biosample", "-id", ",".join(chunk), "-format", "xml"], delay
            )
            root = ET.fromstring(bs_xml)
            for sample in root.findall(".//BioSample"):
                acc = sample.get("accession")
                if not acc:
                    continue
                attrs: Dict[str, str] = {}
                for attr in sample.findall(".//Attribute"):
                    k, v = attr.get("attribute_name"), attr.text
                    if k and v:
                        clean_k = k.strip().replace(" ", "_")
                        attrs[clean_k] = v.strip()
                        all_attr_keys.add(clean_k)
                        kl = k.lower()
                        if any(w in kl for w in ("disease", "phenotype", "condition")):
                            meta["diseases"].add(v.strip())
                        if any(w in kl for w in ("tissue", "source", "body", "cell type")):
                            meta["tissues"].add(v.strip())
                biosample_attrs[acc] = attrs

        # STEP 4: Merge
        print("\n--- STEP 4: Final Merge ---")
        sorted_attr_keys = sorted(all_attr_keys)
        meta["merged_csv_headers"] = list(KEEP_COLS) + sorted_attr_keys
        for run in runs:
            bs_data = biosample_attrs.get(run.get("BioSample", ""), {})
            meta["merged_csv_data"].append(
                {**run, **{k: bs_data.get(k, "") for k in sorted_attr_keys}}
            )

        for key in ("bioprojects", "strategies", "layouts", "platforms", "diseases", "tissues"):
            meta[key] = {x for x in meta[key] if x}

        print("✅ Merge successful.")
        return meta

    except Exception as e:
        print(f"\n⚠️  Pipeline error for {gse_id}: {e}")
        return meta


# ── Markdown writer (identical format to 01_geo_ibd_curator) ─────────────────

def normalize_text(s: str) -> str:
    if not s:
        return ""
    return re.sub(r"\n{3,}", "\n\n", s.replace("\r\n", "\n").replace("\r", "\n")).strip()


def write_markdown_doc(docs_dir: Path, entry: Dict[str, Any]) -> None:
    gse = entry["GSE_ID"]
    meta = entry.get("rich_meta", {})

    bioprojects = ", ".join(meta.get("bioprojects", [])) or "N/A"
    strategies  = ", ".join(meta.get("strategies",  [])) or "N/A"
    layouts     = ", ".join(meta.get("layouts",     [])) or "N/A"
    diseases    = ", ".join(list(meta.get("diseases", []))[:10]) or "None extracted"
    tissues     = ", ".join(list(meta.get("tissues",  []))[:10]) or "None extracted"

    run_table  = "### Run Information Overview\n\n"
    run_table += "| Run | Size (MB) | Strategy | Layout | Disease / Condition | Tissue / Source |\n"
    run_table += "|---|---|---|---|---|---|\n"
    if meta.get("merged_csv_data"):
        for row in meta["merged_csv_data"]:
            disease_val = row.get("disease", row.get("disease_state", row.get("phenotype", "N/A")))
            tissue_val  = row.get("tissue",  row.get("source_name",   row.get("body_site",  "N/A")))
            run_table += (
                f"| {row.get('Run','N/A')} | {row.get('size_MB','0')} "
                f"| {row.get('LibraryStrategy','N/A')} | {row.get('LibraryLayout','N/A')} "
                f"| {disease_val} | {tissue_val} |\n"
            )
    else:
        run_table += "| N/A | N/A | N/A | N/A | N/A | N/A |\n"

    md = f"""# {gse}

## Identification
- **Title:** {entry.get("Title", "").strip()}
- **PMID:** {entry.get("PMID", "").strip()}
- **GEO URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}
- **BioProject(s):** {bioprojects}

## Sequencing Metadata (SRA)
- **SRA Available:** {'Yes' if meta.get('available') else 'No'}
- **Total Runs:** {meta.get('srr_count', 0)}
- **Total Payload:** {meta.get('total_gb', 0):.2f} GB ({meta.get('total_mb', 0):,.0f} MB)
- **Library Strategy:** {strategies}
- **Library Layout:** {layouts}

## BioSample Characteristics
- **Conditions/Diseases:** {diseases}
- **Tissues/Sources:** {tissues}

---
{run_table}
---

## PT Evaluation Notes
- **Tag:** #to_review
- **Overall priority:**
- **Bioinformatics Notes:**

*Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}*

## Abstract/Summary
{normalize_text(entry.get("Summary", ""))}
"""
    out_path = docs_dir / f"{gse}.md"
    out_path.write_text(md, encoding="utf-8")
    print(f"📄 Written: {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def load_id_file(filepath: Path) -> List[str]:
    ids = []
    with filepath.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip().upper()
            if line and not line.startswith("#"):
                ids.append(line)
    return ids


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Targeted curator for explicit GEO accessions (GSE...)."
    )
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--id",    help="Single GSE accession, e.g. GSE114603")
    group.add_argument("--input", help="Text file with one GSE accession per line")
    ap.add_argument("--out", default="../../outputs", help="Master output folder")
    args = ap.parse_args()

    print("\n" + "=" * 55)
    print("👋 GEO ID Curator  —  Targeted GSE Edition")
    print("=" * 55)

    gse_ids: List[str] = (
        [args.id.strip().upper()]
        if args.id
        else load_id_file(Path(args.input).resolve())
    )

    if not gse_ids:
        print("❌ No GSE IDs provided.")
        return 1

    api_key = os.environ.get("NCBI_API_KEY")
    if not api_key:
        print("\n⚠️  No NCBI_API_KEY found in environment.")
        user_key = input("🔑 Enter NCBI API key (or Enter for 3 req/sec):\n> ").strip()
        if user_key:
            os.environ["NCBI_API_KEY"] = user_key
            api_key = user_key
    delay = 0.11 if api_key else 0.34

    out_dir  = Path(args.out).resolve()
    docs_dir = out_dir / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)

    for i, gse_id in enumerate(gse_ids, 1):
        print("\n" + "=" * 60)
        print(f"🧬 PROCESSING {i}/{len(gse_ids)}: {gse_id}")
        print("=" * 60)

        entry = fetch_gse_entry(gse_id, delay)
        if entry is None:
            print(f"⚠️  Skipping {gse_id} — could not fetch GEO entry.")
            continue

        entry["rich_meta"] = get_rich_sra_metadata(gse_id, delay)

        if entry["rich_meta"].get("available") and entry["rich_meta"].get("merged_csv_data"):
            gse_dir = out_dir / gse_id
            gse_dir.mkdir(parents=True, exist_ok=True)
            csv_path = gse_dir / f"{gse_id}_runinfo_merged.csv"
            with csv_path.open("w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=entry["rich_meta"]["merged_csv_headers"])
                w.writeheader()
                w.writerows(entry["rich_meta"]["merged_csv_data"])
            srrs = [r["Run"] for r in entry["rich_meta"]["merged_csv_data"] if "Run" in r]
            (gse_dir / "srr_list.txt").write_text("\n".join(srrs), encoding="utf-8")

        write_markdown_doc(docs_dir, entry)

    print("\n\n🎉 Done. Docs written to:", docs_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
