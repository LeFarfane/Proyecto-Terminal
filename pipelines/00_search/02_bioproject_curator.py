#!/usr/bin/env python3
"""
bioproject_curator.py (The Sniper Edition)

Description:
    A targeted curation script that takes either a single explicit BioProject ID
    (PRJNA...) or a text file of IDs. It performs the SRA extraction, BioSample
    attribute merging, and Markdown generation, completely bypassing the GEO database.

Usage:
    ./bioproject_curator.py --input manual_targets.txt
    ./bioproject_curator.py --id PRJNA904924
"""

import argparse
import csv
import io
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Set

KEEP_COLS = [
    "Run", "avgLength", "size_MB", "LibraryStrategy", "LibrarySelection",
    "LibrarySource", "LibraryLayout", "Platform", "Model", "BioProject",
    "BioSample", "ScientificName"
]


def run_local_cmd(cmd: List[str], input_text: Optional[str] = None) -> str:
    """Executes a local shell command with aggressive logging."""
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
            preview = preview[:500] + "\n... [OUTPUT TRUNCATED FOR READABILITY]"
        print(f"📥 [STDOUT PREVIEW]:\n{preview if preview else '<EMPTY>'}")

        return out
    except FileNotFoundError:
        raise RuntimeError(f"Command not found: {cmd[0]}. Please ensure EDirect is in your PATH.")
    except Exception as e:
        print(f"💥 [EXCEPTION]: {e}")
        return ""


def run_network_cmd(cmd: List[str], delay: float, input_text: Optional[str] = None) -> str:
    """Executes a shell command that interacts with NCBI servers, enforcing a rate limit."""
    time.sleep(delay)
    return run_local_cmd(cmd, input_text)


def load_input_list(filepath: Path) -> List[str]:
    """Reads the target BioProjects from the user's input text file."""
    ids = []
    if filepath.exists():
        with filepath.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip().upper()
                if line and not line.startswith("#"):
                    ids.append(line)
    return ids


def load_ignore_list(filepath: Path) -> Set[str]:
    """Reads the ignore list to permanently bypass rejected datasets."""
    ignore_set = set()
    if filepath.exists():
        with filepath.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip().upper()
                if line and not line.startswith("#"):
                    ignore_set.add(line)
    return ignore_set


def get_bioproject_summary(bp_acc: str, delay: float) -> Dict[str, str]:
    """Fetches the Title and Description directly from the BioProject database."""
    print(f"\n⚙️  Fetching BioProject Summary for {bp_acc}...")
    search_bp = run_network_cmd(["esearch", "-db", "bioproject", "-query", bp_acc], delay)
    docsum_bp = run_network_cmd(["efetch", "-format", "docsum"], delay, input_text=search_bp)

    summary_data = {"Title": "N/A", "Description": "No description available."}
    try:
        root = ET.fromstring(docsum_bp)
        for docsum in root.findall('.//DocumentSummary'):
            title = docsum.findtext('Project_Title')
            if title:
                summary_data["Title"] = title.strip()

            desc = docsum.findtext('Project_Description')
            if desc:
                summary_data["Description"] = desc.strip()
            break
    except ET.ParseError as e:
        print(f"⚠️ XML Parsing Error for BioProject: {e}")

    return summary_data


def get_rich_sra_metadata(bp_acc: str, delay: float) -> Dict[str, Any]:
    """
    Direct SRA Pivot: Skips GEO entirely and directly links the BioProject to SRA.
    """
    meta = {
        "available": False, "srr_count": 0, "total_mb": 0.0, "total_gb": 0.0,
        "strategies": set(), "layouts": set(), "platforms": set(),
        "diseases": set(), "tissues": set(),
        "merged_csv_data": [], "merged_csv_headers": []
    }

    try:
        # STEP 1: BioProject -> SRA RunInfo
        print("\n--- STEP 1: BioProject -> SRA RunInfo ---")
        runs = []
        biosample_ids = set()

        search_bp_sra = run_network_cmd(["esearch", "-db", "bioproject", "-query", bp_acc], delay)
        link_sra = run_network_cmd(["elink", "-target", "sra"], delay, input_text=search_bp_sra)
        runinfo_raw = run_network_cmd(["efetch", "-format", "runinfo"], delay, input_text=link_sra)

        if runinfo_raw.strip():
            reader = csv.DictReader(io.StringIO(runinfo_raw))
            for row in reader:
                if not row.get("Run", "").startswith("SRR"):
                    continue
                filtered_row = {k: row.get(k, "") for k in KEEP_COLS if k in row}
                runs.append(filtered_row)

                if "BioSample" in row and row["BioSample"].startswith("SAMN"):
                    biosample_ids.add(row["BioSample"])

                try:
                    meta["total_mb"] += float(row.get("size_MB") or 0)
                except ValueError:
                    pass
                meta["strategies"].add(row.get("LibraryStrategy", ""))
                meta["layouts"].add(row.get("LibraryLayout", ""))

        print(f"🎯 [RESULT]: Extracted {len(runs)} SRR runs and {len(biosample_ids)} BioSample IDs.")

        if not runs:
            print(f"🛑 [STOP]: No SRA runs found for BioProject: {bp_acc}")
            return meta

        meta["available"] = True
        meta["srr_count"] = len(runs)
        meta["total_gb"] = meta["total_mb"] / 1024

        # STEP 2: BioSample Attributes (Batched natively)
        print("\n--- STEP 2: BioSample Attributes ---")
        biosample_attrs: Dict[str, Dict[str, str]] = {}
        all_attr_keys = set()
        bs_list = list(biosample_ids)

        chunk_size = 200
        for i in range(0, len(bs_list), chunk_size):
            chunk = bs_list[i:i + chunk_size]
            print(f"📦 Fetching BioSample batch {i + 1} to {i + len(chunk)}...")
            bs_xml = run_network_cmd(["efetch", "-db", "biosample", "-id", ",".join(chunk), "-format", "xml"], delay)

            root = ET.fromstring(bs_xml)
            for sample in root.findall('.//BioSample'):
                acc = sample.get('accession')
                if not acc:
                    continue
                attrs = {}
                for attr in sample.findall('.//Attribute'):
                    k = attr.get('attribute_name')
                    v = attr.text
                    if k and v:
                        clean_k = k.strip().replace(" ", "_")
                        attrs[clean_k] = v.strip()
                        all_attr_keys.add(clean_k)

                        k_lower = k.lower()
                        if "disease" in k_lower or "phenotype" in k_lower or "condition" in k_lower:
                            meta["diseases"].add(v.strip())
                        if "tissue" in k_lower or "source" in k_lower or "body" in k_lower or "cell type" in k_lower:
                            meta["tissues"].add(v.strip())
                biosample_attrs[acc] = attrs

        # STEP 3: Final Merge
        print("\n--- STEP 3: Final Metadata Merge ---")
        sorted_attr_keys = sorted(list(all_attr_keys))
        meta["merged_csv_headers"] = [k for k in KEEP_COLS] + sorted_attr_keys

        for run in runs:
            bs_id = run.get("BioSample", "")
            bs_data = biosample_attrs.get(bs_id, {})
            merged_row = run.copy()
            for key in sorted_attr_keys:
                merged_row[key] = bs_data.get(key, "")
            meta["merged_csv_data"].append(merged_row)

        for key in ["strategies", "layouts", "platforms", "diseases", "tissues"]:
            meta[key] = {item for item in meta[key] if item}

        print("✅ Merge successful.")
        return meta

    except Exception as e:
        print(f"\n⚠️ Pipeline Error for {bp_acc}: {e}")
        return meta


def write_markdown_doc(docs_dir: Path, bp_acc: str, bp_summary: Dict[str, str], meta: Dict[str, Any]) -> None:
    """Generates a rich Markdown file anchored around a BioProject ID."""
    platforms = ", ".join(meta.get("platforms", [])) or "N/A"
    strategies = ", ".join(meta.get("strategies", [])) or "N/A"
    layouts = ", ".join(meta.get("layouts", [])) or "N/A"
    diseases = ", ".join(list(meta.get("diseases", []))[:10]) or "None extracted"
    tissues = ", ".join(list(meta.get("tissues", []))[:10]) or "None extracted"

    run_table = "### Run Information Overview\n\n"
    run_table += "| Run | Size (MB) | Strategy | Layout | Disease / Condition | Tissue / Source |\n"
    run_table += "|---|---|---|---|---|---|\n"

    if meta.get("merged_csv_data"):
        for row in meta["merged_csv_data"]:
            disease_val = row.get("disease", row.get("disease_state", row.get("phenotype", "N/A")))
            tissue_val = row.get("tissue", row.get("source_name", row.get("body_site", "N/A")))
            run_table += f"| {row.get('Run', 'N/A')} | {row.get('size_MB', '0')} | {row.get('LibraryStrategy', 'N/A')} | {row.get('LibraryLayout', 'N/A')} | {disease_val} | {tissue_val} |\n"
    else:
        run_table += "| N/A | N/A | N/A | N/A | N/A | N/A |\n"

    md = f"""# {bp_acc}

## Identification
- **Project Title:** {bp_summary.get("Title", "N/A")}
- **NCBI URL:** https://www.ncbi.nlm.nih.gov/bioproject/{bp_acc}

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

## Abstract/Project Description
{bp_summary.get("Description", "N/A")}
"""
    (docs_dir / f"{bp_acc}.md").write_text(md, encoding="utf-8")


def write_run_metadata(log_path: Path, input_source: str, total_processed: int, has_key: bool):
    """Generates a master text file recording the execution."""
    log = f"""---
Targeted BioProject Curator Execution Log
Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
---

## Execution Parameters
- Input Source: {input_source}
- API Rate Limit: {'10 requests/sec (Key Active)' if has_key else '3 requests/sec (No Key)'}

## Results
- Total BioProjects Processed: {total_processed}
"""
    log_path.write_text(log, encoding="utf-8")


def main() -> int:
    try:
        ap = argparse.ArgumentParser(description="Targeted curator for explicit BioProject IDs.")

        # Create a mutually exclusive group so the user must provide exactly one of these options
        group = ap.add_mutually_exclusive_group(required=True)
        group.add_argument("--input", help="Path to a text file containing BioProject IDs (PRJNA...)")
        group.add_argument("--id", help="A single BioProject ID to process directly (e.g., PRJNA904924)")

        # CHANGED: Default is now two directories up -> /Pipelines/outputs
        ap.add_argument("--out", default="../../outputs", help="Master output folder")
        ap.add_argument("--ignore-list", default="ignore_list.txt",
                        help="Filename for the ignore list (saved inside output folder)")
        args = ap.parse_args()

        print("\n" + "=" * 50)
        print("👋 Welcome to the Targeted BioProject Curator (The Sniper Edition)")
        print("=" * 50)

        # Handle the dual-input logic
        target_bp_accs = []
        log_source_name = ""

        if args.input:
            input_path = Path(args.input).resolve()
            if not input_path.exists():
                print(f"❌ Error: Input file '{args.input}' not found.")
                return 1
            target_bp_accs = load_input_list(input_path)
            log_source_name = f"File: {args.input}"
        elif args.id:
            target_bp_accs = [args.id.strip().upper()]
            log_source_name = f"Single ID: {args.id}"

        if not target_bp_accs:
            print(f"❌ No valid BioProject IDs provided.")
            return 1

        api_key = os.environ.get("NCBI_API_KEY")
        if not api_key:
            print("\n⚠️  No NCBI_API_KEY detected in your system environment.")
            user_key = input("🔑 Please enter your NCBI API key (or press Enter for 3 req/sec):\n> ").strip()
            if user_key:
                os.environ["NCBI_API_KEY"] = user_key
                api_key = user_key
        delay_seconds = 0.11 if api_key else 0.34

        # Directory Setup
        out_dir = Path(args.out).resolve()
        docs_dir = out_dir / "docs"
        out_dir.mkdir(parents=True, exist_ok=True)
        docs_dir.mkdir(parents=True, exist_ok=True)

        ignore_filepath = out_dir / args.ignore_list
        failed_sra_log = out_dir / "failed_sra_datasets.txt"

        if not ignore_filepath.exists():
            ignore_filepath.write_text("# Add BioProject/GSE IDs here to skip them in future runs (one per line)\n",
                                       encoding="utf-8")

        ignore_set = load_ignore_list(ignore_filepath)
        if ignore_set:
            print(f"\n🛡️  Loaded {len(ignore_set)} datasets to ignore from {ignore_filepath.name}")

        valid_entries = []
        ignored_count = 0
        for bp in target_bp_accs:
            if bp in ignore_set:
                ignored_count += 1
                continue
            valid_entries.append(bp)

        if ignored_count > 0:
            print(f"\n🚫 Bypassed {ignored_count} BioProjects matched against your ignore list.")

        if not valid_entries:
            print("\n✅ All provided BioProjects were successfully bypassed by your ignore list.")
            return 0

        print(f"\n📊 Processing {len(valid_entries)} Targeted BioProjects...")

        for i, bp_acc in enumerate(valid_entries, 1):

            print("\n" + "=" * 60)
            print(f"🧬 PROCESSING BIOPROJECT {i}/{len(valid_entries)}: {bp_acc}")
            print("=" * 60)

            bp_summary = get_bioproject_summary(bp_acc, delay_seconds)
            rich_meta = get_rich_sra_metadata(bp_acc, delay_seconds)

            if rich_meta.get("available") and rich_meta.get("merged_csv_data"):
                bp_dir = out_dir / bp_acc
                bp_dir.mkdir(parents=True, exist_ok=True)

                csv_path = bp_dir / f"{bp_acc}_runinfo_merged.csv"
                with csv_path.open("w", newline="", encoding="utf-8") as f:
                    w = csv.DictWriter(f, fieldnames=rich_meta["merged_csv_headers"])
                    w.writeheader()
                    w.writerows(rich_meta["merged_csv_data"])

                srr_list_path = bp_dir / "srr_list.txt"
                srrs = [row["Run"] for row in rich_meta["merged_csv_data"] if "Run" in row]
                srr_list_path.write_text("\n".join(srrs), encoding="utf-8")
            else:
                print(f"\n⚠️  No raw SRA runs found for {bp_acc}. Skipping folder creation and logging to file.")
                with failed_sra_log.open("a", encoding="utf-8") as f:
                    f.write(f"{bp_acc}\n")

            write_markdown_doc(docs_dir, bp_acc, bp_summary, rich_meta)

        log_file = out_dir / "run_metadata.txt"
        write_run_metadata(log_file, log_source_name, len(valid_entries), bool(api_key))

        print("\n\n🎉 Targeted Curation Complete! Output mapped to:", out_dir)
        print(f"📄 Run Log located at: {log_file}")
        return 0

    except KeyboardInterrupt:
        print("\n\n⚠️  Script interrupted. Exiting gracefully.")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())