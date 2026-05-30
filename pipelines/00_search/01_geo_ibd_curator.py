#!/usr/bin/env python3
"""
geo_ibd_curator.py (The Fully Automated True Pivot Edition)

Description:
    This script automates the discovery, filtering, and metadata extraction of
    high-throughput sequencing datasets from the NCBI GEO and SRA databases.
    It features aggressive terminal logging for pipeline transparency.

Key Features:
    1. Always-On True Pivot: Uses direct text-search on the BioProject database unconditionally.
    2. Centralized Outputs: All data, ignore lists, and error logs are stored in an 'outputs' folder.
    3. Blacklist Memory: Uses 'ignore_list.txt' to permanently bypass rejected datasets.
    4. Target Output: Generates SRR target lists for downstream FASTQ downloading.
"""

import argparse
import csv
import io
import os
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Set

DEFAULT_QUERY = (
    '(("inflammatory bowel diseases"[MeSH Terms] '
    'OR "crohn disease"[MeSH Terms] '
    'OR "colitis, ulcerative"[MeSH Terms] '
    'OR "inflammatory bowel disease"[All Fields]) '
    'AND ("micrornas"[MeSH Terms] OR "micrornas"[All Fields] OR mirna[All Fields])) '
    'AND ("Homo sapiens"[Organism] '
    'AND "Non-coding RNA profiling by high throughput sequencing"[Filter] '
    'AND ("50"[n_samples] : "10000"[n_samples]))'
)

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


def load_ignore_list(filepath: Path) -> Set[str]:
    """Reads a text file containing GSE IDs that should be permanently skipped."""
    ignore_set = set()
    if filepath.exists():
        with filepath.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip().upper()
                if line and not line.startswith("#"):
                    ignore_set.add(line)
    return ignore_set


def gds_docsum_from_query(query: str, delay: float) -> str:
    """Executes the initial search against GEO and fetches the XML document summaries."""
    print("\n🔍 Executing GEO query...")
    search_out = run_network_cmd(["esearch", "-db", "gds", "-query", query], delay)
    print("📥 Fetching document summaries...")
    return run_network_cmd(["efetch", "-db", "gds", "-format", "docsum"], delay, input_text=search_out)


def parse_gds_docsum(docsum_xml: str) -> List[Dict[str, str]]:
    """Parses the GEO XML natively to prevent nested GSM sample IDs from shifting columns."""
    print("\n⚙️  Parsing XML payload natively...")
    entries = []
    try:
        root = ET.fromstring(docsum_xml)
        for docsum in root.findall('.//DocumentSummary'):
            gds_uid = docsum.findtext('Id', default="").strip()
            acc = docsum.findtext('Accession', default="").strip()
            title = docsum.findtext('title', default="").strip()
            summary = docsum.findtext('summary', default="").strip()

            pmids = ""
            pmid_node = docsum.find('.//PubMedIds/int')
            if pmid_node is not None and pmid_node.text:
                pmids = pmid_node.text.strip()

            if re.match(r"^GSE\d+$", acc):
                entries.append({
                    "gds_uid": gds_uid,
                    "GSE_ID": acc,
                    "Title": title,
                    "Summary": summary,
                    "PMID": pmids,
                })
    except ET.ParseError as e:
        print(f"⚠️ XML Parsing Error: {e}")
    return entries


def get_rich_sra_metadata(gse_id: str, delay: float) -> Dict[str, Any]:
    """
    The True Pivot: Bypasses broken NCBI elink tables by text-searching the
    BioProject database directly for the GSE ID.
    """
    meta = {
        "available": False, "srr_count": 0, "total_mb": 0.0, "total_gb": 0.0,
        "bioprojects": set(), "strategies": set(), "layouts": set(), "platforms": set(),
        "diseases": set(), "tissues": set(),
        "merged_csv_data": [], "merged_csv_headers": []
    }

    try:
        # STEP 1: GSE -> BioProject (Direct Search Method)
        print("\n--- STEP 1: GSE -> BioProject ---")
        search_bp = run_network_cmd(["esearch", "-db", "bioproject", "-query", gse_id], delay)
        docsum_bp = run_network_cmd(["efetch", "-format", "docsum"], delay, input_text=search_bp)
        xtract_bp = run_local_cmd(["xtract", "-pattern", "DocumentSummary", "-element", "Project_Acc"],
                                  input_text=docsum_bp)

        bp_accs = [acc.strip() for acc in xtract_bp.splitlines() if acc.strip().startswith("PRJ")]

        # Fallback: If esearch fails, try the old elink method just in case
        if not bp_accs:
            print(f"\n⚠️  Direct search failed for {gse_id}. Attempting fallback elink...")
            link_bp = run_network_cmd(["elink", "-db", "gds", "-id", gse_id, "-target", "bioproject"], delay)
            docsum_bp = run_network_cmd(["efetch", "-format", "docsum"], delay, input_text=link_bp)
            xtract_bp = run_local_cmd(["xtract", "-pattern", "DocumentSummary", "-element", "Project_Acc"],
                                      input_text=docsum_bp)
            bp_accs = [acc.strip() for acc in xtract_bp.splitlines() if acc.strip().startswith("PRJ")]

        print(f"🎯 [RESULT]: Extracted BioProjects: {bp_accs}")

        if not bp_accs:
            print(f"🛑 [STOP]: No BioProjects found for {gse_id}.")
            return meta

        meta["bioprojects"].update(bp_accs)

        # STEP 2: BioProject -> SRA RunInfo
        print("\n--- STEP 2: BioProject -> SRA RunInfo ---")
        runs = []
        biosample_ids = set()

        for bp_acc in bp_accs:
            search_bp_sra = run_network_cmd(["esearch", "-db", "bioproject", "-query", bp_acc], delay)
            link_sra = run_network_cmd(["elink", "-target", "sra"], delay, input_text=search_bp_sra)
            runinfo_raw = run_network_cmd(["efetch", "-format", "runinfo"], delay, input_text=link_sra)

            if not runinfo_raw.strip():
                continue

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
            print(f"🛑 [STOP]: No SRA runs found for BioProjects: {bp_accs}")
            return meta

        meta["available"] = True
        meta["srr_count"] = len(runs)
        meta["total_gb"] = meta["total_mb"] / 1024

        # STEP 3: BioSample Attributes (Batched for speed, parsed natively)
        print("\n--- STEP 3: BioSample Attributes ---")
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

        # STEP 4: Final Merge
        print("\n--- STEP 4: Final Metadata Merge ---")
        sorted_attr_keys = sorted(list(all_attr_keys))
        meta["merged_csv_headers"] = [k for k in KEEP_COLS] + sorted_attr_keys

        for run in runs:
            bs_id = run.get("BioSample", "")
            bs_data = biosample_attrs.get(bs_id, {})
            merged_row = run.copy()
            for key in sorted_attr_keys:
                merged_row[key] = bs_data.get(key, "")
            meta["merged_csv_data"].append(merged_row)

        for key in ["bioprojects", "strategies", "layouts", "platforms", "diseases", "tissues"]:
            meta[key] = {item for item in meta[key] if item}

        print("✅ Merge successful.")
        return meta

    except Exception as e:
        print(f"\n⚠️ Pipeline Error for {gse_id}: {e}")
        return meta


def normalize_text(s: str) -> str:
    """Removes excessive newlines and carriage returns from text blocks."""
    if not s:
        return ""
    return re.sub(r"\n{3,}", "\n\n", s.replace("\r\n", "\n").replace("\r", "\n")).strip()


def write_markdown_doc(docs_dir: Path, entry: Dict[str, Any]) -> None:
    """Generates a rich Markdown file for a single GSE dataset."""
    gse = entry["GSE_ID"]
    meta = entry.get("rich_meta", {})

    bioprojects = ", ".join(meta.get("bioprojects", [])) or "N/A"
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
    (docs_dir / f"{gse}.md").write_text(md, encoding="utf-8")


def write_run_metadata(log_path: Path, query: str, args: argparse.Namespace, total_processed: int, has_key: bool):
    """Generates a master text file recording the parameters and results of the pipeline execution."""
    log = f"""---
GEO IBD Curator Execution Log
Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
---

## Execution Parameters
- Query: {query}
- Max Requested: {args.max}
- Output Directory: {args.out}
- API Rate Limit: {'10 requests/sec (Key Active)' if has_key else '3 requests/sec (No Key)'}

## Results
- Total Datasets Processed: {total_processed}
"""
    log_path.write_text(log, encoding="utf-8")


def main() -> int:
    try:
        ap = argparse.ArgumentParser()
        ap.add_argument("--max", type=int, default=30, help="Max GSEs to keep after parsing")
        # CHANGED: Default is now two directories up -> /Pipelines/outputs
        ap.add_argument("--out", default="../../outputs", help="Master output folder")
        ap.add_argument("--ignore-list", default="ignore_list.txt",
                        help="Filename for the ignore list (saved inside output folder)")
        args = ap.parse_args()

        print("\n" + "=" * 50)
        print("👋 Welcome to the GEO Data Curator Pipeline (Fully Automated Edition)")
        print("=" * 50)

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

        # File Setup (now centralized in outputs/)
        ignore_filepath = out_dir / args.ignore_list
        failed_sra_log = out_dir / "failed_sra_datasets.txt"

        if not ignore_filepath.exists():
            ignore_filepath.write_text("# Add GSE IDs here to skip them in future runs (one per line)\n",
                                       encoding="utf-8")

        # Initialize the failed SRA log clean for each run
        failed_sra_log.write_text(f"# Failed SRA Links from run on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
                                  encoding="utf-8")

        ignore_set = load_ignore_list(ignore_filepath)
        if ignore_set:
            print(f"\n🛡️  Loaded {len(ignore_set)} datasets to ignore from {ignore_filepath.name}")

        print(f"\n📝 Default Query:\n{DEFAULT_QUERY}\n")
        user_query = input("Enter custom GEO query (or press Enter for default):\n> ").strip()
        final_query = user_query if user_query else DEFAULT_QUERY

        docsum = gds_docsum_from_query(final_query, delay_seconds)
        all_entries = parse_gds_docsum(docsum)

        if not all_entries:
            print("\n❌ No results returned.")
            return 0

        valid_entries = []
        ignored_count = 0
        for e in all_entries:
            if e["GSE_ID"].upper() in ignore_set:
                ignored_count += 1
                continue
            valid_entries.append(e)

        if ignored_count > 0:
            print(f"\n🚫 Bypassed {ignored_count} datasets matched against your ignore list.")

        entries_to_process = valid_entries[:args.max]

        if not entries_to_process:
            print("\n✅ All returned datasets were successfully bypassed by your ignore list.")
            return 0

        print(f"\n📊 Processing top {len(entries_to_process)} datasets...")

        for i, e in enumerate(entries_to_process, 1):

            print("\n" + "=" * 60)
            print(f"🧬 PROCESSING DATASET {i}/{len(entries_to_process)}: {e['GSE_ID']}")
            print("=" * 60)

            # Core extraction always runs now
            rich_meta = get_rich_sra_metadata(e["GSE_ID"], delay_seconds)
            e["rich_meta"] = rich_meta

            if rich_meta.get("available") and rich_meta.get("merged_csv_data"):
                gse_dir = out_dir / e["GSE_ID"]
                gse_dir.mkdir(parents=True, exist_ok=True)

                csv_path = gse_dir / f"{e['GSE_ID']}_runinfo_merged.csv"
                with csv_path.open("w", newline="", encoding="utf-8") as f:
                    w = csv.DictWriter(f, fieldnames=rich_meta["merged_csv_headers"])
                    w.writeheader()
                    w.writerows(rich_meta["merged_csv_data"])

                srr_list_path = gse_dir / "srr_list.txt"
                srrs = [row["Run"] for row in rich_meta["merged_csv_data"] if "Run" in row]
                srr_list_path.write_text("\n".join(srrs), encoding="utf-8")
            else:
                print(f"\n⚠️  No raw SRA runs found for {e['GSE_ID']}. Skipping folder creation and logging to file.")
                with failed_sra_log.open("a", encoding="utf-8") as f:
                    f.write(f"{e['GSE_ID']}\n")

            write_markdown_doc(docs_dir, e)

        log_file = out_dir / "run_metadata.txt"
        write_run_metadata(log_file, final_query, args, len(entries_to_process), bool(api_key))

        print("\n\n🎉 Curation Complete! Output mapped to:", out_dir)
        print(f"📄 Run Log located at: {log_file}")
        print(f"📄 Failed SRA records logged at: {failed_sra_log}")
        print(f"📄 Ignore list located at: {ignore_filepath}")
        return 0

    except KeyboardInterrupt:
        print("\n\n⚠️  Script interrupted. Exiting gracefully.")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())