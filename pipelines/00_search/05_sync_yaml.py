#!/usr/bin/env python3
"""
05_sync_yaml.py  —  Sync datasets.yaml from outputs/docs/*.md

Reads every outputs/docs/*.md, parses the structured metadata fields, and
fills any null values in datasets.yaml.  Existing non-null values are NEVER
overwritten — this script only fills gaps.

Matching strategy (in order):
  1. entry id        == doc filename stem  (e.g. GSE84779)
  2. entry accession == doc filename stem
  3. entry bioproject == doc's BioProject(s) field  (catches GSE↔PRJNA pairs)

With --add-new, entries for docs that have no matching yaml entry are appended
as skeleton blocks (status: in_progress) for you to complete manually.

Requires ruamel.yaml for comment-preserving round-trip editing:
    pip install ruamel.yaml   or   conda install -c conda-forge ruamel.yaml

Usage:
    python 05_sync_yaml.py                       # dry run — show changes only
    python 05_sync_yaml.py --confirm             # apply changes
    python 05_sync_yaml.py --add-new --confirm   # apply + append new entries
    python 05_sync_yaml.py --yaml ../../datasets.yaml --docs ../../outputs/docs
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    from ruamel.yaml import YAML
    from ruamel.yaml.scalarstring import DoubleQuotedScalarString as Q
except ImportError:
    print(
        "❌  ruamel.yaml is required for comment-preserving edits.\n"
        "    Install it with:  pip install ruamel.yaml\n"
        "                 or:  conda install -c conda-forge ruamel.yaml"
    )
    sys.exit(1)


# ── Doc field patterns ────────────────────────────────────────────────────────

_PATTERNS = {
    "title":            re.compile(r"\*\*(?:Project )?Title:\*\*\s*(.+)"),
    "pmid":             re.compile(r"\*\*PMID:\*\*\s*(\d+)"),
    "n_runs_total":     re.compile(r"\*\*Total Runs:\*\*\s*(\d+)"),
    "library_strategy": re.compile(r"\*\*Library Strategy:\*\*\s*(.+)"),
    "library_layout":   re.compile(r"\*\*Library Layout:\*\*\s*(.+)"),
    "_bioproject":      re.compile(r"\*\*BioProject\(s\):\*\*\s*(.+)"),
    "_url_geo":         re.compile(r"\*\*GEO URL:\*\*\s*(https://\S+)"),
    "_url_ncbi":        re.compile(r"\*\*NCBI URL:\*\*\s*(https://\S+)"),
}

# Fields that are synced from doc → yaml (in display order)
SYNC_FIELDS = ["title", "pmid", "bioproject", "url",
               "library_strategy", "library_layout", "n_runs_total"]


def parse_doc(md_path: Path) -> dict:
    """Return a dict of parsed metadata from a doc markdown file."""
    text = md_path.read_text(encoding="utf-8")
    out: dict = {}
    for key, pat in _PATTERNS.items():
        m = pat.search(text)
        if not m:
            continue
        val = m.group(1).strip()
        if not val or val in ("N/A", "0"):
            continue
        if key == "pmid":
            out["pmid"] = int(val)
        elif key == "n_runs_total":
            out["n_runs_total"] = int(val)
        elif key == "_bioproject":
            bps = [x.strip() for x in val.split(",") if x.strip().startswith("PRJ")]
            if bps:
                out["bioproject"] = bps[0]
        elif key == "_url_geo":
            out["url"] = val
        elif key == "_url_ncbi":
            out.setdefault("url", val)   # GEO URL takes priority
        else:
            out[key] = val
    return out


# ── Matching ──────────────────────────────────────────────────────────────────

def find_entry(datasets: list, doc_id: str, doc_fields: dict):
    """Return the yaml entry dict that matches this doc, or None."""
    bp = doc_fields.get("bioproject", "")
    for entry in datasets:
        if entry.get("id") == doc_id:
            return entry
        if entry.get("accession") == doc_id:
            return entry
        # bioproject cross-match (e.g. GSE114591.md ↔ PRJNA471862 entry)
        if bp and entry.get("bioproject") == bp:
            return entry
        if bp and entry.get("id") == bp:
            return entry
        if bp and entry.get("accession") == bp:
            return entry
    return None


# ── Updating ──────────────────────────────────────────────────────────────────

def apply_updates(entry: dict, doc_fields: dict) -> list[tuple[str, object, object]]:
    """Fill null/missing fields. Returns list of (field, old_val, new_val)."""
    changes = []
    for field in SYNC_FIELDS:
        if field not in doc_fields:
            continue
        old = entry.get(field)
        if old is None:
            entry[field] = doc_fields[field]
            changes.append((field, None, doc_fields[field]))
    return changes


def make_new_entry(doc_id: str, doc_fields: dict) -> dict:
    """Build a skeleton yaml entry from a doc with no existing yaml match."""
    bp = doc_fields.get("bioproject")
    return {
        "id":                doc_id,
        "accession":         doc_id,
        "bioproject":        bp,
        "pmid":              doc_fields.get("pmid"),
        "status":            "in_progress",
        "output_dir":        f"outputs/{doc_id}",
        "title":             doc_fields.get("title"),
        "url":               doc_fields.get("url"),
        "library_strategy":  doc_fields.get("library_strategy"),
        "library_layout":    doc_fields.get("library_layout"),
        "tissue":            None,
        "diseases":          None,
        "control_label":     None,
        "comparisons":       None,
        "n_runs_total":      doc_fields.get("n_runs_total"),
        "n_samples_analyzed": None,
        "metadata_csv":      None,
        "final_report":      None,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Fill null fields in datasets.yaml from outputs/docs/*.md"
    )
    ap.add_argument("--yaml", default="../../datasets.yaml",
                    help="Path to datasets.yaml")
    ap.add_argument("--docs", default="../../outputs/docs",
                    help="Path to the docs directory")
    ap.add_argument("--confirm", action="store_true",
                    help="Write changes (omit for dry run)")
    ap.add_argument("--add-new", action="store_true",
                    help="Append skeleton entries for docs with no yaml match")
    args = ap.parse_args()

    yaml_path = Path(args.yaml).resolve()
    docs_dir  = Path(args.docs).resolve()

    if not yaml_path.exists():
        print(f"❌  datasets.yaml not found: {yaml_path}")
        return 1
    if not docs_dir.exists():
        print(f"❌  Docs directory not found: {docs_dir}")
        return 1

    mode = "LIVE" if args.confirm else "DRY RUN"
    print(f"\n{'=' * 58}")
    print(f"  datasets.yaml sync  —  {mode}")
    print(f"  yaml : {yaml_path}")
    print(f"  docs : {docs_dir}")
    print(f"{'=' * 58}\n")

    # Load YAML (round-trip mode preserves all comments & formatting)
    ry = YAML()
    ry.preserve_quotes = True
    ry.width = 120
    with yaml_path.open(encoding="utf-8") as fh:
        data = ry.load(fh)

    datasets = data.get("datasets", [])

    # Process every doc file
    total_updates  = 0
    unmatched_docs = []

    for md_path in sorted(docs_dir.glob("*.md")):
        doc_id     = md_path.stem
        doc_fields = parse_doc(md_path)

        if not doc_fields:
            print(f"  ⚠️   {md_path.name:<28} — no parseable fields, skipping")
            continue

        entry = find_entry(datasets, doc_id, doc_fields)

        if entry is None:
            unmatched_docs.append((doc_id, doc_fields))
            print(f"  ❓  {doc_id:<28} — no matching yaml entry")
            continue

        changes = apply_updates(entry, doc_fields)

        if not changes:
            print(f"  ✅  {doc_id:<28} → matched to [{entry['id']}]  (no nulls to fill)")
        else:
            matched_label = entry["id"]
            for field, old, new in changes:
                print(f"  📝  {doc_id:<28} → [{matched_label}].{field}:  null  →  {new!r}")
            total_updates += len(changes)

    # Handle unmatched docs
    newly_added = []
    if unmatched_docs:
        print()
        if args.add_new:
            for doc_id, doc_fields in unmatched_docs:
                new_entry = make_new_entry(doc_id, doc_fields)
                datasets.append(new_entry)
                newly_added.append(doc_id)
                print(f"  ➕  {doc_id:<28} — new skeleton entry appended")
        else:
            print("  ℹ️   Unmatched docs above have no yaml entry.")
            print("      Re-run with --add-new to append skeleton entries.\n")

    # Summary
    print(f"\n  {'─' * 50}")
    print(f"  Fields to update : {total_updates}")
    print(f"  New entries      : {len(newly_added)}")

    if total_updates == 0 and not newly_added:
        print("\n  Nothing to do — datasets.yaml is already in sync.\n")
        return 0

    if not args.confirm:
        print(f"\n  ⚠️   Dry run — no changes written.")
        print(f"      Re-run with --confirm to apply.\n")
        return 0

    # Write back
    with yaml_path.open("w", encoding="utf-8") as fh:
        ry.dump(data, fh)

    print(f"\n  ✅  datasets.yaml updated successfully.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
