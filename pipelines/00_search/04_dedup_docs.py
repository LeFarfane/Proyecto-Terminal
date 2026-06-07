#!/usr/bin/env python3
"""
04_dedup_docs.py  —  Deduplicate outputs/docs/

For every GSE*.md file, reads the "BioProject(s):" line and checks whether a
PRJNA*.md exists for that same BioProject.  If it does, the PRJNA file is
redundant (the GSE file already contains the GEO title, PMID, and abstract)
and is deleted.

Runs in dry-run mode by default.  Pass --confirm to actually delete.

Usage:
    python 04_dedup_docs.py
    python 04_dedup_docs.py --confirm
    python 04_dedup_docs.py --docs ../../outputs/docs
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


BIOPROJECT_RE = re.compile(r"\*\*BioProject\(s\):\*\*\s*(.+)")


def extract_bioprojects(md_path: Path) -> list[str]:
    """Return all PRJNA* IDs listed on the BioProject(s) line of a doc file."""
    for line in md_path.read_text(encoding="utf-8").splitlines():
        m = BIOPROJECT_RE.search(line)
        if m:
            raw = m.group(1).strip()
            return [tok.strip() for tok in raw.split(",") if tok.strip().startswith("PRJ")]
    return []


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Remove PRJNA doc files that are already covered by a GSE doc file."
    )
    ap.add_argument(
        "--docs", default="../../outputs/docs",
        help="Path to the docs directory (default: ../../outputs/docs)",
    )
    ap.add_argument(
        "--confirm", action="store_true",
        help="Actually delete files (omit for a safe dry run)",
    )
    args = ap.parse_args()

    docs_dir = Path(args.docs).resolve()
    if not docs_dir.exists():
        print(f"❌ Docs directory not found: {docs_dir}")
        return 1

    mode = "LIVE" if args.confirm else "DRY RUN"
    print(f"\n{'=' * 55}")
    print(f"  Doc deduplicator  —  {mode}")
    print(f"  Docs dir: {docs_dir}")
    print(f"{'=' * 55}\n")

    # Build a map: PRJNA_ID -> GSE file that claims it
    covered: dict[str, Path] = {}
    for gse_file in sorted(docs_dir.glob("GSE*.md")):
        for bp_id in extract_bioprojects(gse_file):
            covered[bp_id] = gse_file

    if not covered:
        print("ℹ️  No BioProject references found in any GSE doc file. Nothing to do.")
        return 0

    # Find PRJNA doc files whose ID appears in the covered map
    to_delete: list[tuple[Path, Path]] = []   # (prjna_file, gse_file_that_covers_it)
    for prjna_file in sorted(docs_dir.glob("PRJNA*.md")):
        bp_id = prjna_file.stem          # e.g. "PRJNA331127"
        if bp_id in covered:
            to_delete.append((prjna_file, covered[bp_id]))

    if not to_delete:
        print("✅ No duplicate PRJNA doc files found. Everything is clean.")
        return 0

    print(f"Found {len(to_delete)} redundant PRJNA file(s):\n")
    for prjna_file, gse_file in to_delete:
        print(f"  🗑️  {prjna_file.name:<25}  ←  covered by  {gse_file.name}")

    if not args.confirm:
        print(f"\n⚠️  Dry run — no files deleted.")
        print(f"    Re-run with --confirm to delete.\n")
        return 0

    print()
    deleted = 0
    for prjna_file, _ in to_delete:
        try:
            prjna_file.unlink()
            print(f"  ✅ Deleted: {prjna_file.name}")
            deleted += 1
        except OSError as e:
            print(f"  ❌ Could not delete {prjna_file.name}: {e}")

    print(f"\n🎉 Done — {deleted}/{len(to_delete)} file(s) removed.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
