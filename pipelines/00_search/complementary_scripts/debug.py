#!/usr/bin/env python3
import subprocess
import csv
import io


def run_cmd(cmd, input_text=None):
    print(f"\n🚀 Running: {' '.join(cmd)}")
    try:
        p = subprocess.run(
            cmd,
            input=input_text.encode('utf-8') if input_text else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False  # We want to capture errors, not crash
        )
        if p.returncode != 0:
            print(f"❌ STDERR:\n{p.stderr.decode('utf-8', errors='replace')}")
        return p.stdout.decode('utf-8', errors='replace')
    except Exception as e:
        print(f"💥 EXCEPTION: {e}")
        return ""


def debug_pipeline():
    gse = "GSE114591"
    print(f"=== Debugging Pivot for {gse} ===")

    # STEP 1: elink (GDS -> BioProject)
    elink_out = run_cmd(["elink", "-db", "gds", "-id", gse, "-target", "bioproject"])
    print(f"📥 elink Output Preview (First 300 chars):\n{elink_out[:300]}")

    if not elink_out.strip():
        print("🛑 elink produced no output. Stopping.")
        return

    # STEP 2: efetch docsum
    docsum_out = run_cmd(["efetch", "-format", "docsum"], input_text=elink_out)
    print(f"📥 docsum Output Preview (First 300 chars):\n{docsum_out[:300]}")

    # STEP 3: xtract BioProject
    xtract_out = run_cmd(["xtract", "-pattern", "DocumentSummary", "-element", "Project_Acc"], input_text=docsum_out)
    print(f"📥 xtract Output:\n{xtract_out}")

    bp_acc = xtract_out.strip().split('\n')[0].split('\t')[0] if xtract_out else ""
    print(f"🎯 Extracted BioProject: '{bp_acc}'")

    if not bp_acc.startswith("PRJ"):
        print("🛑 Failed to get a valid PRJ accession. Stopping.")
        return

    # STEP 4: esearch BioProject
    esearch_out = run_cmd(["esearch", "-db", "bioproject", "-query", bp_acc])
    print(f"📥 esearch Output Preview:\n{esearch_out[:300]}")

    # STEP 5: elink to SRA
    elink_sra_out = run_cmd(["elink", "-target", "sra"], input_text=esearch_out)
    print(f"📥 elink SRA Output Preview:\n{elink_sra_out[:300]}")

    # STEP 6: efetch RunInfo
    runinfo_out = run_cmd(["efetch", "-format", "runinfo"], input_text=elink_sra_out)
    print(f"📥 RunInfo Preview (First 500 chars):\n{runinfo_out[:500]}")

    srr_count = sum(1 for line in runinfo_out.splitlines() if line.startswith("SRR"))
    print(f"\n✅ Total SRR runs found: {srr_count}")

    if srr_count == 0:
        print("🛑 No SRR runs found to extract BioSample IDs from. Stopping.")
        return

    # STEP 7: Extract BioSample attributes using the user's xtract logic
    print("\n" + "=" * 40)
    print("=== STEP 7: BioSample Attribute Extraction ===")

    reader = csv.DictReader(io.StringIO(runinfo_out))
    first_samn = None
    for row in reader:
        if "BioSample" in row and row["BioSample"].startswith("SAMN"):
            first_samn = row["BioSample"]
            break

    if not first_samn:
        print("🛑 Could not find a SAMN BioSample ID in the RunInfo data. Stopping.")
        return

    print(f"🎯 Targeting First BioSample: '{first_samn}'")

    fetch_xml = run_cmd(["efetch", "-db", "biosample", "-id", first_samn, "-format", "xml"])
    if not fetch_xml.strip():
        print("🛑 Failed to fetch XML for BioSample.")
        return

    attributes = run_cmd(["xtract", "-pattern", "Attribute", "-element", "@attribute_name", "-element", "."],
                         input_text=fetch_xml)
    print(f"\n📥 BioSample Attributes for {first_samn}:\n{attributes}")


if __name__ == "__main__":
    debug_pipeline()