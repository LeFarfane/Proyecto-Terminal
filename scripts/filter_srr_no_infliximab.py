import csv
from pathlib import Path

input_csv = Path("outputs/GSE270662/GSE270662_runinfo_merged.csv")
output_txt = Path("outputs/GSE270662/srr_list_no_infliximab.txt")

with open(input_csv, newline="") as f:
    reader = csv.DictReader(f)
    runs = [row["Run"] for row in reader if row["treatment"].strip().lower() != "infliximab"]

output_txt.write_text("\n".join(runs) + "\n")
print(f"Wrote {len(runs)} SRR IDs to {output_txt}")
