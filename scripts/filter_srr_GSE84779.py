import csv
from pathlib import Path

exclude = {"CP-P", "NP-P", "CLDP"}

input_csv = Path("outputs/a_GSE84779/GSE84779_runinfo_merged.csv")
output_txt = Path("outputs/a_GSE84779/srr_list_filtered.txt")

with open(input_csv, newline="") as f:
    reader = csv.DictReader(f)
    runs = [row["Run"] for row in reader if row["disease_state"].strip() not in exclude]

output_txt.write_text("\n".join(runs) + "\n")
print(f"Wrote {len(runs)} SRR IDs to {output_txt}")
