#!/usr/bin/env python3

import pandas as pd

csv_path_a = input("runinfo.csv path: ")
csv_path_b = input("biosample_attricute.csv path: ")

df_a = pd.read_csv(csv_path_a)
df_b = pd.read_csv(csv_path_b)
print("Loaded! :)")

#Checkpoint
print("df_a.shape (rows, cols):", df_a.shape)
print("df_b.shape (rows, cols):", df_b.shape)

if len(df_a) != len(df_b):
    raise ValueError(f"Row mismatch: df_a has {len(df_a)} rows, df_b has {len(df_b)} rows")
else:
    print("Row match: df_a and df_b have same length")



keep = [
    "Run",
    "avgLength",
    "size_MB",
    "LibraryStrategy",
    "LibrarySelection",
    "LibrarySource",
    "LibraryLayout",
    "Platform",
    "Model",
    "BioProject",
    "BioSample",
    "ScientificName",
    "Sex",
    "Disease",
    "Tumor",
    "Histological_Type",
    "Body_Site",
    "Submission"
]
#Replacement part
df_a2 = df_a.loc[:, [c for c in keep if c in df_a.columns]]
print("Shape (rows, cols)", df_a2.shape)

#Concatenation part
assert (df_a["BioSample"].astype(str).values == df_b["biosample_id"].astype(str).values).all(), "Run order doesn't match!"
df_out = pd.concat([df_a2.reset_index(drop=True), df_b.reset_index(drop=True)], axis=1)

#Total of runs
total_runs = len(df_out)

#Get MB
size_mb = pd.to_numeric(df_out["size_MB"], errors="coerce").fillna(0)
total_mb = float(size_mb.sum())
total_gb = total_mb / 1024

print("Shape (rows, cols)", df_out.shape)
print("\nColumns:\n", df_out.columns.tolist())

print(df_out.head())

print(f"Total runs: {total_runs}")
print(f"Total size: {total_mb:,.0f} MB")
print(f"Total size: {total_gb:,.2f} GB")

# Save
df_out.to_csv("runinfo.merged.csv", index=False)
print("OK: wrote runinfo.merged.csv with shape:", df_out.shape)