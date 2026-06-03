#!/usr/bin/env python3
import pandas as pd
import glob
import os

def find_files(base_dir, pattern):
    """Searches recursively for a specific file pattern."""
    search_path = os.path.join(base_dir, '**', pattern)
    return glob.glob(search_path, recursive=True)

def main():
    print("=====================================================")
    print("      DESeq2 Auto-Metadata Generator                 ")
    print("=====================================================\n")

    print("This script will automatically find your SRA runinfo file and")
    print("your Master Counts file to generate aligned metadata.")
    
    # 1. Ask for the parent directory (This solves the "two levels deep" problem)
    base_dir = input("\nDrag and drop or paste the Main Dataset Directory here\n(e.g., the a_PRJNA471862 folder) or use . for current directory: ").strip()
    
    if not os.path.isdir(base_dir):
        print("\n[!] Error: That directory does not exist. Exiting.")
        return

    print("\nScanning for required files...")
    
    # 2. Find the files
    runinfo_files = find_files(base_dir, '*runinfo_merged.csv')
    counts_files = find_files(base_dir, 'Master_Filtered_Counts_DESeq2.csv')

    # Safety checks
    if not runinfo_files:
        print("[!] Error: Could not find any file ending in 'runinfo_merged.csv'.")
        return
    if not counts_files:
        print("[!] Error: Could not find 'Master_Filtered_Counts_DESeq2.csv'.")
        return

    # If multiple files are found, default to the first one (or you could add a menu here)
    runinfo_target = runinfo_files[0]
    counts_target = counts_files[0]

    print(f" [+] Found RunInfo: {os.path.basename(runinfo_target)}")
    print(f" [+] Found Counts : {os.path.basename(counts_target)}")

    # 3. Process the data
    try:
        # Load count headers
        counts_df = pd.read_csv(counts_target, nrows=0) 
        sample_ids = counts_df.columns[1:].tolist()

        # Load runinfo
        runinfo_df = pd.read_csv(runinfo_target)

        # Filter
        meta_filtered = runinfo_df[runinfo_df['Run'].isin(sample_ids)].copy()

        # Detect whichever diagnostic column is present and has data
        DIAG_CANDIDATES = ['disease_status', 'disease_state', 'disease', 'Disease',
                           'Diagnosis', 'diagnosis', 'disease state', 'isolate', 'subject']
        col_diag = next(
            (c for c in DIAG_CANDIDATES
             if c in runinfo_df.columns and not runinfo_df[c].isna().all()),
            None
        )
        if col_diag is None:
            print("[!] Warning: no diagnostic column found. Available columns:")
            print(f"    {list(runinfo_df.columns)}")

        # Select columns flexibly (it will keep them if they exist)
        cols_to_keep = ['Run'] + ([col_diag] if col_diag else []) + ['group', 'inflammation', 'tissue', 'sex', 'age', 'source_name']
        meta_filtered = meta_filtered[[c for c in cols_to_keep if c in meta_filtered.columns]]

        # Rename for DESeq2 — always produce a 'disease' column regardless of source name
        rename_map = {'Run': 'sample_id'}
        if col_diag and col_diag != 'disease':
            rename_map[col_diag] = 'disease'
        meta_filtered = meta_filtered.rename(columns=rename_map)

        if col_diag:
            print(f" [+] Diagnostic column: '{col_diag}' -> 'disease'")

        # CRITICAL: Enforce column order
        meta_filtered['sample_id'] = pd.Categorical(meta_filtered['sample_id'], categories=sample_ids, ordered=True)
        meta_filtered = meta_filtered.sort_values('sample_id')

        # 4. Save the output next to the Master Counts file
        # This is the smartest place to save it, so R can easily find both files together
        output_dir = os.path.dirname(counts_target)
        output_path = os.path.join(output_dir, "sample_metadata.csv")
        
        meta_filtered.to_csv(output_path, index=False)

        print("\n=====================================================")
        print(f"Matched {len(sample_ids)} columns to {len(meta_filtered)} metadata rows.")
        if 'disease' in meta_filtered.columns:
            print("Disease Breakdown:")
            for disease, count in meta_filtered['disease'].value_counts().items():
                print(f"  - {disease}: {count}")
        else:
            print("[!] No disease column in output (no diagnostic column matched).")
        print(f"\n[SUCCESS] Saved perfectly aligned metadata to:")
        print(f" -> {output_path}")
        print("=====================================================\n")

    except Exception as e:
        print(f"\n[!] An error occurred during processing: {e}")

if __name__ == "__main__":
    main()