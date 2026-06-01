#!/usr/bin/env python3
import pandas as pd
import glob

def main():
    print("=====================================================")
    print("      miRge3.0 Unifier & DESeq2 Filter Pipeline      ")
    print("=====================================================\n")

    print("Scanning directories for miRge3.0 outputs...")
    
    # Use glob to recursively find all specific files in the subdirectories
    cpm_files = glob.glob('**/miR.RPM.csv', recursive=True)
    count_files = glob.glob('**/miR.Counts.csv', recursive=True)

    if not cpm_files or not count_files:
        print("[!] Error: Could not find both miR.RPM.csv and miR.Counts.csv files.")
        print("Make sure you are running this script in the parent directory of your samples.")
        return

    print(f"Found {len(cpm_files)} RPM files and {len(count_files)} Count files.\n")

    THRESHOLD = 100.0

    # ==========================================
    # STEP 1: PROCESS RPM FILES
    # ==========================================
    print(f"--- STEP 1: Processing RPM Files (Threshold: >= {THRESHOLD}) ---")
    rpm_dfs = []
    for f in cpm_files:
        # Setting index_col=0 sets 'miRNA' as the index for perfect concat alignment
        df = pd.read_csv(f, index_col=0)
        rpm_dfs.append(df)
    
    print("Concatenating RPM matrices...")
    
    # OPTION 1 APPLIED: .copy() defragments memory to prevent performance warnings
    merged_rpm = pd.concat(rpm_dfs, axis=1).copy()
    merged_rpm = merged_rpm.fillna(0)
    
    # Calculate row means across all sample columns (Warning resolved!)
    merged_rpm['global_average'] = merged_rpm.mean(axis=1)
    
    # Identify approved miRNAs
    total_mirnas = len(merged_rpm)
    approved_mirnas = merged_rpm[merged_rpm['global_average'] >= THRESHOLD].index.tolist()
    kept_count = len(approved_mirnas)
    
    print(f" -> Total unique miRNAs: {total_mirnas}")
    print(f" -> miRNAs passing >= {THRESHOLD} CPM: {kept_count}")
    print(f" -> miRNAs filtered out: {total_mirnas - kept_count}")

    # Save approved list
    with open('approved_miRNAs_list.txt', 'w') as f:
        for mirna in approved_mirnas:
            f.write(f"{mirna}\n")

    # ==========================================
    # STEP 2: PROCESS RAW COUNTS
    # ==========================================
    print("\n--- STEP 2: Processing Raw Count Files ---")
    count_dfs = []
    for f in count_files:
        df = pd.read_csv(f, index_col=0)
        count_dfs.append(df)
        
    print("Concatenating Raw Count matrices...")
    # Doing the same defragmentation here for consistency and safety
    merged_counts = pd.concat(count_dfs, axis=1).copy()
    
    # Ensure integers for DESeq2 and fill NAs with 0
    merged_counts = merged_counts.fillna(0).astype(int)
    
    # Filter the concatenated counts to only include our approved rows
    filtered_counts = merged_counts.loc[merged_counts.index.isin(approved_mirnas)].copy()
    
    # Guarantee the rows are exactly in the order of our approved list
    filtered_counts = filtered_counts.loc[approved_mirnas]

    # Reset the index so 'miRNA' becomes a normal column again for the CSV export
    filtered_counts.index.name = 'miRNA'
    filtered_counts.reset_index(inplace=True)

    # Export
    output_filename = 'Master_Filtered_Counts_DESeq2.csv'
    filtered_counts.to_csv(output_filename, index=False)

    print(f" -> Final matrix filtered down to {kept_count} rows.")
    print(f"\n[SUCCESS] Master matrix saved as: '{output_filename}'")
    print("=====================================================\n")

if __name__ == "__main__":
    main()