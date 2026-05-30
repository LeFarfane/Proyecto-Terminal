#!/usr/bin/env bash
# run_mirge_safe_sequential.sh
# Safely iterates through child directories and processes samples one by one.

set -euo pipefail

THREADS=6

# Function from your script to cleanly extract the sample name
get_run_id() {
  local f="$1"
  local b
  b="$(basename "$f")"
  b="${b%.fastq.gz}"
  b="${b%.fastq}"
  b="${b%.fq.gz}"
  b="${b%.fq}"
  echo "$b"
}

echo "Starting nested sequential processing..."

# OUTER LOOP: Iterate through all child directories
for child_dir in */; do
  # Verify it is a valid directory
  [[ -d "$child_dir" ]] || continue

  # Clean up the directory name by removing the trailing slash
  clean_dir="${child_dir%/}"

  echo "=================================================="
  echo "==> Entering folder: $clean_dir"
  echo "=================================================="

  # Use a subshell to jump into the folder safely
  (
    cd "$child_dir" || exit 1

    # Safely look for .fastq.gz files (prevents errors if none exist)
    shopt -s nullglob
    fastq_files=(*.fastq.gz)
    shopt -u nullglob

    if [[ ${#fastq_files[@]} -eq 0 ]]; then
      echo "    No .fastq.gz files found in $clean_dir. Skipping..."
      exit 0 # Exits this folder, outer loop moves to the next one
    fi

    # INNER LOOP: Process each sample ONE by ONE
    for file in "${fastq_files[@]}"; do
      
      run_id="$(get_run_id "$file")"
      
      # Create a dedicated output folder for THIS single sample
      out_dir="mirge3_output/$run_id"
      mkdir -p "$out_dir"

      echo "    ----------------------------------------------"
      echo "    ==> Analyzing sample: $run_id"
      echo "    File: $file"
      echo "    Out : $clean_dir/$out_dir"
      echo "    ----------------------------------------------"

      # Run miRge3.0 on just this one file, using your exact logic
      stdbuf -oL -eL miRge3.0 \
        -s "$file" \
        -db miRBase \
        -lib /home/genesis/miniconda3/envs/mirge3/miRge3_Lib \
        -on human \
        -ex 0.1 \
        -ie \
        -cpu "$THREADS" \
        -o "$out_dir" \
        -spl \
        -a illumina \
        -nxt 20 \
        -q 20 \
        -NX \
        -m 18 \
        -mEC \
        -gff \
        2>&1 | awk -v rid="$run_id" '{ print strftime("[%Y-%m-%d %H:%M:%S]"), "["rid"]", $0 }' | tee "$out_dir/mirge3.log"
      
      echo "    ==> Finished sample: $run_id"
    done
    
    echo "==> Finished all samples in folder: $clean_dir"
  )
done

echo "Done! All directories and samples processed successfully."