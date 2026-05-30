#!/bin/bash

# 1. Define your adapter sequence (Update this if yours is different!)
# This is the standard Illumina Small RNA adapter
ADAPTER="TGGAATTCTCGGGTGCCAAGG"

# 2. Define the output filename
OUTPUT="mirtrace_config.csv"

# 3. Clear the file if it already exists
> $OUTPUT

# 4. Loop through all fastq.gz files
for file in *.fastq.gz; do
    # Get the absolute path of the file
    FULL_PATH=$(realpath "$file")
    
    # Extract the Sample Name (e.g., SRR14062696 from SRR14062696.fastq.gz)
    # This removes everything after the first dot
    SAMPLE_NAME=$(basename "$file" | cut -d. -f1)
    
    # Write to the CSV: path,name,adapter
    echo "$FULL_PATH,$SAMPLE_NAME,$ADAPTER" >> $OUTPUT
done

echo "Done! ✅ Created $OUTPUT with $(wc -l < $OUTPUT) samples."
