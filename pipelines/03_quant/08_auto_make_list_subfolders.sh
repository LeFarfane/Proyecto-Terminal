#!/bin/bash

# Auto-navigate to sra_fastq/muestras_clasificadas from the dataset root
TARGET_DIR="sra_fastq/muestras_clasificadas"
if [[ ! -d "$TARGET_DIR" ]]; then
    echo "❌ No se encontró $TARGET_DIR en: $(pwd)"
    exit 1
fi
cd "$TARGET_DIR"
echo "📂 Trabajando en: $(pwd)"

# Loop through all child folders in the current directory
for dir in */; do
    # Verify it is actually a directory
    if [ -d "$dir" ]; then
        echo "Processing folder: $dir"
        
        # Use a subshell (parentheses) to change directories safely
        (
            cd "$dir" || exit
            
            # Run your exact command
            # The 2>/dev/null hides errors if a folder happens to have no .fastq.gz files
            ls -1 -- *.fastq.gz > fasta_list.txt 2>/dev/null
        )
    fi
done

echo "Done! All child folders have been processed."
