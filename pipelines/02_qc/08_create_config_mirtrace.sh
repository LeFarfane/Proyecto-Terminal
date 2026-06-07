#!/bin/bash

# 1. Pick the 3' adapter from the dataset's library-prep kit (see ../kits.tsv).
#    Set the kit explicitly with KIT=..., or let it resolve from DATASET=<id>.
#    Examples:
#      KIT=nextflex_v3 ./05_create_config_mirtrace.sh
#      DATASET=GSE272890 ./05_create_config_mirtrace.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/kit_params.sh
source "$SCRIPT_DIR/../lib/kit_params.sh"

KIT="$(resolve_kit)"
load_kit_params "$KIT"
ADAPTER="$MIRTRACE_ADAPTER"
echo "Kit: $KIT  ->  miRTrace adapter: $ADAPTER"

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
