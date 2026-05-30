#!/usr/bin/env bash

# Exit on error, undefined vars, or pipe failures
set -euo pipefail

# Define the tools we need to check
EDIRECT_TOOLS=("esearch" "efetch" "elink" "esummary" "xtract")
SRA_TOOLS=("fasterq-dump")

# We'll use an array to keep track of anything missing
MISSING_DEPS=()

# Function to check if a command exists
check_presence() {
    local cmd=$1
    if ! command -v "$cmd" >/dev/null 2>&1; then
        MISSING_DEPS+=("$cmd")
    fi
}

echo "Starting dependency verification..."

# 1. Check EDirect tools
for tool in "${EDIRECT_TOOLS[@]}"; do
    check_presence "$tool"
done

# 2. Check SRA Toolkit
for tool in "${SRA_TOOLS[@]}"; do
    check_presence "$tool"
done

# Logic for results
if [ ${#MISSING_DEPS[@]} -eq 0 ]; then
    # If the array length is 0, everything was found
    echo "------------------------------------------------"
    echo "✅ Success: Both EDirect and SRA Toolkit are installed."
    echo "------------------------------------------------"
else
    # If there are missing dependencies, log them
    LOG_DIR="tools needed log"
    LOG_FILE="$LOG_DIR/dependency_errors.log"

    # Create the directory (no error if it exists)
    mkdir -p "$LOG_DIR"

    {
        echo "--- Dependency Check Failed: $(date) ---"
        echo "The following tools were not found in the PATH:"
        for missing in "${MISSING_DEPS[@]}"; do
            echo "  - $missing"
        done
        echo "------------------------------------------"
    } >> "$LOG_FILE"

    echo "❌ Missing dependencies detected."
    echo "Details have been written to: $LOG_FILE"
    exit 1
fi