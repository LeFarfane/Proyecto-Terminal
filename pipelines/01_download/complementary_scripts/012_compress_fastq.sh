#!/bin/bash

# This script compresses all .fastq files in the current directory
# using pigz (Parallel Implementation of GZip) with 6 threads.

echo "Starting compression of .fastq files using 6 threads..."

# Execute the pigz command
pigz -p 6 *.fastq

echo "Compression complete!"