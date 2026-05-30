#!/bin/bash
mkdir -p fastqc_out && fastqc -t 7 *.fastq.gz -o fastqc_out && multiqc fastqc_out
