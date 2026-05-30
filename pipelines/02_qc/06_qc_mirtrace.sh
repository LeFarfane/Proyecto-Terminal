#!/bin/bash
# Using $1 allows you to pass any config file to the script
mirtrace qc --species hsa -c "$1" --output-dir mirtrace_resultsmirtrace --species hsa --config mirtrace_config.csv --output-dir mirtrace_results
