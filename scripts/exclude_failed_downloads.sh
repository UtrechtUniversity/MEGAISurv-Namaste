#!/usr/bin/env bash

# Run this after running the 'auto-download' variant:
# it finds sample accessions that could not be downloaded,
# and removes these from the input files so that the
# snakemake workflow can continue with the working samples.

# Grab the input file from snakemake's parameters file
batch_file="$(grep "input" config/parameters.yaml | sed -e 's/input: "//g' -e 's/"//g')"

# Save a list of accession IDs from failed downloads
failed_download="results/download_failed.txt"

# Find accession IDs of samples that did not download from their log files
grep "sracha cannot download" log/download_raw_reads/*txt |\
 grep -o -e "SAM[A-Z0-9]*" > ${failed_download}

# Remove the accessions for which download failed from the snakemake input file
mv ${batch_file} "${batch_file}_original"
grep -v -f ${failed_download} "${batch_file}_original" > ${batch_file}
