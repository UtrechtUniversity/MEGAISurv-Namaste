#!/usr/bin/env bash
set -euo pipefail

exec > "${snakemake_log[0]}" 2>&1 # send stdout + stderr to file

json_file="${snakemake_input[0]}"

if [ -z "$json_file" ]
then
    echo "Please provide a (fastplong) JSON file as input argument:"
    echo "bash extract_fastplong_json_summary.sh [example.json]"
    exit 1
elif [ ! -s "$json_file" ]
then
    echo "Input file (${json_file}) does not exist or has size 0"
    exit 1
else
    output_dir="$(dirname ${json_file})/summary"
    mkdir -p ${output_dir}
    output_file="${snakemake_output[0]}"
    (head -n 23 "${json_file}" && echo -e "        }\n}") > ${output_file}
    echo "Extracted just the summary part to ${output_file}"
    exit 0
fi
