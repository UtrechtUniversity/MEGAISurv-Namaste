#!/usr/bin/env bash

# Save the currently produced result files (CSV, TSV)
# in a separate directory.
# (Works from the main directory, from which the Snakemake workflow
#  is also run.)

batch_dir="batch_results"

mkdir -p ${batch_dir}/batch_X

results=$(find results/ -mindepth 1 -maxdepth 1 -regex ".*\\.[ct]sv[\\.gz]*")

echo "Found $(echo ${results} | wc -w) results files"
echo "Copying them to ${batch_dir}..."

for result in ${results}
do
    cp ${result} ${batch_dir}/batch_X/
done

echo "--- Finished! ---"
ls -lh ${batch_dir}/batch_X/

echo
echo "!! Remember to rename the batch before running this script again !!"
