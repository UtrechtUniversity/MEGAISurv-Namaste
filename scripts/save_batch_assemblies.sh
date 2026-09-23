#!/usr/bin/env bash

# Save the currently produced assembly files (FASTA)
# in a separate directory, and combine them in one ZIP
# archive.
# (Works from the main directory, from which the Snakemake workflow
#  is also run.)

batch_dir="batch_results"

mkdir -p ${batch_dir}

assembly_files=$(find results/assembly/ -mindepth 2 -maxdepth 2 -name "assembly.fasta")

echo "Found $(echo ${assembly_files} | wc -w) assembly files"
echo "Copying them to ${batch_dir}..."

for assembly in ${assembly_files}
do
    sample=$(basename $(dirname ${assembly}))
    cp ${assembly} ${batch_dir}/${sample}.fasta
done

echo "Done copying!"
echo "Adding assembly files to ZIP archive"

zip -jm ${batch_dir}/assemblies-batch_X.zip ${batch_dir}/*.fasta

echo "--- Finished! ---"
ls -lh ${batch_dir}/assemblies-batch_X.zip
