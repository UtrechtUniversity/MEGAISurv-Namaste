#!/usr/bin/env bash

# Merge read files from the same samples by concatenating the FASTQ files.
# (Copy and paste the file names and accession IDs as these are few.)

INPUT_DIR="data/raw/"
OUTPUT_DIR="data/deduplicated/"

mkdir -p ${OUTPUT_DIR}

merge_files () {
    if [ ! -s ${OUTPUT_DIR}$3 ]
    then
        echo "Merging $1 and $2 into $3 !"
        cat ${INPUT_DIR}$1 ${INPUT_DIR}$2 > ${OUTPUT_DIR}$3
    else
        echo "Files merged previously:"
        ls -lh ${OUTPUT_DIR}$3
    fi
}

merge_files ERR11557117.fastq.gz ERR11581405.fastq.gz SAMEA112228225.fastq.gz

merge_files ERR11041251.fastq.gz ERR11041272.fastq.gz SAMEA112797804.fastq.gz

merge_files ERR11041253.fastq.gz ERR11041274.fastq.gz SAMEA112797835.fastq.gz

merge_files ERR11041254.fastq.gz ERR11041275.fastq.gz SAMEA112797836.fastq.gz

merge_files ERR11041256.fastq.gz ERR11041277.fastq.gz SAMEA112797838.fastq.gz

merge_files ERR11041257.fastq.gz ERR11041278.fastq.gz SAMEA112797839.fastq.gz

merge_files ERR11041264.fastq.gz ERR11041285.fastq.gz SAMEA112797888.fastq.gz

merge_files ERR11041289.fastq.gz ERR11041268.fastq.gz SAMEA112797924.fastq.gz

merge_files ERR5261868.fastq.gz ERR7927172.fastq.gz SAMEA7999667.fastq.gz

merge_files SRR10394893_1.fastq.gz SRR10394891_1.fastq.gz SAMN13195585.fastq.gz

# To make the analysis a bit easier, link the other read files to the same
# directory.
cd ${OUTPUT_DIR}
ln -s ../../${INPUT_DIR}*.fastq.gz .
# And remember to remove the ones that have been deduplicated...
for run in ERR11557117 ERR11581405 ERR11041251 ERR11041272 ERR11041253 ERR11041274 ERR11041254 ERR11041275 ERR11041256 ERR11041277 ERR11041257 ERR11041278 ERR11041264 ERR11041285 ERR11041289 ERR11041268 ERR5261868 ERR7927172 SRR10394893_1 SRR10394891_1
do
    unlink ${run}*
done
