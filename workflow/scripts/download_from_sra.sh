#!/usr/bin/env bash
set -uo pipefail
set +e # sracha cannot download all accessions, for various reasons like lacking quality values

output_dir="resources/public_metagenomes/"
threads=4

print_help() {
  # Display help for this script
  echo "Download metagenomes from accession IDs in a text file"
  echo
  echo "Syntax: download_metagenome_batches.sh -s [sample_accession] -d [directory] [-h]"
  echo "Options:"
  echo "-s/--sample    Sample accession ID to download (from SRA/ENA)"
  echo "-d/--directory Directory in which to download the files"
  echo "               (default=resources/public_metagenomes/)"
  echo "-t/--threads   Number of CPU threads to use (default=4)"
  echo "-h/--help      Print this help message"
  echo
  exit 0
}

message() {
    # Write a message with timestamp
    echo -e "\e[34m[$(date +"%Y-%m-%d %H:%M:%S")]\e[0m $1"
}

while [[ $# -gt 0 ]]
do
  case "$1" in
    -s|--sample ) sample="$2"
    shift
    shift
    ;;
    -d|--directory ) output_dir="$2"
    shift
    shift
    ;;
    -t|--threads ) threads="$2"
    shift
    shift
    ;;
    -h|--help ) print_help
    ;;
    * )
    echo "Unknown option: $1"
    print_help
  esac
done

sample_file="${output_dir}/${sample}.fastq.gz"

if [ -s ${sample_file} ]
# If the output file already exists
then
    message "Reads for sample ${sample} already downloaded!"
    ls -lh ${sample_file}
else
# If not, start downloading
# Work around bad samples by using an 'download, or do nothing' construction
    message "Looking up run accessions for sample: ${sample}"

    run_accessions=( $(curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${sample}&result=read_run&fields=run_accession" | grep -v "run_accession") )

    message "Found ${#run_accessions[@]} runs: $(echo ${run_accessions[@]})"

    message "Preparing to download metagenomic reads,"
    message "and store them in the directory: ${output_dir}."
    echo "----------"

    mkdir -p ${output_dir}

    {
        sracha get -t 4 --yes --no-progress --output-dir ${output_dir} ${run_accessions[@]} && {
            if [ ${#run_accessions[@]} > 1 ]
            then
                # If there are multiple runs per sample, concatenate to one file
                > ${sample_file}
                for run in ${run_accessions[@]}
                do
                    run_file=${output_dir}/${run}*fastq.gz
                    cat ${run_file} >> ${sample_file}
                    # and remove the separate runs
                    rm ${run_file}
                done
            else
                # Otherwise, just rename the one file
                mv ${output_dir}/${run_accessions[0]}*.fastq.gz "${sample_file}"
            fi
            message "Done! Moved output to ${sample_file}."
        }
        } || {
        message "sracha cannot download ${sample}"
        message "Skipping this sample and moving to the next..."
    }
fi

message "----- Done! -----"
