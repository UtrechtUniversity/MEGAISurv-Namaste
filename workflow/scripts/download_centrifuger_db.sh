#!/usr/bin/env bash
set -euo pipefail

mkdir -p resources/centrifuger_db
cd resources/centrifuger_db
for index in 1 2 3
do
    output_file="cfr_hpv+gbsarscov2.${index}.cfr"
    if [ -e ${output_file} ]
    then
        echo "File ${output_file} has been downloaded before!"
    else
        wget "https://zenodo.org/records/10023239/files/${output_file}"
    fi
done
