#!/usr/bin/env bash
set -euo pipefail

# Create a directory for the taxonomy files
mkdir -p resources/taxdump
cd resources/taxdump

# Download tarball with MD5 checksum
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5

# Unpack tarball if md5sum matches
md5sum -c taxdump.tar.gz.md5 &&\
 tar -zxf taxdump.tar.gz ||\
 echo "Wrong checksum, download again!"

# Remove the tarball and checksum
rm taxdump.tar.gz*