# Extra functions included in the workflow

## Convenience functions

Namaste includes a number of functions that do not generate output by
themselves, but are required to make other analysis steps work.
By including them in the automated analysis, the user can run the
whole workflow in one simple command and does not have to prepare
things manually, such as downloading reference databases.

### Downloading databases

Namaste includes scripts for downloading reference databases for:

- Antibiotic resistance gene screening
[ResFinder](https://bitbucket.org/genomicepidemiology/resfinder_db/)
- Antibiotic resistance mutation screening
[AMRFinder](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest),
    - included in [MetaPointFinder](https://github.com/aldertzomer/metapointfinder)
- Taxonomic classification
[Centrifuger's cfr_hpv+gbsarscov2](https://zenodo.org/records/10023239)
- Translation to scientific names
[NCBI taxdump](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
- Plasmid/virus/chromosome predictions [geNomad](https://zenodo.org/records/14886553)

### Summarising results

Namaste has a couple of custom R scripts embedded in the automated workflow
to parse the output from different tools and write them to easy-to-use
tabular formats (TSV).

### Removing unassembled samples

Namaste relies on [Snakemake](https://snakemake.readthedocs.io/en/stable/)
to run the whole workflow. To generate the final output, it requires all
intermediate steps for all samples to pass. However, some samples may not
yield proper assemblies (for example, samples that have very few reads
or negative control samples). To successfully finish the workflow, these
samples may be moved 'out of the way' with the included Python script
[`exclude_failed_assemblies.py`](https://github.com/UtrechtUniversity/MEGAISurv-Namaste/blob/main/scripts/exclude_failed_assemblies.py).

## Bonus features

The basis of Namaste is to assemble contigs, detect antibiotic resistance
determinants (genes/mutations), and assign taxonomic (species) classifications.

As a little bonus, statistics are collected along the way to determine the
amount of sequence information contained within each step: from raw reads
to high-quality reads, and the assemblies.

Furthermore, development is underway to include binning of contigs to
produce metagenome-assembled genomes (MAGs), predict their completeness
and contamination and assign taxonomy based in [GTDB](https://gtdb.ecogenomic.org/).
