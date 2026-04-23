# Taxonomic classification of contigs

A metagenomic analysis is not complete without taxonomic assignments:
we want to know what species of bacteria are present. Especially those
harbouring antibiotic resistance genes (ARGs).
Namaste takes the [assembled contigs](assembly.md) from which the
[ARGs are masked](arg_screening.md) and classifies those using
[Centrifuger](https://github.com/mourisl/centrifuger) (version 1.0.6)
with the provided database for human, bacteria, viruses and archea:
[cfr_hpv+gbsarscov2](https://zenodo.org/records/10023239).

Centrifuger assigns taxonomy to each contig, based on the lowest
common ancestor of all valid matches in the database.
(See details in [the paper](https://doi.org/10.1186/s13059-024-03244-4).)
By default, it also consider matches with a short hit length, opening
up the possibility of finding matches based on sequencing barcodes that
remain in the input reads and reference database.
To avoid these potential, unwanted hit, Namaste also runs Centrifuger with
a stricter increased minimum hit length of 100bp (option `--min-hitlen 100`).

Centrifuger assigns the taxonomy as taxon ID or 'taxid', which is a number
that corresponds to a taxon from
[NCBI's taxonomy tree](https://www.ncbi.nlm.nih.gov/taxonomy).
To convert this to the scientific name, Namaste uses
[taxonkit reformat](https://bioinf.shenwei.me/taxonkit) (version 0.18.0)
with [NCBI's 'taxdump'](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) file.
and returns names for taxonomic ranks:

- kingdom
- phylum
- class
- order
- family
- genus
- species

These are appended to the tabular output file in a human and machine friendly
TSV format.
Final classifications for each sample are saved per contig, and also
converted to per species tables - each in tab-separated (TSV) format.
This information is also used in the final, overall outputs:
[assembly database](arg_screening.md) and [mutation database](arm_screening.md).

## Output files

Taxonomic classification consists of three steps:

1. Taxon ID assignment by Centrifuger
2. Conversion to scientific name by Taxonkit
3. Calculating per sample microbiota profiles in R (per contig and per species)

Each of these steps is carried out twice: once for Centrifuger's default
settings, and once for the stricter 'minimum hit length = 100' setting.

(Again, note that the ARG-masked cnotigs sequences are used as input.)

```text
results/
  taxonomic_classification/
    {sample}/
      centrifuger_masked.tsv                # Taxonomic classification by Centrifuger
      centrifuger_masked-strict.tsv         # Taxonomic classification with strict setting
      centrifuger_masked+taxa.tsv           # Classification with scientific names by Taxonkit
      centrifuger_masked-strict+taxa.tsv    # Classification + names with strict setting
    microbiota_profile/
      {sample}-per_contig.tsv               # Per sample microbiota profile, summarised per contig
      {sample}-strict-per_contig.tsv        # Per sample and contig profile with strict setting
      {sample}-per_species.tsv              # Per sample microbiota profile, summarised per species
      {sample}-strict-per_species.tsv"      # Per sample and species profile with strict setting
```

For details, please see [output](output_files.md).
