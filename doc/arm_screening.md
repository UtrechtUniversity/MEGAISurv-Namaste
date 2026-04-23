# Screen antibiotic resistance mutations

Metagenomic reads are screened for the presence of mutations that confer
resistance to certain antibiotics. [High-quality reads](preprocess.md)
are mapped to a database of genes and described resistance mutations using
[MetaPointFinder](https://github.com/aldertzomer/metapointfinder) (version 1.01).
In short, MetaPointFinder uses nucleotide matching with the k-mer based tool
[KMA](https://github.com/genomicepidemiology/kma)
and aligns predicted protein (amino acid) sequences with
[DIAMOND](https://github.com/bbuchfink/diamond).
It then uses custom R scripts to score the matches and count the number
of reads with wild-type (WT), resistant (R) and unknown alleles/sequence
variants.

MetaPointFinder is run with a sequence identity threshold of 85% (default).

To match the detected variants on read-level to the assembled contigs,
the read IDs from hits reported by MetaPointFinder are compared to the
reads mapped back to contigs ([BAM file by minimap2](assembly.md)).
These are then combined to provide the number of WT and R variants
per contig and calculate the fraction of resistant mutants on the contig
level.

At the very end, resistance mutations are combined with assembly statistics,
taxonomic classifications and plasmid predictions into one summary table:
the 'mutation database'.

## Output files

Screening metagenomes for antibiotic resistance mutations yields a number
of intermediate results files as well as a file that matches hits to
assembled contigs:

Files starting with the sample name are the
[output of MetaPointFinder](https://github.com/aldertzomer/metapointfinder#output-files).

```text
results/
  resistance_mutations/
    {sample}/
      {sample}.dna.updated_table_with_scores_and_mutations.tsv
      {sample}.prot_updated_with_scores_and_mutations.tsv
      {sample}.class.dna.summary.txt
      {sample}.class.prot.summary.txt
      {sample}.dna.summary.txt
      {sample}.gene.prot.summary.txt
      ...
      matched_to_contigs.tsv                # table matching mutants to contigs
  mutation_database.csv.gz                  # Overall summary of mutations and
                                            #  other contig annotations
```

For details, please see [output](output_files.md).
