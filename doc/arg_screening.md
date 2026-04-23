# Screen antibiotic resistance genes

For detecting antibiotic resistance genes (ARGs) in metagenomic data,
Namaste screens both [high-quality reads](preprocess.md) and
[assembled contigs](assembly.md) using the
[ResFinder database](https://bitbucket.org/genomicepidemiology/resfinder_db/)
and [KMA](https://github.com/genomicepidemiology/kma) (version 1.4.2).

The ResFinder database is a curated database of acquired resistance genes.
Therefore, it does not include genes causing intrinsic resistance in some
bacterial species, and genes that are embedded in the chromosome and are
not known to transfer between bacteria (by horizontal gene transfer).
In other words, Namaste focuses on resistance that is probably
transferrable between bacteria and may spread by, for example, plasmids and
transposons.

The ResFinder database is automatically downloaded before running, ensuring
the latest version is used. (At the time of writing, that is version 2.6.0,
last updated in May 2025.)

KMA is run using the preset options for Nanopore-derived sequences
(options `-bcNano -ont`) and using HMMs to assign templates and provide
additional alignment metrics (option `-hmm`).

In the basis, Namaste uses an assembly-based analysis approach where
assembled contigs are screened for ARGs. However, this may be less
sensitive than a read-based approach where no assembly step is involved
and all the reads may be assigned to ARGs. As KMA is quite fast,
Namaste includes ARG screening for both contigs and reads so that the
outcome may be compared and the relative sensitivity of each method
assessed.

ARGs are often associated with pathogenic species and their presence
in databases may biased by this. To avoid this bias in taxonomic classification,
Namaste masks ARG positions from contigs before classification. For this, it uses
[bedtools maskFastaFromBed](https://bedtools.readthedocs.io/en/latest/index.html)
(version 2.31.1).

At the very end, ARGs are combined with assembly statistics, taxonomic
classifications and plasmid predictions into one summary table:
the 'assembly database'.

## Output files

Screening metagenomes for antibiotic resistance genes yields a number
of intermediate results files of which the most relevant information
is summarised in the 'assembly database':

Files starting with the sample name are the output of KMA.

```text
results/
  assembly/
    assembly_ARG_masked.fasta       # assembly where ARGs are masked with Ns
  resistance_genes/
    {sample}.hmm.aln                # per sample KMA output files (contigs)
    {sample}.hmm.frag.gz
    {sample}.hmm.fra
    {sample}.hmm.res
  resistance_genes-reads/
    {sample}.hmm.aln                # per sample KMA output files (reads)
    {sample}.hmm.frag.gz
    {sample}.hmm.fra
    {sample}.hmm.res
    summary.tsv                     # Overall summary of read-based ARG screening
assembly_database.csv.gz            # Overall summary of ARGs and
                                    #  other contig annotations
```

For details, please see [output](output_files.md).
