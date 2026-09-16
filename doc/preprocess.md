# Preprocessing of raw reads

Raw reads from long-read technologies (that is, Oxford Nanopore or PacBio
SMRT) are quality controlled with
[fastplong](https://github.com/OpenGene/fastplong) (version 0.2.2).
This checks the estimated basecalling accuracy (Phred scores) of each
read, removing reads that do not meet a quality threshold of 15
(the default) in a maximum of 40% of bases (also default).
Additionally, homopolymers are trimmed off the 3'-end of reads
(option `--trim_poly_X`), and the trimming is extended by 10 basepairs
(option `--trimming_extension 10` (default value)), to remove possible
artifacts near the reads ends, yielding slightly cleaner reads.

The number of reads and their lengths are collected before and after
quality filtering, so that the total sequence information (in basepairs)
can be calculated to provide some summary statistics.

## Output files

The preprocessing or quality filtering steps yield high-quality reads
and statistics files. These are written to the directories:

```txt
results/
  filtered_reads/
    {sample}.fastq.gz   # high-quality reads
  read_qc/
    {sample}.json       # report by fastplong, JSON format
    {sample}.html       # report by fastplong, HTML format
    summary/
      {sample}.json     # summarised JSON report
  read_qc_summary.csv   # overall QC summary in table (CSV) format
```

For details, please see [output](output_files.md).

## Next steps

&rarr; [Assembly](assembly.md)

&rarr; [Screen antibiotic resistance _mutations_](arm_screening.md)

&rarr; [Screen antibiotic resistance _genes_](arg_screening.md)
