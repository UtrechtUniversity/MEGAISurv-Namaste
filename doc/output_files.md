# Output files

Namaste includes a handful of processing steps, each using a few
bioinformatics tools to get the job done. All of these tools generate
output files, yielding a plethora of files in your `results/` directory!

Below is a description of the most relevant files.

## 1. Preprocessing

As described in [preprocess](preprocess.md), raw reads are quality filtered,
producing a file with high-quality reads and a report.

The quality-filtered reads are written to the directory `results/filtered_reads/`,
for example:

```bash
$ ls results/filtered_reads/
ERR10188468.fastq.gz  SRR16329143_1.fastq.gz  SRR17232957_1.fastq.gz
```

Each of these files will be slightly smaller than the original, raw read file.

For each sample, there will be a JSON and an HTML report of the filtering
process, of which the 'summary' part is parsed for later use. These
files are saved under `results/read_qc/` and `results/read_qc/summary/`.

```bash
$ ls results/read_qc/summary/
ERR10188468.json  SRR16329143_1.json  SRR17232957_1.json
```

As an example, these files look like:

```bash
$ cat results/read_qc/summary/ERR10188468.json
{
        "summary": {
                "fastplong_version": "0.2.2",
                "before_filtering": {
                        "total_reads":1295276,
                        "total_bases":4898452401,
                        "q20_bases":1985830150,
                        "q30_bases":644728976,
                        "q20_rate":0.405399,
                        "q30_rate":0.131619,
                        "read_mean_length":3781,
                        "gc_content":0.411915
                },
                "after_filtering": {
                        "total_reads":490786,
                        "total_bases":1801486902,
                        "q20_bases":899296433,
                        "q30_bases":331205707,
                        "q20_rate":0.499197,
                        "q30_rate":0.183851,
                        "read_mean_length":3670,
                        "gc_content":0.388144
                }
        }
}
```

Finally, read quality control (QC) statistics are collected for all samples and
stored in a CSV format table: `results/read_qc_summary.csv`. It summarises
the number of reads and basepairs before and after quality filtering, as well
as some other quality indicators:

```bash
$ head -2 results/read_qc_summary.csv
sample,before_total_reads,before_total_bases,before_q20_bases,before_q30_bases,before_q20_rate,before_q30_rate,before_read_mean_length,before_gc_content,after_total_reads,after_total_bases,after_q20_bases,after_q30_bases,after_q20_rate,after_q30_rate,after_read_mean_length,after_gc_content
SRR17232957_1,315586,1285925859,769615133,342928711,0.598491,0.266678,4074,0.507556,263485,1083609370,696137892,320329347,0.642425,0.295613,4112,0.50747
```

## 2. Assembly

The assembly itself generates sequences in fasta format along with a
tab-separated table with basic statistics per contig. These are stored under
`results/assembly/{sample}/` where {sample} indicates the base file name of
the raw reads file that was used to generate it. As described in
[assembly](assembly.md), Namaste uses Flye to assemble reads, which generates
many different output files. For convenience, we only show the relevant ones.

```bash
$ ls -F results/assembly/SRR17232957_1/
assembly.fasta  assembly_info.txt    mapped_back/

$ head -2 results/assembly/SRR17232957_1/assembly_info.txt
#seq_name       length  cov.    circ.   repeat  mult.   alt_group       graph_path
contig_1702     65004   6       N       N       1       *       *,1702,*
```

The `assembly_info.txt` file contains useful information such as contig length,
depth of coverage (cov.) and whether the contig is circular (circ.).
Output files from Flye are also described
[on its GitHub page](https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md#-flye-output).

Furthermore, Namaste maps reads back to contigs, which can be found in the
subdirectory `mapped_back/`.

```bash
$ ls results/assembly/SRR17232957_1/mapped_back/
SRR17232957_1.bam  SRR17232957_1-coverage.tsv

$ head -2 results/assembly/SRR17232957_1/mapped_back/SRR17232957_1-coverage.tsv
#rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
contig_1        1       17948   33      17948   100     4.92378 22.1    57.5
```

Per sample coverage files are generated with `samtools coverage` and are
slightly different from the `assembly_info.txt` by Flye in terms of
values reported (depth of coverage = meandepth).

Finally, these per sample coverage files are combined into one file for all
samples: `results/contig_coverage.tsv`.

```bash
$ head -2 results/contig_coverage.csv
Sample,Contig,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq
SRR17232957_1,contig_1,1,17948,33,17948,100,4.92378,22.1,57.5
```

And the coverage statistics of both Flye and samtools are combined in one
overall statistics file: `results/assembly_stats-concatenated.tsv.gz`.

```bash
$ zless results/assembly_stats-concatenated.tsv.gz | head -3
sample  contig  contig_length   contig_depth_flye       circular        contig_mapped_percent   contig_depth_samtools   contig_mapped_reads
SRR17232957_1   contig_1702     65004   6       N       100     6.54354 207
SRR17232957_1   contig_515      63366   6       N       100     5.86153 150
```

## 3. Antibiotic resistance genes

[Screening](arg_screening.md) of antibiotic resistance genes (ARGs) is done
with KMA, the output files of which are briefly described on
[its GitHub page](https://github.com/genomicepidemiology/kma#result-explanation).
Namaste uses the files with `.frag.gz` and `.res` extension, which are tabular
files of hits per contig and per sample, respectively.

Files for the contig-based screening of ARGs are stored in
`results/resistance_genes/`. Note that the `.frag.gz` file includes complete
sequences, which quickly becomes unwieldy when using contigs as input.

```bash
$ zless results/resistance_genes/SRR17232957_1.hmm.frag.gz | head -2
GTATCTGAACAACTATCTTGTGTGGAATAACCTTGTAAATTACGCCAAAGAAAGCGACATGGAGAAAAGGAACATCTTCTTAACTTTCGTTTTGGCAACATTGAAAACTGCTAAATGCAGAGATTTATCAAACAGACCAGCAGTTCCTCTGGTCGCCTAATTAGAATTTGTGGAGATGATAAGATGGTCAATATAACAGATGTAAAACAGATTCTTCAATTTGCAATAGATGCGGAGATTAAAGTCTTTCTTGATGGTGGCTGGGGTGTAGATGCTCTTCTTGGATATCAGTCAAGAGCCCATAATGATATTGACATTTTTGTAGAAAAGAACGATTATCAGAACTTTATAGAAATAATGAAAGCTAATGGCTTTTATGAGATTAAGATGGAATATACAACATTGAACCATACTGTATGGGAAGATTTGAAAAACAGAATTATTGATTTGCATTGTTTTGAATATACGGACGAAGGTGAAATTCTTTATGATGGGGATTGTTTTCCGGTAGAAACTTTTTCGGGTAAAGGAAGAATTGAGGAAATAGAGGTTTCCTGTATTGAACCATATAGTCAAGTAATGTTCCATCTGGGATACGAGTTTGATGAAAATGATGCACATGATGTGAAGTTATTGTGTGAGACACTTCATATCGAAATTCCAAATGAGTATAGATAACTGCAAATAACAAGTTTGTAGGGGAGTTCTGATACTCCCCTATAAAAATGGCTATTTATCAACTGTTTGTTGTGACATAGCCTTTCATATAATATATCTTTGCATGGATTCCACGGGCCTTACTCCGTGTGTTTTTGCGAATCCACGATAACTCTTCAGGTCGCTCGACCTTCAAAATTACCATGC        1       471     0       495     lnu(C)_1_AY928180       contig_418      4000    4864
AAGATGGTCAATATAACAGATGTAAAACAGATTCTTCAATTTGCAATAGATGCGGAGATTAAAGTCTTTCTTGATGGTGGCTTCTGGGGTGTAGATACTCTTCTTGGATATCAGTCAAGAGCCCATAATGATATTACGACATTTTGTAGAAAAGAACGATTATAAAGGGCACTTTATAGAAATAATGAAGCTAATGGCTTTTATGATTAAGATGGAATATACAACATTGAACCATACTGTATGGGAAGATTTGAAAAACAGACAATGATTTGCATTGTTTTGAATATACGGACGAAGGTGAAATTCTTTATGATGGGGATTGTTTTCCGGTAGAAACTTTTTCGAGGTAAAGGGAAGAATTGAGGAAATAGGGGTTTCCTGTATTGAACCATATAGTCAAGTAATGTTCCATCTGGGATACGGAAGTTTGATGAAAATGATGCACATGATGTGAAGTTATTGTGTGAGACACTTCATATCGAAATTCCCAATGAGTATAGATAACTACCAAATAACAAGTTTGTAGGGGAGTTCTGATACTCCCCTATAAAGAAAACAACATTTATCAACTGTTTGTGTGAGCACAGCCTAAAAAAAATACCTTGCCATCGACAGGGGTATTTTTGTTGACAGTAGCAATACTTAATCTTTGTATATTTCGTGGAGCTCTACCCTCGTAGATAATCGGGGCGTAGCGC  1       406     0       495     lnu(C)_1_AY928180       contig_1422     15104   15802

$ head -2 results/resistance_genes/SRR17232957_1.hmm.res
#Template       Score   Expected        Template_length Template_Identity       Template_Coverage       Query_Identity  Query_Coverage  Depth   q_value p_value
lnu(C)_1_AY928180            901               1             495           99.19          102.42           96.84           97.63            2.02          896.92        1.0e-26
```

Contigs in which the ARGs are masked are written back to the assembly folder:
`results/assembly/{sample}/assembly_ARG_masked.fasta`.

For reads, the results are written to `results/resistance_genes-reads/`.

```bash
$ head -2 results/resistance_genes-reads/SRR17232957_1.hmm.res
#Template       Score   Expected        Template_length Template_Identity       Template_Coverage       Query_Identity  Query_Coverage  Depth   q_value p_value
lnu(C)_1_AY928180           7599              21             495           99.19          100.00           99.19          100.00           19.28         7535.38        1.0e-26

$ zless results/resistance_genes-reads/SRR17232957_1.hmm.frag.gz | head -2
AAAGATGGTCAATATAACAGATGTAAAACAGATTCTTCAATTTGCAATAGATGCGGAGATTAAAGTCTTTCTTGATGGTGGCTTCTGGGGTGTAGATACTCTTCTTGGATATCAGTCAAGAGCCCATAATGATATTACGACATTTTGTAGAAAAGAACGATTATAAAGGGCACTTTATAGAAATAATGAAGCTAATGGCTTTTATGATTAAGATGGAATATACAACATTGAACCATACTGTATGGGAAGATTTGAAAAACAGACAATGATTTGCATTGTTTTGAATATACGGACGAAGGTGAAATTCTTTATGATGGGGATTGTTTTCCGGTAGAAACTTTTTCGAGGTAAAGGGAAGAATTGAGGAAATAGGGGTTTCCTGTATTGAACCATATAGTCAAGTAATGTTCCATCTGGGATACGGAAGTTTGATGAAAATGATGCACATGATGTGAAGTTATTGTGTGAGACACTTCATATCGAAATTCCCAATGAGTATAGATAACTACCAAATAACAAGTTTGTAGGGGAGTTCTGATACTCCCCTATAAAGAAAACAACATTTATCAACTGTTTGTGTGAGCACAGCCTAAAAAAAATACCTTGCCATCGACAGGGGTATTTTTGTTGACAGTAGCAATACTTAATCTTTGTATATTTCGTGGAGCTCTA    1       406     0       495     lnu(C)_1_AY928180       SRR17232957.300718 300718/1     0       672
AACGGTGTTTCTACCAAGTATCTGAACAACTATCTTGTGTGGAATAACCTTGTAAATTACGCCAAAGAAAGCGACATGGAGAAAAGGAACATCTTCTTAACTTTCGTTTTGGCAACATTGAAAACTGCCAGCCGTGTAGATTGATCAAACAGACCAGCAGTTCCTCTGGTCGCTTAATTAGAATTTTGTAGATGATAAGATGGTCAATATAACAGATGTAAAACAGATTTTCAATTTGCAATAGATGCGGAGATTAAAGTCTTTCTTGCAGGCTGGGGTGGATGCTCTTCTGATATCAGTCGAGAGCCCATAATGATATTGACATTTTTGTAGAAAGAACGATTATCAGAACTTTATAGAAATAATGAAAGCTAATGGCTTTTTATGAGATTAAGATGGAATATACAATACGAACCATATCAGGCAAAGATTTTGAAAACAGAATTATTGATTTAAATGTTTTGAATATACGGACGAAGGTGAAATTCTTTATGATGGGGATTGTTTTCCGGTAGAAACTTTTCGGGTAAAGGAAGAATTGAGGAAATAGAGGTTTCCCGTATTGAACCATAGTCAATGGGCTGTTCCATCTGGGATACGAGTTTGATGAAAATGATGCCATGATGTACAGTTATTGTGAGACACTTCATATCGAAATTCCAATGAGGAGTGATAACTGCAAATAACAAGTTTGTAGGGGAGTTCTGATACCCCCTATAAAAATGGCTATTTATCAATCGTTGGTTGTGACATAGCCTTCATATAATATATCTTTGCATGGATTCCACGGCCTTACCCGTGTGTTTTTGCGAATCCACGATAACTCTTCAGGTCGCCTGACCTTCAAAAATTACCATGCAAAAA        1       341     0       495     lnu(C)_1_AY928180       SRR17232957.245932 245932/1     2528    3392
```

And a summary of the read-based ARG screening approach for all samples is
stored as TSV table: `results/resistance_genes-reads/summary.tsv`.

```bash
$ head -2 results/resistance_genes-reads/summary.tsv
Sample  Gene    Score   Expected        Template_length Template_Identity       Template_Coverage       Query_Identity  Query_Coverage  Depth   q_value p_value
SRR17232957_1   lnu(C)_1_AY928180       7599    21      495     99.19   100     99.19   100     19.28   7535.38 1e-26
```

And a complete resistome dataset with ARG info, assembly statistics, taxonomy,
and plasmid/virus predictions (CSV table with 39 columns) is saved as:
`results/assembly_database.csv.gz`.

```bash
$ zless results/assembly_database.csv.gz | head -3
sample,arg,arg_identity,arg_length,arg_depth,contig,start_position,end_position,hit_length,arg_short,arg_coverage,Class,Phenotype,Resistance_mechanism,medical_importance,mean_assembly_depth_flye,mean_assembly_depth_samtools,total_assembly_length,contig_length,contig_depth_flye,circular,contig_mapped_percent,contig_depth_samtools,contig_mapped_reads,Species,Taxonomy,Species_strict,Taxonomy_strict,chromosome_score,plasmid_score,virus_score,plasmid_topology,plasmid_genes,conjugation_genes,amr_genes,virus_topology,virus_genes,virus_taxonomy,genomad_prediction
SRR17232957_1,lnu(C)_1_AY928180,99.19,495,2.02,contig_418,4000,4864,496,lnu(C),100.20202020202021,Lincosamide,Lincomycin,Enzymatic inactivation,HIA,5.768397033656589,7.170406012549915,18310891,22412,5,N,100,6.5833,67,Campylobacter coli,unclassified cellular organisms superkingdom;Campylobacterota;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter;Campylobacter coli,Campylobacter coli,unclassified cellular organisms superkingdom;Campylobacterota;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter;Campylobacter coli,0.7905,0.1618,0.0477,NA,NA,NA,NA,NA,NA,NA,chromosome
SRR17232957_1,lnu(C)_1_AY928180,99.19,495,2.02,contig_1422,15104,15802,496,lnu(C),100.20202020202021,Lincosamide,Lincomycin,Enzymatic inactivation,HIA,5.768397033656589,7.170406012549915,18310891,15802,4,N,100,4.42469,19,unclassified cellular organisms species,unclassified cellular organisms superkingdom;unclassified cellular organisms phylum;unclassified cellular organisms class;unclassified cellular organisms order;unclassified cellular organisms family;unclassified cellular organisms genus;unclassified cellular organisms species,,;;;;;;,0.9389,0.0446,0.0166,NA,NA,NA,NA,NA,NA,NA,chromosome
```

## 4. Antibiotic resistance mutations

Namaste [screens](arm_screening.md) metagenomes for the presence of antibiotic
resistance-conferring mutations using MetaPointFinder.
Its output is explained
[on GitHub](https://github.com/aldertzomer/metapointfinder/tree/main#output-files).

Results are stored under `results/resistance_mutations/{sample}/`:

```bash
$ ls results/resistance_mutations/SRR17232957_1/
accession               SRR17232957_1.class.dna.summary.txt                            SRR17232957_1.frag.gz                SRR17232957_1.prot.input.tsv
class                   SRR17232957_1.class.prot.summary.txt                           SRR17232957_1.gene.dna.summary.txt   SRR17232957_1.prot.updated_table_with_scores_and_mutations.tsv
dna_accession           SRR17232957_1.dna.input.tsv                                    SRR17232957_1.gene.prot.summary.txt  SRR17232957_1.res
dna_class               SRR17232957_1.dna.updated_table_with_scores_and_mutations.tsv  SRR17232957_1.log
matched_to_contigs.tsv  SRR17232957_1.error                                            SRR17232957_1.prot.hits.txt
```

These results are combined in the end to make one big table with all mutations
and contig annotations from the other analysis steps:
`results/mutation_database.csv.gz`. (A table with 33 columns.)

```bash
$ zless results/mutation_database.csv.gz | head -3
sample,contig,class,gene,mutations,Unknown,Wildtype,Resistant,fraction_resistant,mean_assembly_depth_flye,mean_assembly_depth_samtools,total_assembly_length,contig_length,contig_depth_flye,circular,contig_mapped_percent,contig_depth_samtools,contig_mapped_reads,Species,Taxonomy,Species_strict,Taxonomy_strict,chromosome_score,plasmid_score,virus_score,plasmid_topology,plasmid_genes,conjugation_genes,amr_genes,virus_topology,virus_genes,virus_taxonomy,genomad_prediction
ERR10188468,contig_20,MACROLIDE,NZ_CP018347.1@23S_ribosomal_RNA@23S:32347-35250,None,1,NA,0,0,7.500926081831958,7.098058050597744,115259403,370730,26,N,99.9957,24.7572,2715,Streptococcus sp. LPB0220,unclassified cellular organisms superkingdom;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Streptococcus sp. LPB0220,Streptococcus sp. LPB0220,unclassified cellular organisms superkingdom;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Streptococcus sp. LPB0220,0.8733,0.0833,0.0434,NA,NA,NA,NA,NA,NA,NA,chromosome
ERR10188468,contig_20,MACROLIDE_STREPTOGRAMIN,NZ_CP018347.1@23S_ribosomal_RNA@23S:32347-35250,None,1,NA,0,0,7.500926081831958,7.098058050597744,115259403,370730,26,N,99.9957,24.7572,2715,Streptococcus sp. LPB0220,unclassified cellular organisms superkingdom;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Streptococcus sp. LPB0220,Streptococcus sp. LPB0220,unclassified cellular organisms superkingdom;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Streptococcus sp. LPB0220,0.8733,0.0833,0.0434,NA,NA,NA,NA,NA,NA,NA,chromosome
```

## 5. Taxonomic classifications

Namaste [classifies](taxonomy.md) contigs with Centrifuger and its database of
microbiota + human. Centrifuger produces a table in TSV format and has a
`quant` function to calculate microbiota profiles. Furthermore, Namaste uses
both default and 'strict' settings and assigns scientific names to each using
TaxonKit. All of these are stored as TSV files in
`results/taxonomic_classification/{sample}/`.

For details on the file format, consult
[Centrifuger's documentation](https://github.com/mourisl/centrifuger#inputoutput).

```bash
$ ls results/taxonomic_classification/SRR17232957_1/
centrifuger_masked-quant.tsv  centrifuger_masked-strict-quant.tsv  centrifuger_masked-strict+taxa.tsv  centrifuger_masked-strict.tsv  centrifuger_masked+taxa.tsv  centrifuger_masked.tsv

$ head -3 results/taxonomic_classification/SRR17232957_1/*
==> results/taxonomic_classification/SRR17232957_1/centrifuger_masked-quant.tsv <==
name    taxID   taxRank genomeSize      numReads        numUniqueReads  abundance
root    1       no rank 2070143 1696    1018    1.000000
Bacteria        2       superkingdom    4534370 1628    989     0.911828

==> results/taxonomic_classification/SRR17232957_1/centrifuger_masked-strict-quant.tsv <==
name    taxID   taxRank genomeSize      numReads        numUniqueReads  abundance
root    1       no rank 2070143 39      14      1.000000
Bacteria        2       superkingdom    4534370 39      14      1.000000

==> results/taxonomic_classification/SRR17232957_1/centrifuger_masked-strict+taxa.tsv <==
readID  seqID   taxID   score   2ndBestScore    hitLength       queryLength     numMatches
contig_1        unclassified    0       0       0       0       17948   1               ;;;;;;
contig_10       unclassified    0       0       0       0       5781    1               ;;;;;;

==> results/taxonomic_classification/SRR17232957_1/centrifuger_masked-strict.tsv <==
readID  seqID   taxID   score   2ndBestScore    hitLength       queryLength     numMatches
contig_1        unclassified    0       0       0       0       17948   1
contig_10       unclassified    0       0       0       0       5781    1

==> results/taxonomic_classification/SRR17232957_1/centrifuger_masked+taxa.tsv <==
readID  seqID   taxID   score   2ndBestScore    hitLength       queryLength     numMatches
contig_1        NZ_CP035464.1   78344   121     64      26      17948   1       Bifidobacterium pullorum        unclassified cellular organisms superkingdom;Actinomycetota;Actinomycetes;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;Bifidobacterium pullorum
contig_10       NZ_CP102287.1   679190  100     0       25      5781    1       Hoylesella buccalis     unclassified cellular organisms superkingdom;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Hoylesella;Hoylesella buccalis

==> results/taxonomic_classification/SRR17232957_1/centrifuger_masked.tsv <==
readID  seqID   taxID   score   2ndBestScore    hitLength       queryLength     numMatches
contig_1        NZ_CP035464.1   78344   121     64      26      17948   1
contig_10       NZ_CP102287.1   679190  100     0       25      5781    1
```

Also, for each of these Namaste produces microbiota profiles as TSV in
`results/microbiota_profile/`.

```bash
$ ls results/microbiota_profile/SRR17232957_1-*
results/microbiota_profile/SRR17232957_1-per_contig.tsv   results/microbiota_profile/SRR17232957_1-strict-per_contig.tsv
results/microbiota_profile/SRR17232957_1-per_species.tsv  results/microbiota_profile/SRR17232957_1-strict-per_species.tsv

$ head -2  results/microbiota_profile/SRR17232957_1-*
==> results/microbiota_profile/SRR17232957_1-per_contig.tsv <==
readID  seqID   taxID   score   2ndBestScore    hitLength       queryLength     numMatches      Species Lineage length  depth   total_bases     grand_total     Sample  Percentage
contig_282      NZ_LN879430.1   1679721 66113   52062   295     22345   1       Herbinix luporum        unclassified cellular organisms superkingdom;Bacillota;Clostridia;Lachnospirales;Lachnospiraceae;Herbinix;Herbinix luporum  22345   61.431  1372675.6949999998      132539452.71505 SRR17232957_1   1.0356732783189853

==> results/microbiota_profile/SRR17232957_1-per_species.tsv <==
Sample  Species Lineage total_bases     Percentage
SRR17232957_1   unclassified cellular organisms species unclassified cellular organisms superkingdom;unclassified cellular organisms phylum;unclassified cellular organisms class;unclassified cellular organisms order;unclassified cellular organisms family;unclassified cellular organisms genus;unclassified cellular organisms species        8955780.64689   6.757067773732448

==> results/microbiota_profile/SRR17232957_1-strict-per_contig.tsv <==
readID  seqID   taxID   score   2ndBestScore    hitLength       queryLength     numMatches      Species Lineage length  depth   total_bases     grand_total     Sample  Percentage
contig_282      family  186803  59536   59536   259     22345   1       unclassified Lachnospiraceae species    unclassified cellular organisms superkingdom;Bacillota;Clostridia;Lachnospirales;Lachnospiraceae;unclassified Lachnospiraceae genus;unclassified Lachnospiraceae species    22345   61.431  1372675.6949999998      132539452.71505 SRR17232957_1   1.0356732783189853

==> results/microbiota_profile/SRR17232957_1-strict-per_species.tsv <==
Sample  Species Lineage total_bases     Percentage
SRR17232957_1   NA      ;;;;;;  123421656.91637 93.1206930375044
```

In the end, per contig classifications of all samples are collected in two files:
`results/classifications-concatenated.tsv.gz` and
`results/strict_classifications-concatenated.tsv.gz`.

## 6. Plasmid predictions

Namaste uses geNomad to [predict](genomad.md) whether contigs derive from
chromosomal DNA, plasmids or viruses. Per sample results are saved in
`results/plasmid_prediction/{sample}/`.

```bash
$ ls -F results/plasmid_prediction/SRR17232957_1/
assembly_aggregated_classification/     assembly_annotate.log         assembly_marker_classification/     assembly_nn_classification.log  assembly_summary/
assembly_aggregated_classification.log  assembly_find_proviruses/     assembly_marker_classification.log  assembly_score_calibration/     assembly_summary.log
assembly_annotate/                      assembly_find_proviruses.log  assembly_nn_classification/         assembly_score_calibration.log
```

A detailed description of the outputs created by geNomad can be found in
[its documentation](https://portal.nersc.gov/genomad/quickstart.html#understanding-the-outputs).

Output for contigs predicted to derive from plasmids and viruses are
concatenated for all samples and saved under 'results':
`results/plasmid_predictions-concatenated.tsv.gz` and
`results/virus_predictions-concatenated.tsv.gz`.

```bash
$ zless results/plasmid_predictions-concatenated.tsv.gz | head -3
contig  length  plasmid_topology        plasmid_genes   genetic_code    plasmid_score   fdr     n_hallmarks     marker_enrichment       conjugation_genes       amr_genes       sample
contig_1753     4328    No terminal repeats     10      11      0.999   0.001   0       2.1513  NA      NA      SRR17232957_1
contig_1750     7896    No terminal repeats     11      11      0.999   0.001   0       0.8884  NA      NA      SRR17232957_1

$ zless results/virus_predictions-concatenated.tsv.gz | head -3
contig  length  virus_topology  coordinates     virus_genes     genetic_code    virus_score     fdr     n_hallmarks     marker_enrichment       virus_taxonomy  sample
contig_941      8070    No terminal repeats     NA      13      11      0.9989  0.0011  4       13.2017 Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;       SRR17232957_1
contig_837      12674   No terminal repeats     NA      20      11      0.9989  0.0011  6       17.1828 Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;       SRR17232957_1
```
