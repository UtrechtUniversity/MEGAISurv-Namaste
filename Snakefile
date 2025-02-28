"""
Author: Sam Nooij
Organisation: Utrecht University
Department: Clinical Infectiology (KLIF), Infectious Diseases & Immunology,
  Biomolecular Health Sciences, Faculty of Veterinary Medicine
Date: 2025-01-17

Workflow for assembling MinION long-read metagenomes,
taxonomic classification and identification of antibiotic
resistance genes, based on previous work by Aldert Zomer
(from the same department).

Input: (Gzipped) Fastq files of relevant metagenomic samples
Output: (various)

Example use:
    $ snakemake --profile config

N.B. Variables are set in the configuration files under `config`.
"""

from pathlib import Path
import functools
import operator

### Step 1: Import configuration file ###


configfile: Path("config/parameters.yaml")


# Read the input directory and automatically discover
# input files with '.fastq.gz' extension.
INPUT_DIR = Path(config["input_directory"])
INPUT_FILES = list(INPUT_DIR.glob("*.fastq.gz"))
SAMPLES = [file.stem.replace(".fastq", "") for file in INPUT_FILES]

OUTPUT_DIR = config["output_directory"]

### Step 2: Specify output files ###


rule all:
    input:
        expand(OUTPUT_DIR + "flye/{sample}/assembly.fasta", sample=SAMPLES),
        expand(
            OUTPUT_DIR + "kma/{sample}.hmm.{suffix}",
            sample=SAMPLES,
            suffix=["aln", "frag.gz", "fsa", "res"],
        ),
        expand(
            OUTPUT_DIR + "centrifuger/{sample}/centrifuger_masked{extension}",
            sample=SAMPLES,
            extension=[".tsv", "-quant.tsv"],
        ),


### Step 3: Define processing steps that generate the output ###


rule metagenomic_assembly:
    input:
        INPUT_DIR / "{sample}.fastq.gz",
    output:
        OUTPUT_DIR + "flye/{sample}/assembly.fasta",
    params:
        output_dir=OUTPUT_DIR + "flye/{sample}",
        settings="--meta",
    conda:
        "envs/flye.yaml"
    threads: config["flye"]["threads"]
    log:
        "log/flye/{sample}.txt",
    benchmark:
        "log/benchmark/flye/{sample}.txt"
    shell:
        """
flye {params.settings} --threads {threads} --nano-hq {input}\
 --out-dir {params.output_dir}
        """


rule screen_antibiotic_resistance_genes:
    input:
        OUTPUT_DIR + "flye/{sample}/assembly.fasta",
    output:
        aln=OUTPUT_DIR + "kma/{sample}.hmm.aln",
        frag=OUTPUT_DIR + "kma/{sample}.hmm.frag.gz",
        fsa=OUTPUT_DIR + "kma/{sample}.hmm.fsa",
        res=OUTPUT_DIR + "kma/{sample}.hmm.res",
    params:
        db=config["kma"]["database"],
        prefix=OUTPUT_DIR + "kma/{sample}.hmm",
    conda:
        "envs/kma.yaml"
    threads: config["kma"]["threads"]
    log:
        "log/kma/{sample}.txt",
    benchmark:
        "log/benchmark/kma/{sample}.txt"
    shell:
        """
kma -t {threads} -bcNano -t_db {params.db} -i {input} -o {params.prefix} -ont -hmm
        """


rule mask_resistance_gene_positions:
    input:
        frag=OUTPUT_DIR + "kma/{sample}.hmm.frag.gz",
        assembly=OUTPUT_DIR + "flye/{sample}/assembly.fasta",
    output:
        gene_locations=OUTPUT_DIR + "kma/{sample}.locations.txt",
        masked_assembly=OUTPUT_DIR + "flye/{sample}/assembly_ARG_masked.fasta",
    conda:
        "envs/bedtools.yaml"
    threads: 1
    log:
        "log/kma/{sample}-positions.txt",
    benchmark:
        "log/benchmark/kma/{sample}-positions.txt"
    shell:
        """
# Cut contig ID, start and stop positions from file
zcat {input.frag} | cut -f 7-9 | sort | uniq > {output.gene_locations}

maskFastaFromBed -fi {input.assembly} -bed {output.gene_locations} -fo {output.masked_assembly}
        """


rule taxonomic_classification:
    input:
        fasta=OUTPUT_DIR + "flye/{sample}/assembly_ARG_masked.fasta",
        db="/mnt/data/db/centrifuger/cfr_hpv+gbsarscov2.1.cfr",
    output:
        tsv=OUTPUT_DIR + "centrifuger/{sample}/centrifuger_masked.tsv",
        quant=OUTPUT_DIR + "centrifuger/{sample}/centrifuger_masked-quant.tsv",
    params:
        db="/mnt/data/db/centrifuger/cfr_hpv+gbsarscov2",
    conda:
        "envs/centrifuger.yaml"
    threads: config["centrifuger"]["threads"]
    log:
        "log/centrifuger/{sample}.txt",
    benchmark:
        "log/benchmark/centrifuger/{sample}.txt"
    shell:
        """
centrifuger -u {input.fasta} -t {threads} -x {params.db} > {output.tsv}

centrifuger-quant -x {params.db} -c {output.tsv} > {output.quant}
        """
