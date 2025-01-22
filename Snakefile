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


INPUT_DIR = config["input_directory"]
SAMPLES = config["samples"]

### Step 2: Specify output files ###


rule all:
    input:
        expand("assemblies/Quadram/{sample}/assembly.fasta", sample=SAMPLES),
        expand(
            "resistance_genes/{sample}.hmm.{suffix}",
            sample=SAMPLES,
            suffix=["aln", "frag.gz", "fsa", "res"],
        ),
        expand(
            "classifications/Quadram/{sample}/centrifuger_masked{extension}",
            sample=SAMPLES,
            extension=[".tsv", "-quant.tsv", ".kreport.txt"]
        ),


### Step 3: Define processing steps that generate the output ###


rule metagenomic_assembly:
    input:
        INPUT_DIR + "{sample}.fastq.gz",
    output:
        "assemblies/Quadram/{sample}/assembly.fasta",
    params:
        output_dir="assemblies/Quadram/{sample}",
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
        "assemblies/Quadram/{sample}/assembly.fasta",
    output:
        aln="resistance_genes/{sample}.hmm.aln",
        frag="resistance_genes/{sample}.hmm.frag.gz",
        fsa="resistance_genes/{sample}.hmm.fsa",
        res="resistance_genes/{sample}.hmm.res",
    params:
        db=config["kma"]["database"],
        prefix="resistance_genes/{sample}.hmm",
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
        frag="resistance_genes/{sample}.hmm.frag.gz",
        assembly="assemblies/Quadram/{sample}/assembly.fasta",
    output:
        gene_locations="resistance_genes/{sample}.locations.txt",
        masked_assembly="assemblies/Quadram/{sample}/assembly_ARG_masked.fasta",
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
        fasta="assemblies/Quadram/{sample}/assembly_ARG_masked.fasta",
        db="/mnt/data/db/centrifuger/cfr_hpv+gbsarscov2.1.cfr",
    output:
        tsv="classifications/Quadram/{sample}/centrifuger_masked.tsv",
        quant="classifications/Quadram/{sample}/centrifuger_masked-quant.tsv",
        report="classifications/Quadram/{sample}/centrifuger_masked.kreport.txt",
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
centrifuger-kreport -x {params.db} {output.tsv} > {output.report}
        """
