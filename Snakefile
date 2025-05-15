"""
Author: Sam Nooij
Organisation: Utrecht University
Department: Clinical Infectiology (KLIF), Infectious Diseases & Immunology,
  Biomolecular Health Sciences, Faculty of Veterinary Medicine
Date: 2025-01-17

Workflow for assembling MinION long-read metagenomes, taxonomic classification
and identification of antibiotic resistance genes, based on previous work by
Aldert Zomer (from the same department). Now extended with raw read quality-
control and filtering.

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
        # Metagenomic assemblies (by Flye)
        expand(OUTPUT_DIR + "assembly/{sample}/assembly.fasta", sample=SAMPLES),
        # Assembly assessment reports (by metaQUAST)
        expand(OUTPUT_DIR + "quast/{sample}/metaquast.log", sample=SAMPLES),
        # Simple assembly statistics (by seqkit)
        OUTPUT_DIR + "assembly/assembly_statistics-seqkit.tsv",
        # Antibiotic resistance gene screening (by KMA)
        expand(
            OUTPUT_DIR + "kma/{sample}.hmm.{suffix}",
            sample=SAMPLES,
            suffix=["aln", "frag.gz", "fsa", "res"],
        ),
        # Taxonomic classifications (by Centrifuger)
        expand(
            OUTPUT_DIR + "centrifuger/{sample}/centrifuger_masked{extension}",
            sample=SAMPLES,
            extension=[".tsv", "-quant.tsv"],
        ),
        # Microbiota profiles (TaxonKit + custom R script)
        expand(
            OUTPUT_DIR + "microbiota_profiles/{sample}-per_{unit}.tsv",
            sample=SAMPLES,
            unit=["contig", "species"],
        ),
        # Chromosome/plasmid/virus predictions (geNomad)
        expand(
            OUTPUT_DIR
            + "genomad/{sample}/assembly_aggregated_classification/assembly_aggregated_classification.tsv",
            sample=SAMPLES,
        ),


### Step 3: Define processing steps that generate the output ###


rule read_quality_control:
    input:
        INPUT_DIR / "{sample}.fastq.gz",
    output:
        filtered=OUTPUT_DIR + "filtered/{sample}.fastq.gz",
        json=OUTPUT_DIR + "read_qc/{sample}.json",
        html=OUTPUT_DIR + "read_qc/{sample}.html",
    conda:
        "envs/fastplong.yaml"
    threads: config["fastplong"]["threads"]
    log:
        "log/read_qc/{sample}.txt",
    benchmark:
        "log/benchmark/read_qc/{sample}.txt"
    shell:
        """
fastplong --in {input} --out {output.filtered}\
 --json {output.json} --html {output.html}\
 --thread {threads} > {log} 2>&1
        """


rule metagenomic_assembly:
    input:
        OUTPUT_DIR + "filtered/{sample}.fastq.gz",
    output:
        assembly=OUTPUT_DIR + "assembly/{sample}/assembly.fasta",
        info=OUTPUT_DIR + "assembly/{sample}/assembly_info.txt",
    params:
        output_dir=OUTPUT_DIR + "assembly/{sample}",
        settings="--meta",
    conda:
        "envs/flye.yaml"
    threads: config["flye"]["threads"]
    log:
        "log/assembly/{sample}.txt",
    benchmark:
        "log/benchmark/assembly/{sample}.txt"
    shell:
        """
flye {params.settings} --threads {threads} --nano-hq {input}\
 --out-dir {params.output_dir} > {log} 2>&1
        """


rule assess_assembly:
    input:
        OUTPUT_DIR + "assembly/{sample}/assembly.fasta",
    output:
        report=OUTPUT_DIR + "quast/{sample}/report.html",
        icarus=OUTPUT_DIR + "quast/{sample}/icarus.html",
        log=OUTPUT_DIR + "quast/{sample}/metaquast.log",
    params:
        output_dir=OUTPUT_DIR + "quast/{sample}",
    conda:
        "envs/quast.yaml"
    threads: config["metaquast"]["threads"]
    log:
        "log/assess_assembly/{sample}.txt",
    benchmark:
        "log/benchmark/assess_assembly/{sample}.txt"
    shell:
        """
metaquast.py -o {params.output_dir} -t {threads} {input} > {log} 2>&1
        """


rule simple_assembly_statistics:
    input:
        expand(OUTPUT_DIR + "assembly/{sample}/assembly.fasta", sample=SAMPLES),
    output:
        OUTPUT_DIR + "assembly/assembly_statistics-seqkit.tsv",
    conda:
        "envs/seqkit.yaml"
    threads: config["simple_stats"]["threads"]
    log:
        "log/simple_assembly_statistics.txt",
    benchmark:
        "log/benchmark/simple_assembly_statistics.txt"
    shell:
        """
seqkit stats -Ta -j {threads} > {output} 2> {log}
        """


rule screen_antibiotic_resistance_genes:
    input:
        OUTPUT_DIR + "assembly/{sample}/assembly.fasta",
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
kma -t {threads} -bcNano -t_db {params.db} -i {input} -o {params.prefix}\
 -ont -hmm > {log} 2>&1
        """


rule mask_resistance_gene_positions:
    input:
        frag=OUTPUT_DIR + "kma/{sample}.hmm.frag.gz",
        assembly=OUTPUT_DIR + "assembly/{sample}/assembly.fasta",
    output:
        gene_locations=OUTPUT_DIR + "kma/{sample}.locations.txt",
        masked_assembly=OUTPUT_DIR + "assembly/{sample}/assembly_ARG_masked.fasta",
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

maskFastaFromBed -fi {input.assembly} -bed {output.gene_locations}\
 -fo {output.masked_assembly} > {log} 2>&1
        """


rule taxonomic_classification:
    input:
        fasta=OUTPUT_DIR + "assembly/{sample}/assembly_ARG_masked.fasta",
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
centrifuger -u {input.fasta} -t {threads} -x {params.db} > {output.tsv}\
 2> {log}

centrifuger-quant -x {params.db} -c {output.tsv} > {output.quant} 2> {log}
        """


rule lookup_taxids:
    input:
        OUTPUT_DIR + "centrifuger/{sample}/centrifuger_masked.tsv",
    output:
        OUTPUT_DIR + "centrifuger/{sample}/centrifuger_masked+taxa.tsv",
    threads: 1
    conda:
        "envs/taxonkit.yaml"
    log:
        "log/lookup_taxids/{sample}.txt",
    benchmark:
        "log/benchmark/lookup_taxids/{sample}.txt"
    shell:
        """
taxonkit reformat {input} -I 3\
 -f '{{s}}\t{{k}};{{p}};{{c}};{{o}};{{f}};{{g}};{{s}}' -F\
 > {output} 2> {log}
        """


rule generate_microbiota_profiles:
    input:
        assembly_info=OUTPUT_DIR + "assembly/{sample}/assembly_info.txt",
        classifications=OUTPUT_DIR + "centrifuger/{sample}/centrifuger_masked+taxa.tsv",
    output:
        per_contig=OUTPUT_DIR + "microbiota_profiles/{sample}-per_contig.tsv",
        per_species=OUTPUT_DIR + "microbiota_profiles/{sample}-per_species.tsv",
    params:
        sample="{sample}",
    conda:
        "envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/generate_profiles/{sample}.txt",
    benchmark:
        "log/benchmark/generate_profiles/{sample}.txt"
    script:
        "scripts/generate_profiles.R"


rule genomad:
    input:
        fasta=OUTPUT_DIR + "assembly/{sample}/assembly.fasta",
        db=config["genomad"]["database"],
    output:
        aggregated_classification=OUTPUT_DIR
        + "genomad/{sample}/assembly_aggregated_classification/assembly_aggregated_classification.tsv",
        plasmid_summary=OUTPUT_DIR
        + "genomad/{sample}/assembly_summary/assembly_plasmid_summary.tsv",
        virus_summary=OUTPUT_DIR
        + "genomad/{sample}/assembly_summary/assembly_virus_summary.tsv",
    params:
        output_dir=OUTPUT_DIR + "genomad/{sample}/",
    conda:
        "envs/genomad.yaml"
    threads: config["genomad"]["threads"]
    log:
        "log/genomad/{sample}.txt",
    benchmark:
        "log/benchmark/genomad/{sample}.txt"
    shell:
        """
genomad end-to-end -t {threads} --cleanup --enable-score-calibration\
 {input.fasta} {params.output_dir} {input.db} > {log} 2>&1
        """
