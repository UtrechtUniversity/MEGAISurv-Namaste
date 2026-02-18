## General helper functions: define input samples
from pathlib import Path

# Read the input directory and automatically discover
# input files with '.fastq.gz' extension.
INPUT_DIR = Path(config["input_directory"])
INPUT_FILES = list(INPUT_DIR.glob("*.fastq.gz"))
SAMPLES = [file.stem.replace(".fastq", "") for file in INPUT_FILES]


rule make_assembly_database:
    input:
        assembly_info=expand(
            "results/assembly/{sample}/assembly_info.txt", sample=SAMPLES
        ),
        arg_hits=expand("results/resistance_genes/{sample}.hmm.frag.gz", sample=SAMPLES),
        arg_results=expand("results/resistance_genes/{sample}.hmm.res", sample=SAMPLES),
        classification=expand(
            "results/taxonomic_classification/{sample}/centrifuger_masked+taxa.tsv",
            sample=SAMPLES,
        ),
        genomad_scores=expand(
            "results/plasmid_prediction/{sample}/assembly_aggregated_classification/assembly_aggregated_classification.tsv",
            sample=SAMPLES,
        ),
        genomad_plasmid=expand(
            "results/plasmid_prediction/{sample}/assembly_summary/assembly_plasmid_summary.tsv",
            sample=SAMPLES,
        ),
        genomad_virus=expand(
            "results/plasmid_prediction/{sample}/assembly_summary/assembly_virus_summary.tsv",
            sample=SAMPLES,
        ),
    output:
        assembly_stats="results/assembly_stats-concatenated.tsv.gz",
        taxonomic_classification="results/classifications-concatenated.tsv.gz",
        plasmid_prediction="results/plasmid_predictions-concatenated.tsv.gz",
        virus_prediction="results/virus_predections-concatenated.tsv.gz",
        assembly_database="results/assembly_database.csv.gz",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/make_assembly_database.txt",
    benchmark:
        "log/benchmark/make_assembly_database.txt"
    script:
        "../scripts/create_assembly_database.R"
