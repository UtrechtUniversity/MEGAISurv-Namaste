# Read accession IDs from a text file (list with one accession per line)
INPUT_LIST = config["input_list"]
SAMPLES = []
with open(INPUT_LIST, "r") as input_list:
    for line in input_list:
        SAMPLES.append(line.strip())

assert len(SAMPLES) > 0, (
    f"-----\nNo input samples found in {INPUT_LIST}.\n"
    "Please make sure that there are SRA (sample) accession IDs in this file!\n"
    "-----\n"
)


rule download_raw_reads:
    output:
        temp("resources/public_metagenomes/{sample}.fastq.gz"),
    conda:
        "../envs/sracha.yaml"
    threads: config["download_raw_reads"]["threads"]
    log:
        "log/download_raw_reads/{sample}.txt",
    benchmark:
        "log/benchmark/download_raw_reads/{sample}.txt"
    script:
        "../scripts/download_from_sra.sh"


rule make_assembly_database:
    input:
        assembly_info=expand(
            "results/assembly/{sample}/assembly_info.txt", sample=SAMPLES
        ),
        mapped_coverage="results/contig_coverage.csv",
        arg_hits=expand("results/resistance_genes/{sample}.hmm.frag.gz", sample=SAMPLES),
        arg_results=expand("results/resistance_genes/{sample}.hmm.res", sample=SAMPLES),
        classification=expand(
            "results/taxonomic_classification/{sample}/centrifuger_masked+taxa.tsv",
            sample=SAMPLES,
        ),
        strict_classification=expand(
            "results/taxonomic_classification/{sample}/centrifuger_masked-strict+taxa.tsv",
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
        strict_classification="results/strict_classifications-concatenated.tsv.gz",
        genomad_scores="results/plasmid_prediction/aggregated_classification_scores-concatenated.tsv.gz",
        plasmid_prediction="results/plasmid_predictions-concatenated.tsv.gz",
        virus_prediction="results/virus_predictions-concatenated.tsv.gz",
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


rule make_mutation_database:
    input:
        arm_results=expand(
            "results/resistance_mutations/{sample}/{sample}.dna.updated_table_with_scores_and_mutations.tsv",
            sample=SAMPLES,
        ),
        arm_contigs=expand(
            "results/resistance_mutations/{sample}/matched_to_contigs.tsv",
            sample=SAMPLES,
        ),
        assembly_stats="results/assembly_stats-concatenated.tsv.gz",
        classification="results/classifications-concatenated.tsv.gz",
        strict_classification="results/strict_classifications-concatenated.tsv.gz",
        genomad_scores="results/plasmid_prediction/aggregated_classification_scores-concatenated.tsv.gz",
        plasmid="results/plasmid_predictions-concatenated.tsv.gz",
        virus="results/virus_predictions-concatenated.tsv.gz",
    output:
        mutation_database="results/mutation_database.csv.gz",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/make_mutation_database.txt",
    benchmark:
        "log/benchmark/make_mutation_database.txt"
    script:
        "../scripts/create_mutation_database.R"
