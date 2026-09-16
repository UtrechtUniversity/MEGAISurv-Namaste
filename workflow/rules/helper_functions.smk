from pathlib import Path

# Check if input is text file or directory
INPUT_PATH = Path(config["input"])

if INPUT_PATH.is_dir() and not INPUT_PATH.samefile(""):
    # If directory, list all files with .fastq.gz extension
    INPUT_DIR = INPUT_PATH
    INPUT_FILES = list(INPUT_DIR.glob("*.fastq.gz"))
    SAMPLES = [file.stem.replace(".fastq", "") for file in INPUT_FILES]

elif INPUT_PATH.is_file():
    # If a file, try to read its content as input accessions
    # Read accession IDs from a text file (list with one accession per line)
    INPUT_DIR = Path("resources/public_metagenomes/")
    SAMPLES = []
    with open(INPUT_PATH, "r") as input_list:
        for line in input_list:
            SAMPLES.append(line.strip())

else:
    print(
        f"No valid input found in {INPUT_PATH}.\n"
        "This workflow requires the user to set an input file or directory.\n"
        "Please provide one in 'config/parameters.yaml'."
    )
    exit(1)


# Check if there are input samples
assert len(SAMPLES) > 0, (
    f"-----\nNo input samples found in {INPUT_PATH}.\n"
    "Please make sure that the input is either one of:"
    "1) a directory with gzipped FASTQ files (must have '.fastq.gz' extension)\n",
    "2) a text file with SRA accession numbers (one per line).\n-----\n",
)


rule download_raw_reads:
    output:
        temp("resources/public_metagenomes/{sample}.fastq.gz"),
    params:
        out_dir=subpath(output[0], parent=True),
    conda:
        "../envs/sracha.yaml"
    threads: config["download_raw_reads"]["threads"]
    log:
        "log/download_raw_reads/{sample}.txt",
    benchmark:
        "log/benchmark/download_raw_reads/{sample}.txt"
    shell:
        """
bash workflow/scripts/download_from_sra.sh -s {wildcards.sample}\
 -d {params.out_dir} -t {threads} > {log} 2>&1
        """


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
