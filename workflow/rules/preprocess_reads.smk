### Preprocess reads


rule read_quality_control:
    input:
        "resources/public_metagenomes/{sample}.fastq.gz",
    output:
        filtered=temp("results/filtered_reads/{sample}.fastq.gz"),
        json="results/read_qc/{sample}.json",
        html="results/read_qc/{sample}.html",
    wildcard_constraints:
        sample="SAM[A-Z]+[0-9]+",
    conda:
        "../envs/fastplong.yaml"
    threads: config["fastplong"]["threads"]
    log:
        "log/read_qc/{sample}.txt",
    benchmark:
        "log/benchmark/read_qc/{sample}.txt"
    shell:
        """
fastplong --in {input} --out {output.filtered}\
 --trim_poly_x --trimming_extension 10\
 --json {output.json} --html {output.html}\
 --thread {threads} > {log} 2>&1
        """


rule extract_read_qc_summaries:
    input:
        rules.read_quality_control.output.json,
    output:
        "results/read_qc/summary/{sample}.json",
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/extract_read_qc_summaries/{sample}.txt",
    benchmark:
        "log/benchmark/extract_read_qc_summaries/{sample}.txt"
    script:
        "../scripts/extract_fastplong_json_summary.sh"


rule summarise_read_qc_data:
    input:
        json=expand("results/read_qc/summary/{sample}.json", sample=SAMPLES),
    output:
        "results/read_qc_summary.csv",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/summarise_read_qc_data.txt",
    benchmark:
        "log/benchmark/summarise_read_qc_data.txt"
    script:
        "../scripts/collect_read_qc_data.R"


rule subsample_large_samples:
    input:
        "results/filtered_reads/{sample}.fastq.gz",
    output:
        "results/downsampled_reads/{sample}.fastq.gz",
    conda:
        "../envs/seqkit.yaml"
    threads: config["subsample_samples"]["threads"]
    resources:
        mem_mb=int(20000),
        runtime=int(30),
    log:
        "log/subsample_samples/{sample}.txt",
    benchmark:
        "log/benchmark/subsample_samples/{sample}.txt"
    shell:
        """
seqkit sample2 -2 -n 1000000 -o {output} --threads {threads} {input}
        """
