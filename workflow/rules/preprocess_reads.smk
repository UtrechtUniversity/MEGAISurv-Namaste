### Preprocess reads


rule read_quality_control:
    input:
        INPUT_DIR / "{sample}.fastq.gz",
    output:
        filtered="data/tmp/filtered/{sample}.fastq.gz",
        json="data/tmp/read_qc/{sample}.json",
        html="data/tmp/read_qc/{sample}.html",
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
        expand("data/tmp/read_qc/{sample}.json", sample=SAMPLES),
    output:
        expand("data/tmp/read_qc/summary/{sample}.json", sample=SAMPLES),
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/extract_read_qc_summaries.txt",
    benchmark:
        "log/benchmark/extract_read_qc_summaries.txt"
    shell:
        """
for inputfile in {input}
do
    bash scripts/extract_fastplong_json_summary.sh ${{inputfile}}
done > {log} 2>&1
        """


rule summarise_read_qc_data:
    input:
        expand("data/tmp/read_qc/summary/{sample}.json", sample=SAMPLES),
    output:
        "data/processed/read_qc_summary.csv",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/summarise_read_qc_data.txt",
    benchmark:
        "log/benchmark/summarise_read_qc_data.txt"
    script:
        "../scripts/collect_read_qc_data.R"
