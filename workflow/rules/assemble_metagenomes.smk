### Assemble metagenomes


rule metagenomic_assembly:
    input:
        "results/filtered_reads/{sample}.fastq.gz",
    output:
        assembly="results/assembly/{sample}/assembly.fasta",
        info="results/assembly/{sample}/assembly_info.txt",
    params:
        output_dir=subpath(output.assembly, parent=True),
        settings="--meta",
    conda:
        "../envs/flye.yaml"
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
        "results/assembly/{sample}/assembly.fasta",
    output:
        report="results/contig_qc/{sample}/report.html",
        icarus="results/contig_qc/{sample}/icarus.html",
        log="results/contig_qc/{sample}/metaquast.log",
    params:
        output_dir=subpath(output.report, parent=True),
    conda:
        "../envs/quast.yaml"
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
        expand("results/assembly/{sample}/assembly.fasta", sample=SAMPLES),
    output:
        "results/assembly/assembly_statistics-seqkit.tsv",
    conda:
        "../envs/seqkit.yaml"
    threads: config["simple_stats"]["threads"]
    log:
        "log/simple_assembly_statistics.txt",
    benchmark:
        "log/benchmark/simple_assembly_statistics.txt"
    shell:
        """
seqkit stats -Ta -j {threads} > {output} 2> {log}
        """


rule map_reads_to_assembly:
    input:
        reads="results/filtered_reads/{sample}.fastq.gz",
        assembly="results/assembly/{sample}/assembly.fasta",
    output:
        bam="results/assembly/{sample}/mapped_back/{sample}.bam",
        cov="results/assembly/{sample}/mapped_back/{sample}-coverage.tsv",
    conda:
        "../envs/minimap2.yaml"
    threads: config["minimap2"]["threads"]
    log:
        "log/map_reads_to_assembly/{sample}.txt",
    benchmark:
        "log/benchmark/map_reads_to_assembly/{sample}.txt"
    shell:
        """
minimap2 -x map-ont -t {threads} -a {input.assembly} {input.reads} 2> {log} |\
 samtools sort -o {output.bam} --write-index - > {log} 2>&1
samtools coverage -o {output.cov} {output.bam} >> {log} 2>&1
        """


rule summarise_read_mapping:
    input:
        cov_files=expand(
            "results/assembly/{sample}/mapped_back/{sample}-coverage.tsv",
            sample=SAMPLES,
        ),
    output:
        "results/contig_coverage.csv",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/summarise_read_mapping.txt",
    benchmark:
        "log/benchmark/summarise_read_mapping.txt"
    script:
        "../scripts/collect_read_mapping_data.R"
