## Screen antibiotic resistance mutations


rule screen_antibiotic_resistance_mutations:
    input:
        reads="results/filtered_reads/{sample}.fastq.gz",
    output:
        mpf_basics=expand(
            "results/resistance_mutations/{{sample}}/{file}",
            file=["accession", "class", "dna_accession", "dna_class"],
        ),
        mpf_summary=multiext(
            "results/resistance_mutations/{sample}/{sample}.",
            class_dna="class.dna.summary.txt",
            class_prot="class.prot.summary.txt",
            gene_dna="gene.dna.summary.txt",
            gene_prot="gene.prot.summary.txt",
        ),
        mpf_log=multiext(
            "results/resistance_mutations/{sample}/{sample}.", error="error", log="log"
        ),
        mpf_kma=multiext(
            "results/resistance_mutations/{sample}/{sample}.", "frag.gz", "res"
        ),
        mpf_input=multiext(
            "results/resistance_mutations/{sample}/{sample}.",
            "dna.input.tsv",
            "prot.input.tsv",
        ),
        mpf_table=multiext(
            "results/resistance_mutations/{sample}/{sample}.",
            dna="dna.updated_table_with_scores_and_mutations.tsv",
            prot="prot.updated_table_with_scores_and_mutations.tsv",
        ),
    params:
        db="resources/AMRFinder",
        output=subpath(output[0], parent=True),
    conda:
        "../envs/metapointfinder.yaml"
    threads: config["metapointfinder"]["threads"]
    log:
        "log/metapointfinder/{sample}.txt",
    benchmark:
        "log/benchmark/metapointfinder/{sample}.txt"
    shell:
        """
metapointfinder.py --force --input {input} --db {params.db}\
 --output {params.output} --identity 85 --threads {threads} > {log} 2>&1
        """


rule match_mutations_to_contigs:
    input:
        bam="results/assembly/{sample}/mapped_back/{sample}.bam",
        table="results/resistance_mutations/{sample}/{sample}.dna.updated_table_with_scores_and_mutations.tsv",
    output:
        "results/resistance_mutations/{sample}/matched_to_contigs.tsv",
    conda:
        "../envs/minimap2.yaml"
    threads: 1
    log:
        "log/match_mutations_to_contigs/{sample}.txt",
    benchmark:
        "log/benchmark/match_mutations_to_contigs/{sample}.txt"
    shell:
        """
printf "read\tflag\tcontig\tms_score\tAS_score\n" > {output} 2> {log}
samtools view -F 4 -F 256 -F 512 -F 2048\
 -N <(cut -f 3 {input.table}) {input.bam} | cut -f 1-3,13,14\
 >> {output} 2>> {log}
        """
