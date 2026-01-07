### Classify contigs
## 1: taxonomic classification (Centrifuger + TaxonKit)


rule taxonomic_classification:
    input:
        fasta="data/tmp/assembly/{sample}/assembly_ARG_masked.fasta",
        db=config["centrifuger"]["database"] + ".2.cfr",
    output:
        tsv="data/tmp/centrifuger/{sample}/centrifuger_masked.tsv",
        quant="data/tmp/centrifuger/{sample}/centrifuger_masked-quant.tsv",
    params:
        db=subpath(input.db, strip_suffix=".2.cfr"),
    conda:
        "../envs/centrifuger.yaml"
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
        "data/tmp/centrifuger/{sample}/centrifuger_masked.tsv",
    output:
        "data/tmp/centrifuger/{sample}/centrifuger_masked+taxa.tsv",
    params:
        taxondb=config["taxon_database"],
    threads: 1
    conda:
        "../envs/taxonkit.yaml"
    log:
        "log/lookup_taxids/{sample}.txt",
    benchmark:
        "log/benchmark/lookup_taxids/{sample}.txt"
    shell:
        """
taxonkit reformat {input} -I 3 --data-dir {params.taxondb}\
 -f '{{s}}\t{{k}};{{p}};{{c}};{{o}};{{f}};{{g}};{{s}}' -F\
 > {output} 2> {log}
        """


rule generate_microbiota_profiles:
    input:
        assembly_info="data/tmp/assembly/{sample}/assembly_info.txt",
        classifications="data/tmp/centrifuger/{sample}/centrifuger_masked+taxa.tsv",
    output:
        per_contig="data/tmp/microbiota_profiles/{sample}-per_contig.tsv",
        per_species="data/tmp/microbiota_profiles/{sample}-per_species.tsv",
    params:
        sample="{sample}",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/generate_profiles/{sample}.txt",
    benchmark:
        "log/benchmark/generate_profiles/{sample}.txt"
    script:
        "../scripts/generate_profiles.R"


## 2: chromosome/plasmid/virus (geNomad)


rule genomad:
    input:
        fasta="data/tmp/assembly/{sample}/assembly.fasta",
        db=config["genomad"]["database"],
    output:
        aggregated_classification="data/tmp/genomad/{sample}/assembly_aggregated_classification/assembly_aggregated_classification.tsv",
        plasmid_summary="data/tmp/genomad/{sample}/assembly_summary/assembly_plasmid_summary.tsv",
        virus_summary="data/tmp/genomad/{sample}/assembly_summary/assembly_virus_summary.tsv",
    params:
        output_dir=subpath(output.plasmid_summary, ancestor=2),
    conda:
        "../envs/genomad.yaml"
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
