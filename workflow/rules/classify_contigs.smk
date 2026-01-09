### Classify contigs
## 1: taxonomic classification (Centrifuger + TaxonKit)


rule download_centrifuger_database:
    output:
        multiext(
            "resources/centrifuger_db/cfr_hpv+gbsarscov2.",
            out1="1.cfr",
            out2="2.cfr",
            out3="3.cfr",
        ),
    conda:
        "../envs/bash.yaml"
    threads: 1
    log:
        "log/download_centrifuger_database.txt",
    benchmark:
        "log/benchmark/download_centrifuger_database.txt"
    shell:
        """
bash workflow/scripts/download_centrifuger_db.sh > {log} 2>&1
        """


rule taxonomic_classification:
    input:
        fasta="data/tmp/assembly/{sample}/assembly_ARG_masked.fasta",
        db=collect(
            "resources/centrifuger_db/cfr_hpv+gbsarscov2.{suffix}",
            suffix=["1.cfr", "2.cfr", "3.cfr"],
        ),
    output:
        tsv="data/tmp/centrifuger/{sample}/centrifuger_masked.tsv",
        quant="data/tmp/centrifuger/{sample}/centrifuger_masked-quant.tsv",
    params:
        db="resources/centrifuger_db/cfr_hpv+gbsarscov2"
    conda:
        "../envs/centrifuger.yaml"
    threads: config["centrifuger"]["threads"]
    log:
        "log/centrifuger/{sample}.txt",
    benchmark:
        "log/benchmark/centrifuger/{sample}.txt"
    shell:
        """
centrifuger -u {input.fasta} -t {threads}\
 -x {params.db} > {output.tsv} 2> {log}

centrifuger-quant -x {params.db} -c {output.tsv} > {output.quant} 2> {log}
        """


rule download_taxdump:
    output:
        multiext(
            "resources/taxdump/",
            names="names.dmp",
            nodes="nodes.dmp",
            delnodes="delnodes.dmp",
            merged="merged.dmp",
        ),
    threads: 1
    conda:
        "../envs/bash.yaml"
    log:
        "log/download_taxdump.txt",
    benchmark:
        "log/benchmark/download_taxdump.txt"
    shell:
        """
bash workflow/scripts/download_taxdump.sh > {log} 2>&1
        """


rule lookup_taxids:
    input:
        assembly="data/tmp/centrifuger/{sample}/centrifuger_masked.tsv",
        db=collect(
            "resources/taxdump/{dmp}",
            dmp=["names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp"],
        ),
    output:
        "data/tmp/centrifuger/{sample}/centrifuger_masked+taxa.tsv",
    threads: 1
    conda:
        "../envs/taxonkit.yaml"
    log:
        "log/lookup_taxids/{sample}.txt",
    benchmark:
        "log/benchmark/lookup_taxids/{sample}.txt"
    shell:
        """
taxonkit reformat {input.assembly} -I 3 --data-dir resources/taxdump\
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


rule download_genomad_database:
    output:
        db=directory("resources/genomad_db"),
    conda:
        "../envs/genomad.yaml"
    threads: 1
    log:
        "log/download_genomad_database.txt",
    benchmark:
        "log/benchmark/download_genomad_database.txt"
    shell:
        """
mkdir -p $(dirname {output.db})
genomad download-database $(dirname {output.db}) > {log} 2>&1
        """


rule genomad:
    input:
        fasta="data/tmp/assembly/{sample}/assembly.fasta",
        db="resources/genomad_db",
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
