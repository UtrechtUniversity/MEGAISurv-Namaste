### Classify contigs
## taxonomic classification (Centrifuger + TaxonKit) *using stricter settings*


rule strict_taxonomic_classification:
    input:
        fasta="results/assembly/{sample}/assembly_ARG_masked.fasta",
        db=collect(
            "resources/centrifuger_db/cfr_hpv+gbsarscov2.{suffix}",
            suffix=["1.cfr", "2.cfr", "3.cfr"],
        ),
    output:
        tsv="results/taxonomic_classification/{sample}/centrifuger_masked-strict.tsv",
        quant="results/taxonomic_classification/{sample}/centrifuger_masked-strict-quant.tsv",
    params:
        db=subpath(input[1], strip_suffix=".1.cfr"),
    conda:
        "../envs/centrifuger.yaml"
    threads: config["centrifuger"]["threads"]
    log:
        "log/strict-centrifuger/{sample}.txt",
    benchmark:
        "log/benchmark/strict-centrifuger/{sample}.txt"
    shell:
        """
centrifuger -u {input.fasta} -t {threads} --min-hitlen 100\
 -x {params.db} > {output.tsv} 2> {log}

centrifuger-quant -x {params.db} -c {output.tsv} > {output.quant} 2> {log}
        """


rule strict_lookup_taxids:
    input:
        assembly="results/taxonomic_classification/{sample}/centrifuger_masked-strict.tsv",
        db=collect(
            "resources/taxdump/{dmp}",
            dmp=["names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp"],
        ),
    output:
        "results/taxonomic_classification/{sample}/centrifuger_masked-strict+taxa.tsv",
    threads: 1
    conda:
        "../envs/taxonkit.yaml"
    log:
        "log/strict-lookup_taxids/{sample}.txt",
    benchmark:
        "log/benchmark/strict-lookup_taxids/{sample}.txt"
    shell:
        """
taxonkit reformat {input.assembly} -I 3 --data-dir resources/taxdump\
 -f '{{s}}\t{{k}};{{p}};{{c}};{{o}};{{f}};{{g}};{{s}}' -F\
 > {output} 2> {log}
        """


rule strict_generate_microbiota_profiles:
    input:
        coverage_info="results/assembly/{sample}/mapped_back/{sample}-coverage.tsv",
        classifications="results/taxonomic_classification/{sample}/centrifuger_masked-strict+taxa.tsv",
    output:
        per_contig="results/microbiota_profile/{sample}-strict-per_contig.tsv",
        per_species="results/microbiota_profile/{sample}-strict-per_species.tsv",
    params:
        sample="{sample}",
    conda:
        "../envs/R_tidyverse.yaml"
    threads: 1
    log:
        "log/strict-generate_profiles/{sample}.txt",
    benchmark:
        "log/benchmark/strict-generate_profiles/{sample}.txt"
    script:
        "../scripts/generate_profiles.R"
