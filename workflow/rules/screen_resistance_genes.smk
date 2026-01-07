### Screen antibiotic resistance genes


rule create_resfinder_database:
    output:
        multiext(
            "resources/resfinder_db/all",
            out1=".comp.b",
            out2=".fsa",
            out3=".length.b",
            out4=".seq.b",
        ),
    conda:
        "../envs/kma.yaml"
    threads: 1
    log:
        "log/create_resfinder_database.txt",
    benchmark:
        "log/benchmark/create_resfinder_database.txt"
    shell:
        """
rm -rf $(dirname {output.out1})
bash workflow/scripts/prepare_resfinder.sh > {log} 2>&1

(cd $(dirname {output.out1}) && python3 INSTALL.py) >> {log} 2>&1
        """


rule screen_antibiotic_resistance_genes:
    input:
        db=expand(
            "resources/resfinder_db/all" + ".{extension}",
            extension=["comp.b", "fsa", "length.b", "seq.b"],
        ),
        assembly="data/tmp/assembly/{sample}/assembly.fasta",
    output:
        aln="data/tmp/kma/{sample}.hmm.aln",
        frag="data/tmp/kma/{sample}.hmm.frag.gz",
        fsa="data/tmp/kma/{sample}.hmm.fsa",
        res="data/tmp/kma/{sample}.hmm.res",
    params:
        db="resources/resfinder_db/all",
        prefix=subpath(output.aln, strip_suffix=".aln"),
    conda:
        "../envs/kma.yaml"
    threads: config["kma"]["threads"]
    log:
        "log/kma/{sample}.txt",
    benchmark:
        "log/benchmark/kma/{sample}.txt"
    shell:
        """
kma -t {threads} -bcNano -t_db {params.db} -i {input.assembly} -o {params.prefix}\
 -ont -hmm > {log} 2>&1
        """


rule mask_resistance_gene_positions:
    input:
        frag="data/tmp/kma/{sample}.hmm.frag.gz",
        assembly="data/tmp/assembly/{sample}/assembly.fasta",
    output:
        gene_locations="data/tmp/kma/{sample}.locations.txt",
        masked_assembly="data/tmp/assembly/{sample}/assembly_ARG_masked.fasta",
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    log:
        "log/kma/{sample}-positions.txt",
    benchmark:
        "log/benchmark/kma/{sample}-positions.txt"
    shell:
        """
# Cut contig ID, start and stop positions from file
zcat {input.frag} | cut -f 7-9 | sort | uniq > {output.gene_locations}

maskFastaFromBed -fi {input.assembly} -bed {output.gene_locations}\
 -fo {output.masked_assembly} > {log} 2>&1
        """
