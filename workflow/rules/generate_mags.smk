### Generate Metaegenome-Assembled Genomes (MAGs)
## 1: bin contigs using different binners


rule bin_assemblies_semibin:
    input:
        masked_assembly="data/tmp/assembly/{sample}/assembly_ARG_masked.fasta",
        bam="data/tmp/assembly/{sample}/mapped_back/{sample}.bam",
    output:
        info="data/tmp/assembly/{sample}/semibin2/bins_info.tsv",
        contig_bins="data/tmp/assembly/{sample}/semibin2/contig_bins.tsv",
        bin_dir=directory("data/tmp/assembly/{sample}/semibin2/output_bins"),
    conda:
        "../envs/semibin.yaml"
    threads: config["semibin"]["threads"]
    log:
        "log/bin_assemblies_semibin/{sample}.txt",
    benchmark:
        "log/benchmark/bin_assemblies_semibin/{sample}.txt"
    shell:
        """
SemiBin2 single_easy_bin -t {threads}\
 -i {input.masked_assembly} -b {input.bam} -o $(dirname {output.info})\
  --sequencing-type=long_read --environment=global > {log} 2>&1
        """


rule bin_assemblies_metabat:
    input:
        masked_assembly="data/tmp/assembly/{sample}/assembly_ARG_masked.fasta",
        bam="data/tmp/assembly/{sample}/mapped_back/{sample}.bam",
    output:
        depth=temp("data/tmp/assembly/{sample}/metabat/{sample}-metabat_depth.txt"),
        info="data/tmp/assembly/{sample}/metabat/{sample}-metabat.BinInfo.txt",
    params:
        prefix=subpath(output.info, strip_suffix=".BinInfo.txt"),
    conda:
        "../envs/metabat2.yaml"
    threads: config["metabat2"]["threads"]
    log:
        "log/bin_assemblies_metabat/{sample}.txt",
    benchmark:
        "log/benchmark/bin_assemblies_metabat/{sample}.txt"
    shell:
        """
jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam} > {log} 2>&1

metabat2 -i {input.masked_assembly} -a {output.depth} -o {params.prefix} >> {log} 2>&1
        """


rule bin_assemblies_vamb:
    input:
        raw_assembly="data/tmp/assembly/{sample}/assembly.fasta",
        bam="data/tmp/assembly/{sample}/mapped_back/{sample}.bam",
    output:
        bins=directory("data/tmp/assembly/{sample}/vamb/bins"),
    params:
        prefix=subpath(output.bins, parent=True),
    conda:
        "../envs/vamb.yaml"
    threads: config["vamb"]["threads"]
    log:
        "log/bin_assemblies_vamb/{sample}.txt",
    benchmark:
        "log/benchmark/bin_assemblies_vamb/{sample}.txt"
    shell:
        """
rm -r {params.prefix}

vamb bin default --fasta {input.raw_assembly}\
 --bamdir $(dirname {input.bam}) --outdir {params.prefix}\
 -p {threads} --minfasta 200000 --seed 42 > {log} 2>&1
        """


## 2: classify bins taxonomically (GTDB-Tk)


rule classify_bins:
    input:
        bin_dir="data/tmp/assembly/{sample}/semibin2/output_bins",
    output:
        "data/tmp/assembly/{sample}/semibin2/gtdbtk/gtdbtk.bac120.summary.tsv",
    conda:
        "../envs/gtdbtk.yaml"
    threads: config["gtdbtk"]["threads"]
    log:
        "log/classify_bins/{sample}.txt",
    benchmark:
        "log/benchmark/classify_bins/{sample}.txt"
    shell:
        """
gtdbtk classify_wf --cpus {threads} --genome_dir {input.bin_dir}\
 --extension .fa.gz --out_dir $(dirname {output}) > {log} 2>&1
        """
