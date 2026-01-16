#!/usr/bin/env python3

"""
Create a list of samples that passed and failed assembly,
and move the reads of those that failed to a new directory
so that Snakemake will not fail while attempting to analyse
them further.
"""

from pathlib import Path
import pandas as pd
from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

# Read Snakemake's configuration file as YAML, using PyYAML following
# its documentation: https://pyyaml.org/wiki/PyYAMLDocumentation


def parse_yaml(config_file=str):
    """
    Read Snakemake's config/parameters.yaml file to extract the input files
    and their sample names.
    Also extract the output file path.

    Return: 1) a dictionary with sample names and reads files,
            2) a Path object of the output directory
    """
    parameters_dict = load(open(config_file, "r"), Loader=Loader)
    input_directory = Path(parameters_dict["input_directory"])
    output_directory = Path(parameters_dict["output_directory"])
    input_files = list(input_directory.glob("*.fastq.gz"))
    samples = [file.stem.replace(".fastq", "") for file in input_files]

    samples_and_reads = {"Samples": samples, "Input_reads": input_files}

    return samples_and_reads, output_directory


def find_assembly_files(samples=list, output_directory=str):
    """
    Given a list of sample names, look for their corresponding assembly output
    files in the directory that Snakemake uses (default=data/tmp/assembly)
    """
    file_list = []
    assembly_qc = []

    for sample in samples:
        assembly_file = output_directory / f"assembly/{sample}/assembly.fasta"

        # If the assembly file exists
        if assembly_file.is_file():
            # Record its name and include it in further analyses,
            file_list.append(assembly_file)
            assembly_qc.append("include")
        else:
            # Otherwise, record 'NA' and exclude it
            file_list.append("NA")
            assembly_qc.append("exclude")

    assembly_dict = {"Assembly_file": file_list, "Assembly_QC": assembly_qc}

    return assembly_dict


def move_samples_without_assembly(qc_dict=dict):
    """
    Given the input path to the reads files and a dictionary with
    assembly files and QC data (include/exclude), move the reads
    from samples that did not assemble (= exclude) to a new
    subdirectory 'cannot_assemble'.
    """
    # Path.mkdir(input_path / "cannot_assemble", parents=True)

    for index in range(len(qc_dict["Samples"])):
        sample = qc_dict["Samples"][index]
        reads_file = qc_dict["Input_reads"][index]
        qc_verdict = qc_dict["Assembly_QC"][index]

        if qc_verdict == "exclude":
            target_dir = reads_file.parent / "cannot_assemble"
            Path.mkdir(target_dir, parents=True, exist_ok=True)
            # method 'rename' actually moves files:
            reads_file.rename(target_dir / reads_file.name)
            print(
                f"Assembly of sample {sample} failed. Moving input to subdirectory 'cannot_assemble'"
            )
        else:
            print(f"Assembly of sample {sample} passed.")
            continue


def write_assembly_qc_table(qc_dict=dict, outputfile=str):
    """
    Given the dictionary with sample names, input reads files,
    assembly files (if present) and assembly QC, write
    all that information to a tab-separated table file (TSV).
    """
    qc_df = pd.DataFrame(qc_dict)
    qc_df.to_csv(path_or_buf=outputfile, sep="\t", header=True, index=False)


def main():
    samples_and_reads_dict, output_directory = parse_yaml(
        config_file="config/parameters.yaml"
    )

    assembly_dict = find_assembly_files(
        samples=samples_and_reads_dict["Samples"], output_directory=output_directory
    )

    combined_dict = {**samples_and_reads_dict, **assembly_dict}

    move_samples_without_assembly(qc_dict=combined_dict)

    write_assembly_qc_table(
        qc_dict=combined_dict, outputfile="data/tmp/assembly_qc.tsv"
    )

    exit(0)


if __name__ == "__main__":
    exit(main())
