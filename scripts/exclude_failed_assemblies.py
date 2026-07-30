#!/usr/bin/env python3

"""
Create a list of samples that passed and failed assembly,
and move the reads of those that failed to a new directory
so that Snakemake will not fail while attempting to analyse
them further.
"""

from pathlib import Path
import pandas as pd
from yaml import load

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

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
    try:
        input_directory = Path(parameters_dict["input_directory"])
        input_files = list(input_directory.glob("*.fastq.gz"))
        samples = [file.stem.replace(".fastq", "") for file in input_files]

        samples_and_reads = {"Samples": samples, "Input_reads": input_files}

        return "read_dir", samples_and_reads

    except KeyError:
        # If not given a directory, it's a text file for downloading!
        input_list = Path(parameters_dict["input_list"])
        print("Found list of accession IDs instead of input files.")

        samples = []
        with open(input_list, "r") as infile:
            for line in infile:
                samples.append(line.strip())

        samples_dict = {"Samples": samples}

        return "accession_list", samples_dict


def find_assembly_files(samples=list):
    """
    Given a list of sample names, look for their corresponding assembly output
    files in the directory that Snakemake uses (default=results/assembly/)
    """
    file_list = []
    assembly_qc = []

    for sample in samples:
        assembly_file = Path(f"results/assembly/{sample}/assembly.fasta")

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
    for index in range(len(qc_dict["Samples"])):
        sample = qc_dict["Samples"][index]
        reads_file = qc_dict["Input_reads"][index]
        qc_verdict = qc_dict["Assembly_QC"][index]

        if qc_verdict == "exclude":
            target_dir = reads_file.parent / "cannot_assemble"
            target_dir.mkdir(parents=True, exist_ok=True)
            # method 'rename' actually moves files:
            reads_file.rename(target_dir / reads_file.name)
            print(
                f"Assembly of sample {sample} failed."
                "Moving input to subdirectory 'cannot_assemble'"
            )
        else:
            print(f"Assembly of sample {sample} passed.")
            continue


def update_accession_list(config_file=str, qc_dict=dict):
    """
    Given the input path to the assembly files and QC data (include/exclude),
    look up the list with input accessions, and update it so that
    it contains only samples that assembled correctly. Move the other
    accession IDs to a different file for backup.
    """
    parameters_dict = load(open(config_file, "r"), Loader=Loader)
    input_list = Path(parameters_dict["input_list"])
    updated_list = input_list.parent / (str(input_list.stem) + "-updated.txt")
    failed_list = input_list.parent / (str(input_list.stem) + "-failed_assembly.txt")

    with open(updated_list, "w") as updated:
        with open(failed_list, "w") as failed:
            for index in range(len(qc_dict["Samples"])):
                sample = qc_dict["Samples"][index]
                qc_verdict = qc_dict["Assembly_QC"][index]

                if qc_verdict == "include":
                    updated.write(f"{sample}\n")
                elif qc_verdict == "exclude":
                    failed.write(f"{sample}\n")
                else:
                    print(f"Found an unexpected QC result: {sample} - {qc_verdict}")
                    exit(1)

    # Move the old input list to "backup"
    input_list.rename(input_list.parent / (str(input_list.stem) + "-BACKUP.txt"))
    # And make the updated list the new 'current'
    updated_list.rename(input_list)


def write_assembly_qc_table(qc_dict=dict, outputfile=str):
    """
    Given the dictionary with sample names, input reads files,
    assembly files (if present) and assembly QC, write
    all that information to a tab-separated table file (TSV).
    """
    qc_df = pd.DataFrame(qc_dict)
    qc_df.to_csv(path_or_buf=outputfile, sep="\t", header=True, index=False)


def main():
    print("Looking for input reads by reading Snakemake parameters (YAML) file")
    method, samples_and_dict = parse_yaml(config_file="config/parameters.yaml")

    print("Found:\n", pd.DataFrame.from_dict(samples_and_dict), "\n")

    print("Looking up assemblies")
    assembly_dict = find_assembly_files(samples=samples_and_dict["Samples"])
    print("Found:\n", pd.DataFrame.from_dict(assembly_dict), "\n")

    combined_dict = {**samples_and_dict, **assembly_dict}

    if method == "read_dir":
        print("Moving reads that could not assemble")
        move_samples_without_assembly(qc_dict=combined_dict)

    elif method == "accession_list":
        print("Moving sample accessions that could not assemble from input list")
        update_accession_list(
            config_file="config/parameters.yaml", qc_dict=combined_dict
        )

    else:
        print("Error:")
        print(
            "Something went wrong in reading Snakemake's input"
            "from config/parameters.yaml"
        )
        exit(1)

    print("\nWriting summary report")
    write_assembly_qc_table(qc_dict=combined_dict, outputfile="results/assembly_qc.tsv")

    exit(0)


if __name__ == "__main__":
    exit(main())
