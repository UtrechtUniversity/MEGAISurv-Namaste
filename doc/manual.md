# Namaste user manual

## Quick start :material-run-fast:

Install dependencies:
[git](https://git-scm.com/downloads/),
[mamba](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install)
and
[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation).

Download the repository:

``` bash
git clone https://github.com/UtrechtUniversity/MEGAISurv-Namaste.git
```

Move into the downloaded directory:

```bash
cd MEGAISurv-Namaste
```

Collect long-read metagenomes for input, for example:
(using [`sracha`](https://rnabioco.github.io/sracha-rs/)
to download public metagenomes from the European Nucleotide Archive
([ENA](https://www.ebi.ac.uk/ena/browser/home)).)

```bash
sracha get --output-dir data SRR28879900
sracha get --output-dir data SRR28879905
sracha get --output-dir data SRR28879907
```

Or insert a link to files in a different location:

```bash
cd data
ln -s /path/to/metagenomes/*fastq.gz .
```

(Where you replace "/path/to/metagenomes/" with the actual path on your system!)

Input files in the directory `data/` should be automatically recognised.
Test this by doing a dry-run:

```bash
snakemake --profile config -n
```

If that returns no errors, proceed with running the actual workflow:

```bash
snakemake --profile config
```

## 1. Before you start :octicons-info-16:

The Namaste workflow processes long-read metagenomes from the specified
input folder (default='data/'). It is based on the
[Snakemake](https://snakemake.readthedocs.io/en/stable/)
workflow management system and uses
[mamba](https://mamba.readthedocs.io/en/latest/index.html)
for installing dependencies.
Furthermore, you will need [git](https://git-scm.com/) as that is currently
the only available method to install Namaste. (Download from GitHub.)

Input files are detected automatically as long as they are in the specified
input folder, which is defined in
[`config/parameters.yaml`](https://github.com/UtrechtUniversity/MEGAISurv-Namaste/blob/main/config/parameters.yaml).

### Estimated disk use :material-harddisk:

Besides the input metagenomes that may be big, Namaste needs a number of
databases to work. These include:

- AMRFinder: ~10MB
- Centrifuger 'cfr_hpv+gbsarscov2': 43GB
- geNomad default database: 1.4GB
- MetaPointFinder (=AMRFinder)
- ResFinder: ~100MB
- NCBI Blast taxonomy (taxdump): ~500MB

Total: ~45GB

### Download and install software :material-download:

Before you begin, you need to install:
(follow these links to find installation instructions)

1. [git](https://git-scm.com/downloads/)

2. [mamba](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install)

3. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation)

We recommend Snakemake is installed via mamba. This is also the default
and linked above. Namaste has been tested with Snakemake version 9.3.0
and is expected to work with any version >=6.

When you have these tools installed, you can download Namaste:

``` bash
git clone https://github.com/UtrechtUniversity/MEGAISurv-Namaste.git
```

Change directory into the newly downloaded folder to get started:

```bash
cd MEGAISurv-Namaste
```

You may rename this folder if you want to, for example:

```bash
mv MEGAISurv-Namaste namaste
cd namaste
```

### Adjusting parameters :material-tune:

Namaste has a few options that may be modified by the user.
These are listed in two configuration files:

```txt
config
├── config.yaml
└── parameters.yaml
```

The most important are the input directory
('input_directory' in `config/parameters.yaml`, default='data/')
and the number of CPU threads to use
('jobs' in `config/config.yaml`, default=72).

Please modify these in your favourite text editor to fit your setup.

## 2. Running the workflow :material-run:

The workflow is fully automated and should complete with one command.
For details on what happens under the hood, see the tab 'Workflow details'.

One can do a 'dry-run' to test if all preparations have been satisfied:

``` bash
snakemake --profile config -n
```

To run the actual workflow:

``` bash
snakemake --profile config
```

### Exceptions: failed assembly :warning:

Sometimes, a metagenome may not contain sufficient reads to generate a
de novo assembly. For example, when negative controls (blanks) are included.
The workflow cannot successfully complete the analysis of these samples and
returns errors for steps downstream of the assembly. There is a script
included to automatically flag these samples and move them to a subdirectory,
so that the workflow may exit successfully.

```bash
python scripts/exclude_failed_assemblies.py
```

This script reads the `config/parameters.yaml` file to determine the correct
input directory. Input reads are moved to a subdirectory `cannot_assemble`.
It also generates a simple QC report listing which samples did and which
did not yield a working assembly (fasta) file.

## 3. Interpreting results :material-magnify:

After running the workflow, the user is presented with a number of output files.
These are described in detail under the tab '[Output files](output_files.md)'.
