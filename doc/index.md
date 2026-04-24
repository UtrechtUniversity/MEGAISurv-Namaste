# Namaste :pray:

View project on [:simple-github: GitHub](https://github.com/UtrechtUniversity/MEGAISurv-Namaste) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [![Built with MkDocs](https://cdn.jsdelivr.net/npm/@intergrav/devins-badges@3/assets/cozy/built-with/mkdocs_vector.svg)](https://www.mkdocs.org/) &nbsp;+ &nbsp;[![Using Material for MkDocs](https://cdn.jsdelivr.net/gh/Andre601/devins-badges@v3.x-mkdocs-material/assets/compact-minimal/built-with/mkdocs-material_vector.svg)](https://squidfunk.github.io/mkdocs-material)

-----

## How to use :book:

To learn more about how to use the workflow, see the
[user manual](manual.md).

## What it does :simple-xyflow:

```mermaid
flowchart LR
    A[Raw reads] -->|Preprocess:  fastplong| B(High-quality reads)
    B -->|Assembly: metaFlye| C(Contigs)
    B -->|Screen ARMs: MetaPointFinder| M(Reads with resistance mutations)
    M -->|Map to: minimap2| C(Contigs)
    C -->|Screen ARGs: KMA/ResFinder| D(ARG-containing contigs)
    D -->|Mask: BEDtools| E(Masked contigs)
    E -->|Classify: Centrifuger| F(Taxonomy-assigned contigs)
    C -->|Predict plasmids: geNomad| G(Plasmid/virus/chromosome contigs)
```

The workflow works as follows:

1. [Preprocess reads](preprocess.md)

2. [Assemble metagenomes](assembly.md)

3. [Screen for presence of antibiotic resistance genes](arg_screening.md)

4. [Also screen for antibiotic resistance mutations](arm_screening.md)

5. [Classify contigs taxonomically](taxonomy.md)

6. [Predict whether contigs derive from chromosomes, plasmids or viruses](genomad.md)

And comes with several [additional functions](extra_functions.md), like
automated downloading of reference databases and combining resistome data
in comprehensive CSV-formatted tables.

## What it creates :simple-files:

For an overview of the output files that are created in the process, please
consult [output files](output_files.md).
