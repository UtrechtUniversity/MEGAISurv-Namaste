# Plasmid/virus/chromosome prediction

To analyse the potential of antibiotic resistance genes (ARGs) to spread to
other species, it helps to know whether they reside on the bacterial chromosome,
or on a plasmid or a virus-derived DNA element (mobile genetic elements).
To predict the nature of each sequence (contig), Namaste uses
[geNomad](https://portal.nersc.gov/genomad/) (version 1.8.0).

geNomad uses a combined marker gene and neural network (machine learning)
approach to deliver state-of-the-art predictions of mobile elements in
metagenomic data. The results are appended to the database with
ARGs, assembly statistics and taxonomic classifications to provide a
complete overview of the resistome and data-driven information on the
spreading potential of ARGs.

## Output files

geNomad comprises a complete workflow that generates many output files,
the most important of which are:

```txt
results/
  plasmid_prediction/
    {sample}/
      assembly_aggregated_classification.tsv    # Predictions with weighted
                                                # (marked-based and neural net)
                                                # classification scores
      assembly_plasmid_summary.tsv              # Summary of plasmid predictions
      assembly_virus_summary.tsv                # Summary of virus predictions
```

These are the files that are used to generate the final 'database' files.

For details, please see [output](output_files.md).
