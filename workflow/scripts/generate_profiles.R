#! /usr/bin/env Rscript

# For samples provided by Snakemake, read the Flye assembly statistics
# (`*_info.txt` file) and the Centrifuger output with
# (TaxonKit-)assigned taxonomy names to create tables that list for each
# sample the normalised/quantified taxa.

library(tidyverse)

# Prepare the output directory if it does not exist yet:
output_dir <- snakemake@output["per_contig"] %>%
  as.character() %>%
  dirname()
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read the 'sample' variable from Snakemake
sample <- snakemake@params[["sample"]]

# Read the Flye assembly statistics file
info_df <- read_delim(snakemake@input[["assembly_info"]],
  delim = "\t", show_col_types = FALSE
)
# And the taxonomic classifications from Centrifuger,
# interpreted by TaxonKit:
classification_df <- read_delim(snakemake@input[["classifications"]],
  delim = "\t", col_names = FALSE, skip = 1,
  show_col_types = FALSE
)

# Manually adjust the column names (to simplify further scripting)
colnames(info_df) <- c("seq_name", "length", "depth", "circular", "repeat", "multiples", "alt_group", "graph_path")
colnames(classification_df) <- c("readID", "seqID", "taxID", "score", "2ndBestScore", "hitLength", "queryLength", "numMatches", "Species", "Lineage")

# Join (merge) the two together
quantified_df <- left_join(
  x = classification_df,
  y = info_df %>%
    select(seq_name, length, depth),
  by = c("readID" = "seq_name")
)

# Calculate the total bases per contig (length * depth)
quantified_df <- quantified_df %>%
  mutate(total_bases = length * depth)

# Calculate the grand total per sample (sum up all basepairs)
quantified_df$grand_total <- colSums(quantified_df %>% select(total_bases))
# And also add the sample names to the dataframe (also to simplify downstream processing)
quantified_df$Sample <- sample

# Calculate percentage of total bases per contig
quantified_df <- quantified_df %>%
  mutate(Percentage = total_bases / grand_total * 100) %>%
  arrange(desc(Percentage))

# And save as a TSV file
write_delim(x = quantified_df, file = as.character(snakemake@output["per_contig"]), delim = "\t")

# Calculate percentages per species
quantified_per_species <- quantified_df %>%
  group_by(Sample, Species, Lineage) %>%
  summarise(
    total_bases = sum(total_bases),
    Percentage = sum(Percentage)
  ) %>%
  arrange(desc(Percentage))

# And also save as a TSV file
write_delim(x = quantified_per_species, file = as.character(snakemake@output["per_species"]), delim = "\t")
