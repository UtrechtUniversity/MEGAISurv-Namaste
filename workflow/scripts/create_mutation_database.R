#!/usr/bin/env Rscript

suppressPackageStartupMessages(
  library(tidyverse)
)

sink(
  file = file(snakemake@log[[1]], open = "wt"),
  type = "message"
)

################################################################################
# Pseudocode/work description:
#
# 1.a for each read, lookup in the file
#     `results/resistance_mutations/{sample}/{sample}.dna.updated_table_with_scores_and_mutations.tsv`
#    - if they contain the WT or R genotype (column 10 = 'Status' ('Wildtype' or 'Resistant'))
#    - what gene they mapped to (column 2 = 'gene')
#    - mutations (column 8 = 'DetectedMutations')
#    - antibiotic class (column 1 = 'class')
#    - write to separate file,
# 2.b calculate sums per gene, per contig: x WT reads, y R reads
#    - attach gene and mutations
#
# 2. match 1.b. with other 'assembly db' data:
#    - taxonomy, plasmid prediction, ... basically re-use that script,
#      or better yet, use the prepared concatenated files!
################################################################################


arm_results_files <- snakemake@input[["arm_results"]] # MetaPointFinder's '.dna.updated_table_with_scores_and_mutations.tsv' files
arm_to_contig_files <- snakemake@input[["arm_contigs"]] # Match reads back to contigs using minimap2's BAM files + samtools view

read_stats <- function(filename, name_position) {
  # Cut the sample name from the file path using the name's position in the
  # absolute file path (counting from the end: 1 = last, 2 = second last, etc.)
  sample_name <- str_split_1(string = filename, pattern = "/") %>%
    tail(name_position) %>%
    head(1)

  df <- read_delim(
    file = filename,
    show_col_types = FALSE
  ) %>%
    mutate(sample = sample_name)

  return(df)
}

# Read MetaPointFinder's updated DNA tables
arm_table <- do.call(
  rbind,
  lapply(X = arm_results_files, FUN = read_stats, name_position = 2)
)

# Read the read-to-contig mapping files
arm_mapping <- do.call(
  rbind,
  lapply(X = arm_to_contig_files, FUN = read_stats, name_position = 2)
)

# Combine the two
mapped_arm_df <- left_join(
  x = arm_table,
  y = arm_mapping,
  by = c("sample", "read")
)

# Summarise the number of reads with wildtype and resistant genotype
columns <- c(Wildtype = 0, Resistant = 0)

mutations <- mapped_arm_df %>%
  group_by(sample, contig, class, gene) %>%
  mutate(mutations = paste0(DetectedMutations, collapse = ",")) %>%
  group_by(sample, contig, class, gene, mutations) %>%
  count(Status) %>%
  pivot_wider(names_from = Status, values_from = n) %>%
  add_column(., !!!columns[setdiff(names(columns), colnames(.))]) %>%
  mutate(fraction_resistant = Resistant / Wildtype) %>%
  mutate(fraction_resistant = if_else(!is.na(fraction_resistant),
    fraction_resistant,
    0
  ))


# Now append other information, starting with assembly statistics
assembly_stats <- read_delim(
  file = snakemake@input[["assembly_stats"]],
  delim = "\t", show_col_types = FALSE
)


assembly_stats_summary <- assembly_stats %>%
  group_by(sample) %>%
  summarise(
    mean_assembly_depth_flye = mean(contig_depth_flye),
    mean_assembly_depth_samtools = mean(contig_depth_samtools),
    total_assembly_length = sum(contig_length)
  )

arm_and_assembly <- left_join(
  x = mutations,
  y = assembly_stats_summary,
  by = "sample"
) %>%
  left_join(
    x = .,
    y = assembly_stats,
    by = c("sample", "contig")
  )


# Then taxonomic classifications
classifications <- read_delim(
  file = snakemake@input[["classification"]],
  delim = "\t", show_col_types = FALSE
)

strict_classifications <- read_delim(
  file = snakemake@input[["strict_classification"]],
  delim = "\t", show_col_types = FALSE
)


db_with_classifications <- left_join(
  x = arm_and_assembly,
  y = full_join(
    x = classifications %>%
      select(sample, readID, Species, Taxonomy),
    y = strict_classifications %>%
      select(sample, readID, Species, Taxonomy),
    by = c("sample", "readID"),
    suffix = c("", "_strict")
  ),
  by = c("sample", "contig" = "readID")
)


# And finally geNomad scores and plasmid/virus predictions
genomad_scores <- read_delim(
  file = snakemake@input[["genomad_scores"]],
  delim = "\t", show_col_types = FALSE
)

plasmid_predictions <- read_delim(
  file = snakemake@input[["plasmid"]],
  delim = "\t", show_col_types = FALSE
)

virus_predictions <- read_delim(
  file = snakemake@input[["virus"]],
  delim = "\t", show_col_types = FALSE
)


db_with_mge <- left_join(
  x = db_with_classifications,
  y = genomad_scores,
  by = c("sample", "contig" = "seq_name")
) %>%
  left_join(
    x = .,
    y = plasmid_predictions %>%
      select(
        sample, contig, plasmid_topology,
        plasmid_genes, conjugation_genes,
        amr_genes
      ),
    by = c("sample", "contig")
  ) %>%
  left_join(
    x = .,
    y = virus_predictions %>%
      select(
        sample, contig, virus_topology,
        virus_genes, virus_taxonomy
      ),
    by = c("sample", "contig")
  ) %>%
  mutate(genomad_prediction = case_when(
    !is.na(plasmid_topology) ~ "plasmid",
    !is.na(virus_topology) ~ "virus",
    TRUE ~ "chromosome"
  ))


# Write the final output to a (gzipped) CSV file
write_csv(
  x = db_with_mge %>% distinct(),
  file = snakemake@output[["mutation_database"]],
)
