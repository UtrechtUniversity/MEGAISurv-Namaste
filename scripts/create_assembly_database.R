#!/usr/bin/env Rscript

# Look up information from KMA, Flye, Resfinder, geNomad, and Centrifuger
# to create a comprehensive overview table of all the ARG-containing contigs.

library(tidyverse)
library(here)

assembly_stats_files <- Sys.glob(paths = here("data", "tmp", "assembly", "*", "assembly_info.txt"))
arg_hits_files <- Sys.glob(paths = here("data", "tmp", "kma", "*.hmm.frag.gz"))
classification_files <- Sys.glob(paths = here("data", "tmp", "centrifuger", "*", "centrifuger_masked+taxa.tsv"))
genomad_scores_files <- Sys.glob(paths = here("data", "tmp", "genomad", "*", "assembly_aggregated_classification", "assembly_aggregated_classification.tsv"))
genomad_plasmid_files <- Sys.glob(paths = here("data", "tmp", "genomad", "*", "assembly_summary", "assembly_plasmid_summary.tsv"))
genomad_virus_files <- Sys.glob(paths = here("data", "tmp", "genomad", "*", "assembly_summary", "assembly_virus_summary.tsv"))

duplicates <- read_delim(here("data", "samples_with_multiple_runs.tsv"))


## 1. Resistance gene information (gene name, contig position, hit length?)

read_arg_stats <- function(filename){
  # Cut the sample name from the file path
  sample_name <- basename(filename) %>% gsub(
    pattern = ".hmm.frag.gz",
    replacement = "",
    x = .
  )
  
  df <- read_delim(
    file = filename,
    col_names = FALSE,
    show_col_types = FALSE
  ) %>%
    mutate(sample = sample_name)
  
  return(df)
}

# Read antibiotic resistance gene information from KMA output files
arg_stats <- do.call(
  rbind,
  lapply(X = arg_hits_files, FUN = read_arg_stats)
) %>%
  select(sample, X7, X6, X8, X9, X4, X5)

# Make the dataframe nicer by adding column names, and add columns 'hit_length'
#  and 'arg_short'
colnames(arg_stats) <- c("sample", "contig", "arg", "start_position", "end_position", "hit_start", "hit_end")
arg_stats <- arg_stats %>%
  mutate(hit_length = hit_end - hit_start + 1,
         arg_short = gsub(pattern = "_.*",
                          replacement = "",
                          x = arg)) %>%
  select(-c("hit_start", "hit_end"))

# Deduplicate samples by removing runs that derive from the same biosample
deduplicated_args <- arg_stats %>%
  filter(! sample %in% duplicates$run)
# (Remember that the the dataframe contained both entries for runs and the
#  samples: by removing runs we keep only the samples, which are the two
#  corresponding runs concatenated.)

## Look up antibiotic classes from ResFinder's `phenotypes.txt` file
arg_annotation <- read_delim(
  file = "https://bitbucket.org/genomicepidemiology/resfinder_db/raw/cf9bbc7b13f04de987f7dd4a3a1440c7af0b1ce0/phenotypes.txt",
  show_col_types = FALSE,
  delim = "\t"
) %>%
  rename("Resistance_mechanism" = "Mechanism of resistance")

deduplicated_args <- left_join(
  x = deduplicated_args,
  y = arg_annotation %>%
    select(`Gene_accession no.`, Class, Phenotype, Resistance_mechanism),
  by = c("arg" = "Gene_accession no.")
)


## 2. Assembly stats (length, depth and circularity)

read_stats <- function(filename, name_position){
  # Cut the sample name from the file path using the name's position in the
  # absolute file path (counting from the end: 1 = last, 2 = second last, etc.)
  sample_name <- str_split_1(string = filename, pattern = "/") %>% tail(name_position) %>% head(1)
  
  df <- read_delim(
    file = filename,
    show_col_types = FALSE
    ) %>%
    mutate(sample = sample_name)

  return(df)
}

# Read Flye assembly info output files and concatenate in one dataframe
assembly_stats <- do.call(
  rbind,
  lapply(X = assembly_stats_files, FUN = read_stats, name_position = 2)) %>%
  select(sample, '#seq_name', length, 'cov.', 'circ.')
colnames(assembly_stats) <- c("sample", "contig", "contig_length", "contig_depth", "circular")
write_delim(
  x = assembly_stats,
  file = here("data", "processed", "assembly_stats-concatenated.tsv"),
  delim = "\t"
)

assembly_stats_summary <- assembly_stats %>%
  group_by(sample) %>%
  summarise(mean_contig_depth = mean(contig_depth),
            total_assembly_length = sum(contig_length))


arg_and_assembly <- left_join(
  x = deduplicated_args,
  y = assembly_stats_summary,
  by = "sample"
) %>%
  left_join(
    x = .,
    y = assembly_stats,
    by = c("sample", "contig")
  )

# Remove dataframes that are no longer necessary to save memory
rm(arg_stats, deduplicated_args, assembly_stats, assembly_stats_summary)

## 3. Taxonomic classifications

read_classification_files <- function(filename) {
  # Cut the sample name from the file path
  sample_name <- str_split_1(string = filename, pattern = "/") %>% tail(2) %>% head(1)

  df <- read_delim(
    file = filename,
    show_col_types = FALSE
  ) %>%
    rename("taxonomy" = `...9`) %>%
    separate_wider_delim(
      cols = taxonomy,
      delim = "\t",
      names = c("Species", "Taxonomy")) %>%
    mutate(sample = sample_name)
  
  return(df)
}

classifications <- do.call(
  rbind,
  lapply(X = classification_files, FUN = read_classification_files)
)
write_delim(
  x = classifications,
  file = here("data", "processed", "classifications-concatenated.tsv"),
  delim = "\t"
)

db_with_classifications <- left_join(
  x = arg_and_assembly,
  y = classifications %>%
    select(sample, readID, Species, Taxonomy),
  by = c("sample", "contig" = "readID")
)

rm(classifications, arg_and_assembly)

## 4. Chromosome/virus/plasmid classifications

genomad_scores <- do.call(
  rbind,
  lapply(X = genomad_scores_files, FUN = read_stats, name_position = 3)
)

plasmid_classifications <- do.call(
  rbind,
  lapply(X = genomad_plasmid_files, FUN = read_stats, name_position = 3)
) %>%
  rename("contig" = "seq_name",
         "plasmid_topology" = "topology",
         "plasmid_genes" = "n_genes")
write_delim(
  x = plasmid_classifications,
  file = here("data", "processed", "plasmid_predictions-concatenated.tsv"),
  delim = "\t"
)

virus_classifications <- do.call(
  rbind,
  lapply(X = genomad_virus_files, FUN = read_stats, name_position = 3)
) %>%
  rename("contig" = "seq_name",
         "virus_topology" = "topology",
         "virus_genes" = "n_genes",
         "virus_taxonomy" = "taxonomy")
write_delim(
  x = virus_classifications,
  file = here("data", "processed", "virus_predictions.tsv"),
  delim = "\t"
)

db_with_mge <- left_join(
  x = db_with_classifications,
  y = genomad_scores,
  by = c("sample", "contig" = "seq_name")
) %>%
  left_join(
    x = .,
    y = plasmid_classifications %>%
      select(sample, contig, plasmid_topology,
             plasmid_genes, conjugation_genes,
             amr_genes),
    by = c("sample", "contig")
  ) %>%
  left_join(
    x = .,
    y = virus_classifications %>%
      select(sample, contig, virus_topology,
             virus_genes, virus_taxonomy)
  ) %>%
  mutate(genomad_prediction = case_when(
    ! is.na(plasmid_topology) ~ "plasmid",
    ! is.na(virus_topology) ~ "virus",
    TRUE ~ "chromosome"
  ))

write_csv(
  x = db_with_mge,
  file = here("data", "processed", "assembly_database.csv")
)
