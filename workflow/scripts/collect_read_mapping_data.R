#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

coverage_files <- snakemake@input[["cov_files"]]

read_coverage_file <- function(coverage_file) {
  # Strip sample name from file name
  sample <- basename(coverage_file) %>% gsub(pattern = "-coverage.tsv", replacement = "", x = .)

  # Read the dataframe from tab-separated file
  coverage_df <- read_delim(file = coverage_file, show_col_types = FALSE)

  # Insert sample names
  if (nrow(coverage_df) == 0) {
    return(data.frame())
  } else {
    coverage_df <- cbind(Sample = sample, coverage_df %>%
      rename("Contig" = "#rname"))
  }

  return(coverage_df)
}

complete_coverage_df <- do.call(
  rbind,
  lapply(X = coverage_files, FUN = read_coverage_file)
)

write_csv(x = complete_coverage_df, file = snakemake@output[[1]])
