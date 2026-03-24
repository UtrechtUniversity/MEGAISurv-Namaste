#!/usr/bin/env Rscript

suppressPackageStartupMessages(
  library(tidyverse)
)

sink(
  file = file(snakemake@log[[1]], open = "wt"),
  type = "message"
)

kma_reads_files <- snakemake@input[["results"]]

read_kma_output <- function(kma_file) {
  # Extract the sample name from the file name
  sample <- basename(kma_file) %>%
    gsub(
      pattern = ".hmm.res", # excluding the (double) suffix
      replacement = "",
      x = .
    )

  df <- read_delim(
    kma_file,
    show_col_types = FALSE,
    delim = "\t",
    trim_ws = TRUE
  ) %>%
    rename("Gene" = "#Template") %>%
    mutate(Sample = sample) %>% # Fill in the sample
    select(Sample, everything()) # And make sample the first column

  return(df)
}

kma_reads_df <- do.call(
  rbind,
  lapply(X = kma_reads_files, FUN = read_kma_output)
)

# Apply the same selection as for the assembly-based database
filtered_args <- kma_reads_df %>%
  filter(Template_Coverage >= 60)


write_delim(
  x = filtered_args,
  file = snakemake@output[[1]],
  delim = "\t"
)
