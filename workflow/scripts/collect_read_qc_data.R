#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

sink(
  file = file(snakemake@log[[1]], open = "wt"),
  type = "message"
)

read_qc_files <- snakemake@input[["json"]]

json_to_df <- function(json_file) {
  # Get the sample name from the file name (remove '.json' extension)
  sample <- basename(json_file) %>% gsub(pattern = ".json", replacement = "", x = .)

  # Read the JSON,
  json_data <- read_json(path = json_file)
  df <- cbind(
    json_data %>%
      # extract the 'summary' part from before...
      pluck("summary", 2) %>%
      as_tibble() %>%
      rename_with(~ paste0("before_", .x)),
    json_data %>%
      # ... and after filtering.
      pluck("summary", 3) %>%
      as_tibble() %>%
      rename_with(~ paste0("after_", .x))
  )

  # Bind the columns together in one dataframe
  df <- cbind(sample = sample, df)

  return(df)
}

complete_qc_df <- do.call(
  rbind,
  lapply(X = read_qc_files, FUN = json_to_df)
)

write_csv(x = complete_qc_df, file = snakemake@output[[1]])
