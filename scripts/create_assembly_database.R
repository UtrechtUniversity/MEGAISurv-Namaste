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

read_arg_stats <- function(filename) {
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
  mutate(
    hit_length = hit_end - hit_start + 1,
    arg_short = gsub(
      pattern = "_.*",
      replacement = "",
      x = arg
    )
  ) %>%
  select(-c("hit_start", "hit_end"))

# Deduplicate samples by removing runs that derive from the same biosample
deduplicated_args <- arg_stats %>%
  filter(!sample %in% duplicates$run)
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

# Further annotate antibiotics resistances using medically important classes
# assigned by the WHO:
# https://cdn.who.int/media/docs/default-source/gcp/who-mia-list-2024-lv.pdf
# (In particular tables 1-3)
hpcia_antibiotics <- c("cefcapene", "cefdinir", "cefditoren", "cefepime", "cefetamet", "cefixime", "cefmenoxime", "cefodizime", "cefoperazone", "cefoselis", "cefotaxime", "cefovecin", "cefozopran", "cefpiramide", "cefpirome", "cefpodoxime", "cefquinome", "cefsulodin", "ceftazidime", "cefteram pivoxil", "ceftibuten", "ceftiofur", "ceftizoxime", "ceftolozane", "ceftriaxone", "latamoxef", "besifloxacin", "cinoxacin", "ciprofloxacin", "danofloxacin", "delafloxacin", "difloxacin", "enoxacin", "enrofloxacin", "fleroxacin", "flumequine", "garenoxacin", "gatifloxacin", "gemifloxacin", "grepafloxacin", "ibafloxacin", "lascufloxacin", "levofloxacin", "levonadifloxacin", "lomefloxacin", "marbofloxacin", "moxifloxacin", "nadifloxacin", "nalidixic acid", "nemonoxacin", "norfloxacin", "ofloxacin", "orbifloxacin", "oxolinic acid", "ozenoxacin pazufloxacin", "pefloxacin", "pipemidic acid", "piromidic acid", "pradofloxacin", "prulifloxacin", "rosoxacin", "rufloxacin", "sitafloxacin", "sparfloxacin", "temafloxacin", "trovafloxacin", "colistin", "polymyxin B", "fosfomycin")

cia_antibiotics <- c("amikacin", "apramycin", "arbekacin", "astromicin", "bekanamycin", "dibekacin", "dihydrostreptomycin", "framycetin", "gentamicin", "isepamicin", "kanamycin", "micronomicin neomycin", "netilmicin", "paromomycin", "ribostamycin", "sisomicin", "streptoduocin", "streptomycin", "tobramycin", "rifabutin", "rifampicin", "rifamycin", "rifapentine", "rifaximin", "azithromycin", "cethromycin", "clarithromycin", "dirithromycin", "erythromycin", "flurithromycin", "gamithromycin", "josamycin", "kitasamycin", "midecamycin", "miocamycin", "oleandomycin", "rokitamycin", "roxithromycin", "spiramycin", "tildipirosin", "tilmicosin", "troleandomycin", "tulathromycin", "tylosin", "tylvalosin")

hia_antibiotics <- c("chloramphenicol", "florfenicol", "thiamphenicol", "cefacetrile", "cefaclor", "cefadroxil", "cefalexin", "cefalonium", "cefaloridine", "cefalotin", "cefamandole", "cefapirin", "cefatrizine", "cefazedone", "cefazolin", "cefbuperazone", "cefmetazole", "cefminox", "cefonicid", "ceforanide", "cefotetan", "cefotiam", "cefoxitin", "cefprozil", "cefradine", "cefroxadine", "ceftezole", "cefuroxime", "flomoxef", "loracarbef", "clindamycin", "lincomycin", "pirlimycin", "metronidazole", "ornidazole", "ronidazole", "secnidazole", "tinidazole", "mecillinam", "pivmecillinam", "amoxicillin", "ampicillin", "azidocillin", "bacampicillin", "epicillin", "hetacillin", "metampicillin", "pivampicillin", "sultamicillin", "talampicillin", "temocillin", "amoxicillin-clavulanic acid", "ampicillin-sulbactam", "cloxacillin", "flucloxacillin", "meticillin", "methicillin", "nafcillin", "oxacillin", "benethamine-benzylpenicillin", "benzathine-benzylpenicillin", "penicillin G", "clometocillin", "penamecillin", "penethamate", "hydriodide", "pheneticillin", "phenoxymethylpenicillin", "penicillin V", "procaine benzylpenicillin", "propicillin", "pristinamycin", "quinupristin-", "dalfopristin", "virginiamycin", "brodimoprim", "formosulfathiazole", "iclaprim", "ormetoprim", "phthalylsulfathiazole", "pyrimethamine", "sulfachlorpyridazine", "sulfadiazine", "sulfadimethoxine", "sulfadimidine", "sulfafurazole (=", "sulfisoxazole)", "sulfaguanidin", "sulfaisodimidine", "sulfalene", "sulfamazone", "sulfamerazine", "sulfamethazine", "sulfamethizole", "sulfamethoxazole", "sulfamethoxypyridazine", "sulfametomidine", "sulfametoxydiazine", "sulfametrole", "sulfamoxole", "sulfanilamide", "sulfaperin", "sulfaphenazole", "sulfapyridine", "sulfaquinoxaline", "sulfathiazole", "sulfathiourea", "tetroxoprim", "trimethoprim", "fusidic acid", "chlortetracycline", "clomocycline", "demeclocycline", "doxycycline", "lymecycline", "metacycline", "minocycline", "oxytetracycline", "penimepicycline", "rolitetracycline", "sarecycline", "tetracycline")

ia_antibiotics <- c("spectinomycin", "bacitracin", "enramycin", "enduramycin", "methenamine hippurate", "methenamine", "mandelate", "furaltadone", "furazidine", "furazolidone", "nifuroxazide", "nifurtoinol", "nitrofural", "nitrofurantoin", "lefamulin retapamulin", "tiamulin", "valnemulin")

not_important <- c("novobiocin", "nitarsone", "roxarsone", "bicozamycin", "halquinol", "laidlomycin", "lasalocid", "maduramicin", "monensin", "narasin", "salinomycin", "semduramicin", "avilamycin", "bambermycin", "flavomycin", "flavophospholipol", "moenomycin", "carbadox", "olaquindox")

find_any_match <- function(phenotype, class) {
  return(
    any(
      str_split(tolower(phenotype), pattern = ", ") %>%
        unlist() %in% class
    )
  )
}

# Add WHO classifications one by one
deduplicated_args <- deduplicated_args %>%
  mutate(
    hpcia = sapply(
      X = deduplicated_args$Phenotype,
      FUN = find_any_match,
      class = hpcia_antibiotics
    ),
    cia = sapply(
      X = Phenotype,
      FUN = find_any_match,
      class = cia_antibiotics
    ),
    hia = sapply(
      X = Phenotype,
      FUN = find_any_match,
      class = hia_antibiotics
    ),
    ia = sapply(
      X = Phenotype,
      FUN = find_any_match,
      class = ia_antibiotics
    ),
    not_important = sapply(
      X = Phenotype,
      FUN = find_any_match,
      class = not_important
    ),

    # And finally combine them into one column:
    medical_importance = case_when(
      hpcia ~ "HPCIA",
      cia ~ "CIA",
      hia ~ "HIA",
      ia ~ "IA",
      not_important ~ "not_important",
      TRUE ~ "unknown"
    )
  ) %>%
  # And remove intermediate columns
  select(
    -c(hpcia, cia, hia, ia, not_important)
  )


## 2. Assembly stats (length, depth and circularity)

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

# Read Flye assembly info output files and concatenate in one dataframe
assembly_stats <- do.call(
  rbind,
  lapply(X = assembly_stats_files, FUN = read_stats, name_position = 2)
) %>%
  select(sample, "#seq_name", length, "cov.", "circ.")
colnames(assembly_stats) <- c("sample", "contig", "contig_length", "contig_depth", "circular")
write_delim(
  x = assembly_stats,
  file = here("data", "processed", "assembly_stats-concatenated.tsv"),
  delim = "\t"
)

assembly_stats_summary <- assembly_stats %>%
  group_by(sample) %>%
  summarise(
    mean_contig_depth = mean(contig_depth),
    total_assembly_length = sum(contig_length)
  )


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
  sample_name <- str_split_1(string = filename, pattern = "/") %>%
    tail(2) %>%
    head(1)

  df <- read_delim(
    file = filename,
    show_col_types = FALSE
  ) %>%
    rename("taxonomy" = `...9`) %>%
    separate_wider_delim(
      cols = taxonomy,
      delim = "\t",
      names = c("Species", "Taxonomy")
    ) %>%
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
  rename(
    "contig" = "seq_name",
    "plasmid_topology" = "topology",
    "plasmid_genes" = "n_genes"
  )
write_delim(
  x = plasmid_classifications,
  file = here("data", "processed", "plasmid_predictions-concatenated.tsv"),
  delim = "\t"
)

virus_classifications <- do.call(
  rbind,
  lapply(X = genomad_virus_files, FUN = read_stats, name_position = 3)
) %>%
  rename(
    "contig" = "seq_name",
    "virus_topology" = "topology",
    "virus_genes" = "n_genes",
    "virus_taxonomy" = "taxonomy"
  )
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
      select(
        sample, contig, plasmid_topology,
        plasmid_genes, conjugation_genes,
        amr_genes
      ),
    by = c("sample", "contig")
  ) %>%
  left_join(
    x = .,
    y = virus_classifications %>%
      select(
        sample, contig, virus_topology,
        virus_genes, virus_taxonomy
      )
  ) %>%
  mutate(genomad_prediction = case_when(
    !is.na(plasmid_topology) ~ "plasmid",
    !is.na(virus_topology) ~ "virus",
    TRUE ~ "chromosome"
  ))

# Somehow, duplicated rows appear in the dataframe,
# so save the 'distinct' version:
write_csv(
  x = db_with_mge %>% distinct(),
  file = here("data", "processed", "assembly_database.csv")
)
