#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# -----------------------------
# CLI Options
# -----------------------------
option_list <- list(
  make_option(c("-m", "--min_support"), type = "double", default = 0,
              help = "Minimum reads support to report SV [default: %default]", metavar = "NUMBER")
)

# ---- Parse options ----
opt <- parse_args(OptionParser(option_list = option_list))

min_support <- opt$min_support

# Get list of all stellerator files
input_files <- list.files(pattern = "stellerator.*\\.tsv$", full.names = TRUE)

# Extract sample ID from first file (all files are from same sample)
sample_id <- unique(sub("-stellerator.*", "", basename(input_files[1])))

clamp0 <- function(x) {
  pmax(as.numeric(x), 0)
}

process_table <- function(input) {
  table <- read_tsv(input) %>%
    mutate(FUSION = paste(query_gene, matched_partner_gene, sep = "-"),
           BREAKPOINT1 = str_split_i(breakpoint_estimate, "/", 1),
           BREAKPOINT2 = str_split_i(breakpoint_estimate, "/", 2)) %>%
    rename(CHR1 = reference_name, CHR2 = inferred_partner_reference) %>%
    group_by(FUSION, CHR1, BREAKPOINT1, CHR2, BREAKPOINT2) %>%
    summarise(SUPPORT = n(), .groups = "drop") %>%
    filter(SUPPORT >= min_support) %>%
    mutate(TYPE = "SUPPLEMENTARY READS",
           DIRECTION = NA) %>%
    select(FUSION, CHR1, BREAKPOINT1, CHR2, BREAKPOINT2, TYPE, DIRECTION, SUPPORT)

  return(table)
}

process_figeno <- function(input,table) {
  figeno <- read_tsv(input) %>%
    group_by(query_gene) %>%
    summarise(BREAKPOINT1 = min(alignment_start),
              BREAKPOINT2 = max(inferred_partner_start),
              query_gene = unique(query_gene),
              matched_partner_gene = unique(matched_partner_gene),
              reference_name = unique(reference_name),
              inferred_partner_reference = unique(inferred_partner_reference),
              .groups = "drop") %>%
    mutate(alteration = paste(query_gene, matched_partner_gene, sep = "-")) %>%
    filter(alteration %in% table$FUSION) %>%
    mutate(
      pos  = paste0(reference_name, ":", clamp0(BREAKPOINT1 - 30000), "-", clamp0(BREAKPOINT1 + 30000)),
      pos2 = paste0(inferred_partner_reference, ":", clamp0(BREAKPOINT2 - 30000), "-", clamp0(BREAKPOINT2 + 30000)),
      FUSION = paste(sample_id, alteration, sep = "_")
    ) %>%
    select(FUSION, pos, pos2)

    return(figeno)
}

tables <- lapply(input_files,process_table)
all_tables <- bind_rows(tables)
figenos <- lapply(input_files,process_figeno, all_tables)
all_figeno <- bind_rows(figenos)

# Write combined tables with sample ID in filename
write_tsv(all_tables, paste(sample_id, "stellerator_table_fusions.tsv", sep = "_"), col_names = FALSE, quote = "none")
write_tsv(all_figeno, paste(sample_id, "stellerator_figeno_fusions.txt", sep = "_"), col_names = FALSE, quote = "none")
