#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# -----------------------------
# CLI Options
# -----------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Path to input TSV file (required)", metavar = "file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Please provide an input TSV file with --input", call. = FALSE)
}

input <- opt$input
sample_id <- sub("_filt\\.tsv$", "", basename(input))

# -----------------------------
# Helper Functions
# -----------------------------
make_range <- function(x, window = 30000) {
  parsed <- str_match(x, "CHR(\\w+):(\\d+)")
  if (is.na(parsed[1])) return(NA_character_)
  paste0("chr", parsed[,2], ":", as.numeric(parsed[,3]) - window, "-", as.numeric(parsed[,3]) + window)
}

clean_genes <- function(x) {
  x <- str_remove_all(x, "(&?LOC\\d+&?)")
  x <- str_replace_all(x, "^&|&$", "")
  str_replace_all(x, "&{2,}", "&")
}

extract_gene <- function(x) {
  str_extract(x, "(?<=\\|(HIGH|MODERATE)\\|)[^|]+")
}

process_variant <- function(df, type = c("BND", "DEL_INS"), window_bnd = 30000, window_delins = 20000) {
  type <- match.arg(type)

  df <- df %>%
    mutate(
      GENE = extract_gene(X8),
      START = as.numeric(X2)
    )

  if (type == "BND") {
    df <- df %>%
      filter(str_detect(X3, "BND")) %>%
      mutate(
        pos = paste0(X1, ":", START - window_bnd, "-", START + window_bnd),
        pos2 = make_range(X5),
        GENE = paste(sample_id, str_replace_all(GENE, "&", "-"), sep = "-")
      ) %>%
      select(GENE, pos, pos2)

  } else {
    df <- df %>%
      filter(!str_detect(X3, "BND")) %>%
      mutate(
        END = as.numeric(str_extract(X8, "(?<=END=)\\d+")),
        GENE = clean_genes(GENE),
        GENE = ifelse(GENE == "", X1, GENE),
        pos = paste0(X1, ":", START - window_delins, "-", END + window_delins),
        GENE = paste(sample_id, GENE, sep = "-")
      ) %>%
      select(GENE, pos)
  }

  return(df)
}

# -----------------------------
# Read Input
# -----------------------------
vcf <- read_tsv(input, col_names = FALSE)

# -----------------------------
# Process Variants
# -----------------------------
vcf_bnd     <- process_variant(vcf, type = "BND")
vcf_del_ins <- process_variant(vcf, type = "DEL_INS")

# -----------------------------
# Save to output
# -----------------------------
safe_write <- function(df, file) {
  if (nrow(df) > 0) {
    write_tsv(df, file, col_names = FALSE, quote = "none")
  } else {
    message(paste("Skipping", file, "- dataframe is empty"))
  }
}

safe_write(vcf_bnd, "region_fusions.txt")
safe_write(vcf_del_ins, "region_indel.txt")
