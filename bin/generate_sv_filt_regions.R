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
              help = "Path to input TSV file (required)", metavar = "file"),
  make_option(c("-t", "--target"), type = "character",
              help = "Path to target csv file", metavar = "file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Please provide an input TSV file with --input", call. = FALSE)
}

input <- opt$input
target_list <- opt$target
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
  # remove LOC followed by digits
  x <- str_remove_all(x, "LOC\\d+")
  # remove any double/trailing separators created by removal
  x <- str_replace_all(x, "^&|&$", "")
  x <- str_replace_all(x, "&{2,}", "&")
  # keep only the first gene if there are multiple
  x <- str_replace(x, "&.*$", "")
  x
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
        GENE = paste(sample_id, str_replace_all(GENE, "&", "-"), sep = "_")
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
        GENE = paste(sample_id, GENE, sep = "_")
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
if (nrow(vcf) > 0) {
    vcf_bnd     <- process_variant(vcf, type = "BND")
    vcf_del_ins <- process_variant(vcf, type = "DEL_INS")
} else {
    vcf_bnd <- vcf
    vcf_del_ins <- vcf
}

# -----------------------------
# Extract genes from sample_id column
# -----------------------------
extract_genes_bnd <- function(x) {
  # Remove sample_id prefix
  genes_part <- sub(".*_", "", x)
  # Split on "-" if there are multiple genes
  str_split(genes_part, "-", simplify = FALSE) %>% unlist()
}

extract_genes_delins <- function(x) {
  # For DEL/INS, only one gene per row
  sub(".*_", "", x)
}

# Initialize detected_genes vector
detected_genes <- character(0)

# Only extract genes if dataframes are not empty
if (nrow(vcf_bnd) > 0) {
  detected_genes <- c(detected_genes, extract_genes_bnd(vcf_bnd[[1]]))
}

if (nrow(vcf_del_ins) > 0) {
  detected_genes <- c(detected_genes, extract_genes_delins(vcf_del_ins[[1]]))
}

detected_genes <- unique(detected_genes)

# Identify missing genes only if target list not empty:

if (nrow(targets_df > 0)) {
  missing_genes <- setdiff(targets_df$GENE, detected_genes)


  missing_rows <- targets_df %>%
    filter(GENE %in% missing_genes)
} else {
  missing_rows <- data.frame()
}

# Append missing rows to vcf_del_ins if any and match GENE to vcf_del_ins
if (nrow(missing_rows) > 0) {
  missing_rows <- missing_rows %>%
    mutate(GENE=paste(sample_id, GENE, sep = "_"))
  vcf_del_ins <- bind_rows(vcf_del_ins, missing_rows)
  message("Appended ", nrow(missing_rows), " missing target genes to vcf_del_ins")
}


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

safe_write(vcf_bnd, paste(sample_id, "region_fusions.txt", sep = "_"))
safe_write(vcf_del_ins, paste(sample_id, "region_indel.txt", sep = "_"))
