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
              help = "Path to target csv file", metavar = "file"),
  make_option(c("-e", "--exclude"), type = "character",
              help = "Path blacklist of artefact", metavar = "file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Please provide an input TSV file with --input", call. = FALSE)
}

input <- opt$input
target_list <- opt$target
sample_id <- sub("_filt\\.tsv$", "", basename(input))
blacklist <- opt$exclude

# -----------------------------
# Helper Functions
# -----------------------------

clean_genes <- function(x) {
  # Remove leading/trailing or double separators
  x <- str_replace_all(x, "^&|&$", "")
  x <- str_replace_all(x, "&{2,}", "&")

  # Split by '&' and keep first and last gene
  parts <- str_split(x, "&", simplify = FALSE)

  sapply(parts, function(p) {
    p <- p[p != ""]  # remove empty
    if (length(p) == 0) return(NA_character_)
    if (length(p) == 1) return(p)
    paste0(p[1], "-", p[length(p)])
  })
}

clamp0 <- function(x) {
  pmax(as.numeric(x), 0)
}

extract_gene <- function(x) {
  str_extract(x, "(?<=\\|(HIGH|MODERATE)\\|)[^|]+")
}

# Detect direction based on ALT
get_fusion_direction <- function(alt) {
  dplyr::case_when(
    str_detect(alt, "^[ATCGN]\\[") ~ "->->",    # like A[CHR2:123[
    str_detect(alt, "^[ATCGN]\\]") ~ "-><-",    # like A]CHR2:123]
    str_detect(alt, "^\\]") ~ "<-<-",           # like ]CHR2:123]A
    str_detect(alt, "^\\[") ~ "<-->",           # like [CHR2:123[A
    TRUE ~ NA_character_
  )
}

# Extract the chromosome from ALT (e.g. ]CHR7:152401117]T → chr7)
get_chr_from_alt <- function(alt) {
  chr <- str_extract(alt, "CHR\\d+|CHRX|CHRY")

  case_when(
    chr == "CHRX" ~ "chrX",
    chr == "CHRY" ~ "chrY",
    TRUE ~ str_to_lower(chr)   # CHR1 -> chr1, CHR12 -> chr12, etc.
  )
}

process_variant <- function(df, type = c("BND", "DEL_INS"), window_bnd = 30000, window_delins = 20000) {
  type <- match.arg(type)

  df <- df %>%
    mutate(
      GENE = extract_gene(X8),
      START = as.numeric(X2),
      SUPPORT = as.numeric(str_extract(X8, "(?<=SUPPORT=)\\d+")),
      ALT = X5
    )

  if (type == "BND") {
    df <- df %>%
      filter(str_detect(X3, "BND")) %>%
      mutate(
        direction = get_fusion_direction(ALT),
        chr_alt = get_chr_from_alt(ALT),
        END = as.numeric(str_extract(ALT, "(?<=:)\\d+(?=[]\\[])")),
        CHR1 = X1,
        CHR2 = chr_alt,
        BREAKPOINT1 = START,
        BREAKPOINT2 = END
      ) %>%
      mutate(
        pos = paste0(
            CHR1, ":",
            clamp0(BREAKPOINT1 - window_bnd), "-",
            clamp0(BREAKPOINT1 + window_bnd)
        ),
        pos2 = paste0(
            CHR2, ":",
            clamp0(BREAKPOINT2 - window_bnd), "-",
            clamp0(BREAKPOINT2 + window_bnd)
        ),
        TYPE = str_extract(X8, "(?<=\\|)[^|]+(?=\\|(HIGH|MODERATE)\\|)"),
        GENE = clean_genes(GENE)
      )
  } else {
    df <- df %>%
      filter(!str_detect(X3, "BND")) %>%
      mutate(
        END = as.numeric(str_extract(X8, "(?<=END=)\\d+")),
        LEN = abs(END-START),
        TYPE = str_extract(X8, "(?<=SVTYPE=)[^;]+"),
        GENE = clean_genes(GENE),
        GENE = ifelse(GENE == "", X1, GENE),
        pos = paste0(
            X1, ":",
            clamp0(START - window_bnd), "-",
            clamp0(END + window_bnd)
        ),
      ) %>%
      filter(LEN <= 1000000)
  }

  return(df)
}


# -----------------------------
# Read Input
# -----------------------------
vcf <- read_tsv(input, col_names = FALSE)
targets_df <- read.csv(target_list)

# -----------------------------
# Process Variants
# -----------------------------
if (nrow(vcf) > 0) {
  vcf_bnd <- process_variant(vcf, type = "BND") %>%
    # Collapse duplicate fusions: keep lowest BREAKPOINT1, highest BREAKPOINT2
    group_by(GENE, CHR1, CHR2) %>%
    summarise(
      BREAKPOINT1 = min(BREAKPOINT1),
      BREAKPOINT2 = max(BREAKPOINT2),
      SUPPORT     = max(SUPPORT),
      direction   = first(direction),
      TYPE        = first(TYPE),
      .groups     = "drop"
    ) %>%
    mutate(
      pos  = paste0(CHR1, ":", clamp0(BREAKPOINT1 - 30000), "-", clamp0(BREAKPOINT1 + 30000)),
      pos2 = paste0(CHR2, ":", clamp0(BREAKPOINT2 - 30000), "-", clamp0(BREAKPOINT2 + 30000)),
      GENE = paste(sample_id, GENE, sep = "_")
    ) %>%
    select(GENE, pos, pos2)

  vcf_del_ins <- process_variant(vcf, type = "DEL_INS") %>%
    # Collapse duplicate genes: keep lowest START, highest END
    group_by(GENE, X1, TYPE) %>%
    summarise(
      START   = min(START),
      END     = max(END),
      SUPPORT = max(SUPPORT),
      .groups = "drop"
    ) %>%
    mutate(
      LEN  = abs(END - START),
      pos  = paste0(X1, ":", clamp0(START - 20000), "-", clamp0(END + 20000)),
      GENE = paste(sample_id, GENE, sep = "_")
    ) %>%
    filter(LEN < 1000000) %>%
    select(GENE, pos)

  delin_table <- process_variant(vcf, type = "DEL_INS") %>%
    group_by(GENE, X1, TYPE) %>%
    summarise(
      START   = min(START),
      END     = max(END),
      SUPPORT = max(SUPPORT),
      .groups = "drop"
    ) %>%
    mutate(LEN = abs(END - START)) %>%
    filter(LEN < 1000000) %>%
    rename(CHR = X1) %>%
    select(CHR, GENE, TYPE, SUPPORT, START, END, LEN)

  bnd_table <- process_variant(vcf, type = "BND") %>%
    group_by(GENE, CHR1, CHR2) %>%
    summarise(
      BREAKPOINT1 = min(BREAKPOINT1),
      BREAKPOINT2 = max(BREAKPOINT2),
      SUPPORT     = max(SUPPORT),
      direction   = first(direction),
      TYPE        = first(TYPE),
      .groups     = "drop"
    ) %>%
    rename(FUSION = GENE) %>%
    select(FUSION, CHR1, BREAKPOINT1, CHR2, BREAKPOINT2, TYPE, direction, SUPPORT) %>%
    arrange(FUSION, CHR1)

} else {
  vcf_bnd     <- vcf
  vcf_del_ins <- vcf
  delin_table <- vcf
  bnd_table   <- vcf
}

# -----------------------------
# Add gene to sample_id column
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
  missing_rows <- targets_df %>%
    rowwise() %>%
    filter(!any(str_detect(detected_genes, paste0("\\b", GENE, "\\b")))) %>%
    ungroup()
} else {
  missing_rows <- data.frame()
}

# Append missing rows to vcf_del_ins if any and match GENE to vcf_del_ins
if (nrow(missing_rows) > 0) {
  target_final <- missing_rows %>%
    mutate(GENE=paste(sample_id, GENE, sep = "_"))
  message("Added ", nrow(missing_rows), " missing target genes")
} else {
    target_final <- data.frame()
}

# --------------------------------------------
# Exclude blacklisted genes that are artefacts
# ---------------------------------------------

# Define SV blacklist suffixes

unwanted_calls <- readLines(blacklist)
unwanted_pattern <- "(^|_|-)LOC|(^|_|-)LINC"

# -----------------------------
# Save to output
# -----------------------------
safe_write_figeno <- function(df, file) {
  if (nrow(df) > 0) {
    
    if ("GENE" %in% colnames(df)) {
      df <- df %>%
        filter(!str_detect(GENE, unwanted_pattern)) %>%
        filter(!GENE %in% paste(sample_id, unwanted_calls, sep = "_"))
    } else if ("FUSION" %in% colnames(df)) {
      df <- df %>%
        filter(!str_detect(FUSION, unwanted_pattern)) %>%
        filter(!FUSION %in% paste(sample_id, unwanted_calls, sep = "_"))
    }
    if (nrow(df > 0)) {
      write_tsv(df, file, col_names = FALSE, quote = "none")
    } else {
      message(paste("Skipping", file, "No SV remaining after filtering out blacklist"))
    }
  } else {
    message(paste("Skipping", file, "- dataframe is empty"))
  }
}

safe_write_tables <- function(df, file) {
  if (nrow(df) > 0) {
    
    loc_fusion_pattern <- paste0("^", sample_id, "_LOC|-LOC|_LINC|-LINC")
    
    if ("GENE" %in% colnames(df)) {
      df <- df %>%
        filter(!str_detect(GENE, unwanted_pattern)) %>%
        filter(!GENE %in% unwanted_calls)
    } else if ("FUSION" %in% colnames(df)) {
      df <- df %>%
        filter(!str_detect(FUSION, unwanted_pattern)) %>%
        filter(!FUSION %in% unwanted_calls)
    }

    write_tsv(df, file, col_names = FALSE, quote = "none")
  } else {
    write_tsv(df, file, col_names = FALSE, quote = "none")
  }
}


safe_write_figeno(vcf_bnd, paste(sample_id, "region_fusions.txt", sep = "_"))
safe_write_figeno(vcf_del_ins, paste(sample_id, "region_indel.txt", sep = "_"))
safe_write_figeno(target_final, paste(sample_id, "targets_nohit.txt", sep = "_"))
safe_write_tables(bnd_table, paste(sample_id, "table_fusions.tsv", sep = "_"))
safe_write_tables(delin_table, paste(sample_id, "table_indel.tsv", sep = "_"))
