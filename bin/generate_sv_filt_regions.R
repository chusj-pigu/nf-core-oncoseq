#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# -----------------------------
# CLI Options
# -----------------------------
option_list <- list(
  make_option(c("-t", "--target"), type = "character",
              help = "Path to target csv file", metavar = "file"),
  make_option(c("-e", "--exclude"), type = "character",
              help = "Path blacklist of artefact", metavar = "file"),
  make_option(c("-p", "--panel"), type = "character",
              help = "Path adaptive panel", metavar = "file")
)

opt <- parse_args(OptionParser(option_list = option_list))

input <- list.files(path = ".", pattern = "\\.tsv$")
target_list <- opt$target
sample_id <- unique(sub("_(severus|sniffles)_filt\\.tsv$","",basename(input)))
blacklist <- opt$exclude
panel <- opt$panel

# -----------------------------
# Helper Functions
# -----------------------------

process_bed <- function(input_bed) {
  bed <- read.delim(input_bed, header = FALSE) %>%
    dplyr::rename(chr = V1, start = V2, end = V3, gene = V4) %>%
    mutate(gene = gsub("^\\d{4}_.+?_", "", gene)) %>%
    mutate(gene = ifelse(duplicated(gene), paste(gene, chr, sep = "_"), gene)) %>%
    pull(gene)

  return(bed)
}

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

extract_support <- function(x,y) {
  index_dv = str_which(str_split(x, pattern = ":")[[1]], pattern = "DV")
  str_split_i(y, pattern = ":", i = index_dv)
}

extract_partner <- function(x) {
  gsub("MATE_ID=", "", str_split_i(x, pattern = ";", i = 4))
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
  chr <- str_extract(str_to_lower(alt), "chr\\d+|chrX|chrY")
}

# extract sv type according to program used
get_type_from_alt <- function(id,alt) {
  case_when(str_detect(id, "Sniffles") ~ str_extract(alt, "(?<=\\|)[^|]+(?=\\|(HIGH|MODERATE)\\|)"),
            str_detect(id, "severus") ~ gsub(".*BND_TYPE=([^;]+);.*", "\\1", alt))
}

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

process_vcf <- function(vcf) {
    df <- vcf %>%
        mutate(
        GENE = extract_gene(X8),
        ID = gsub("(_BND\\d+)_\\d+$", "\\1", X3),
        START = as.numeric(X2),
        SUPPORT = as.numeric(extract_support(X9,X10)),
        ALT = X5,
        INFO2 = ifelse(str_detect(X3, "severus"),str_split_i(X8, pattern = ",", i = 2), GENE),
        GENE2 = ifelse(str_detect(X3, "Sniffles"), str_split_i(GENE, pattern = "&", i = 2),
                        ifelse(is.na(INFO2), NA, extract_gene(INFO2)))
        ) %>%
        mutate(
        GENE    = ifelse(str_detect(GENE, "&"), str_split_i(GENE, pattern = "&", i = 1), GENE),
        strand1 = str_split_i(gsub(".*STRANDS?=([^;]+);.*", "\\1", X8), pattern = "", i=1),
        strand2 = str_split_i(gsub(".*STRANDS?=([^;]+);.*", "\\1", X8), pattern = "", i=2)) %>%
        filter(str_detect(X8, "HIGH|MODERATE")) %>%
        filter(str_detect(GENE, pattern_start) | str_detect(GENE, pattern_end) |
                str_detect(GENE2, pattern_start) | str_detect(GENE2, pattern_end))
    return(df)
}

parse_bnd <- function(df) {
    bnd <- df %>%
        filter(str_detect(ID, "BND")) %>%
        mutate(
            direction = get_fusion_direction(ALT),
            chr_alt = get_chr_from_alt(ALT),
            END = as.numeric(str_extract(ALT, "(?<=:)\\d+(?=[]\\[])")),
            CHR1 = X1,
            CHR2 = chr_alt,
            BREAKPOINT1 = START,
            BREAKPOINT2 = END,
            TYPE = get_type_from_alt(ID,X8),
            FUSION = ifelse(is.na(GENE2), GENE, paste(GENE,GENE2, sep = "-"))) %>%
        filter(!FUSION %in% unwanted_calls) %>%
        ## Collapse rows that represent same fusion
        rowwise() %>%
        mutate(
            chr_min  = min(CHR1, CHR2),
            chr_max  = max(CHR1, CHR2),
            bp_min   = min(BREAKPOINT1, BREAKPOINT2),
            bp_max   = max(BREAKPOINT1, BREAKPOINT2)
        ) %>%
        ungroup() %>%
        distinct(FUSION, chr_min, chr_max, bp_min, bp_max, .keep_all = TRUE) %>%
        select(-chr_min, -chr_max, -bp_min, -bp_max)
    return(bnd)
}

parse_other_sv <- function(df) {
    other <- df %>%
        filter(!str_detect(ID, "BND")) %>%
        mutate(
        END = as.numeric(str_extract(X8, "(?<=END=)\\d+")),
        LEN = abs(END-START),
        TYPE = str_extract(X8, "(?<=SVTYPE=)[^;]+"),
        GENE = clean_genes(GENE),
        GENE = ifelse(GENE == "", X1, GENE)) %>%
        filter(!GENE %in% unwanted_calls)
    return(other)
}

summary_table_bnd <- function(df) {
    bnd <- df %>%
        select(FUSION,CHR1,BREAKPOINT1,CHR2,BREAKPOINT2,TYPE,direction,SUPPORT)
    return(bnd)
}

summary_table_other <- function(df) {
    other <- df %>%
        select(X1,GENE,TYPE,START,END,LEN,SUPPORT)
    return(other)
}

region_figeno_bnd <- function(bnd) {
    figeno_bnd <- bnd %>%
        group_by(FUSION, CHR1, CHR2) %>%
        summarise(
            BREAKPOINT1 = min(BREAKPOINT1),
            BREAKPOINT2 = max(BREAKPOINT2),
            .groups     = "drop"
        ) %>%
        mutate(
            pos  = paste0(CHR1, ":", clamp0(BREAKPOINT1 - 30000), "-", clamp0(BREAKPOINT1 + 30000)),
            pos2 = paste0(CHR2, ":", clamp0(BREAKPOINT2 - 30000), "-", clamp0(BREAKPOINT2 + 30000)),
            FUSION = paste(sample_id, FUSION, sep = "_")
        ) %>%
        select(FUSION, pos, pos2)
    return(figeno_bnd)
}

region_figeno_other <- function(other) {
    figeno_other <- other %>%
        group_by(X1,GENE) %>%
        summarise(
            START       = min(START),
            END         = max(END),
            .groups     = "drop"
        ) %>%
        mutate(
            pos = paste0(
            X1, ":",
            clamp0(START - 20000), "-",
            clamp0(END + 20000)
            ),
            GENE = paste(sample_id, GENE, sep = "_"),
            LEN  = END - START
        ) %>%
        filter(LEN <= 1000000) %>%
        select(GENE,pos)
    return(figeno_other)
}

input_sv_figeno <- function(df) {
    figeno_table <- df %>%
        mutate(
            chr2 = case_when(is.na(get_chr_from_alt(ALT)) ~ X1,
                                TRUE ~ get_chr_from_alt(ALT)),
            pos2 = case_when(is.na(str_extract(ALT, "(?<=:)\\d+(?=[]\\[])")) ~ as.numeric(str_extract(X8, "(?<=END=)\\d+")),
                            TRUE ~ as.numeric(str_extract(ALT, "(?<=:)\\d+(?=[]\\[])"))),
            svtype = str_extract(X8, "(?<=SVTYPE=)[^;]+")
        ) %>%
        mutate(
            color = case_when(svtype == "BND" ~ "#27ae60",
                                svtype == "DEL" ~ "#4a69bd",
                                svtype == "INV" ~ "#8e44ad",
                                svtype == "DUP" ~ "#e55039",
                                TRUE ~ "#4a69bd")
        ) %>%
        select(X1, START, chr2, pos2, strand1, strand2, color, svtype) %>%
        rename(chr1 = X1, pos1 = START)
    return(figeno_table)
}

# -----------------------------
# Process Variants
# -----------------------------

vcf <- lapply(input, function(f) {
    read_tsv(f, col_names = FALSE, skip = 1, show_col_types = FALSE)
})

names(vcf) <- sapply(input, function(f) {
    if (grepl("severus", basename(f))) {
        "severus"
    } else if (grepl("sniffles", basename(f))) {
        "sniffles"
    }
})

targets_df <- read.csv(target_list)
unwanted_calls <- readLines(blacklist)
genes_panel <- process_bed(panel)

# Pattern to select genes in panel
pattern_end <- paste0(
  "-(", paste(genes_panel, collapse = "|"), ")$"  # gene at the end: x-GENE
)

pattern_start <- paste0(
  "^(",
  paste(genes_panel, collapse = "|"),
  ")(-|$|_)"
)

vcf <- vcf[vapply(vcf, nrow, integer(1)) > 0]

if (length(vcf) > 0) {

    df <- lapply(vcf, process_vcf)
    df <- df[vapply(df, nrow, integer(1)) > 0]

    bnd <- lapply(df, parse_bnd)
    other <- lapply(df, parse_other_sv)

    bnd_merged   <- bind_rows(bnd)
    other_merged <- bind_rows(other)

    if (nrow(bnd_merged) > 0 & nrow(other_merged) > 0) {
        figeno_bnd <- region_figeno_bnd(bnd_merged)
        figeno_other <- region_figeno_other(other_merged)

        figeno <- bnd_merged %>%
            select(ALT,X1,X8,START,strand1,strand2) %>%
            rbind(select(other_merged,ALT,X1,X8,START,strand1,strand2))
        figeno_table <- input_sv_figeno(figeno)
    } else if (nrow(bnd_merged) > 0) {
        figeno_bnd <- region_figeno_bnd(bnd_merged)
        figeno_other <- other_merged

        figeno <- bnd_merged %>%
            select(ALT,X1,X8,START,strand1,strand2)
        figeno_table <- input_sv_figeno(figeno)
    } else if (nrow(other_merged) > 0) {
        figeno_other <- region_figeno_other(other_merged)
        figeno_bnd <- bnd_merged

        figeno <- other_merged %>%
            select(ALT,X1,X8,START,strand1,strand2)
        figeno_table <- input_sv_figeno(figeno)
    } else {
        figeno_bnd <- bnd_merged
        figeno_other <- other_merged

        figeno_table <- bnd_merged
    }

    table_bnd <- lapply(bnd, summary_table_bnd)
    table_other <- lapply(other, summary_table_other)

} else {
    figeno_bnd <- vcf
    figeno_other <- vcf
    table_bnd <- vcf
    table_other <- vcf
    figeno_table <- vcf
}

# -----------------------------
# Add data frame for gene of interest not present in vcf
# -----------------------------

# Initialize detected_genes vector
detected_genes <- character(0)

# Only extract genes if dataframes are not empty
if (nrow(figeno_bnd) > 0) {
  detected_genes <- c(detected_genes, extract_genes_bnd(figeno_bnd[[1]]))
}

if (nrow(figeno_other) > 0) {
  detected_genes <- c(detected_genes, extract_genes_delins(figeno_other[[1]]))
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

# -----------------------------
# Save to output
# -----------------------------
safe_write_figeno <- function(df, file) {
  if (nrow(df) > 0) {
      write_tsv(df, file, col_names = FALSE, quote = "none")
  } else {
    message(paste("Skipping", file, "- dataframe is empty"))
  }
}

safe_write_table <- function(lst, suffix) {

    expected <- c("severus", "sniffles")

    empty_df <- data.frame()

    for (nm in expected) {

        outfile <- paste(sample_id, nm, suffix, sep = "_")

        df <- lst[[nm]]

        # replace NULL or missing with empty df
        if (is.null(df) || !is.data.frame(df)) {
            df <- empty_df
        }

        # write (empty or not)
        write_tsv(
            df,
            outfile,
            col_names = FALSE,
            quote = "none"
        )
    }
}

safe_write_figeno_table <- function(df, file) {
  if (nrow(df) > 0) {
      write_tsv(df, file, col_names = FALSE, quote = "none")
  } else {
      df <- data.frame(chr1 = "", pos1 = "", chr2 = "",
        pos2 = "", strand1 = "", strand2 = "", color = "", svtype = "")
    write_tsv(df, file, col_names = TRUE, quote = "none")
  }
}

safe_write_figeno(figeno_bnd, paste(sample_id, "region_fusions.txt", sep = "_"))
safe_write_figeno(figeno_other, paste(sample_id, "region_other.txt", sep = "_"))
safe_write_figeno(target_final, paste(sample_id, "targets_nohit.txt", sep = "_"))
safe_write_table(table_bnd, "table_fusions.tsv")
safe_write_table(table_other, "table_other.tsv")
safe_write_figeno_table(figeno_table, paste(sample_id, "table_figeno.tsv", sep = "_"))

