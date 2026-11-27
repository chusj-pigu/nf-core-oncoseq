#!/usr/bin/env Rscript

# Copyright (c) the nf-core community under an open-source MIT license.

suppressPackageStartupMessages({
    library(optparse)
    library(tidyverse)
})

option_list <- list(
    make_option(c("-n", "--nofilt"), type = "character", default = NULL,
                help = "Path to the bed file with coverage calculated with no filters [default: %default]", metavar = "FILE"),
    make_option(c("-p", "--primary"), type = "character", default = NULL,
                help = "Path to the bed file with coverage calculated with filter for primary alignements [default: %default]", metavar = "FILE"),
    make_option(c("-u", "--unique"), type = "character", default = NULL,
                help = "Path to the bed file with coverage calculated with filter for unique alignments [default: %default]", metavar = "FILE"),
    make_option(c("-b", "--bcov"), type = "double", default = NULL,
                help = "Background coverage [default: %default]", metavar = "NUMBER"),
    make_option(c("-l", "--lowgenes"), type = "character", default = NULL,
                help = "Path to text file containing list of low fidelity genes separated by a newline [default: %default]", metavar = "FILE"),
    make_option(c("-c", "--chrmap"), type = "character", default = NULL,
                help = "Path to TSV mapping BAM chromosome names to UCSC names (from chrom_name_map.sh) [default: %default]", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = "output.pdf",
                help = "Output pdf name [default: %default]", metavar = "FILE")
)

# ---- Parse options ----
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Inputs ----
full_bed_file   <- opt$nofilt
prim_bed_file   <- opt$primary
mapq60_bed_file <- opt$unique
bg_cov <- as.numeric(opt$bcov)
genes_low_fidelity <- readLines(opt$lowgenes)

# ---- Functions ----

# "HG003_panel_nofilter.regions.bed" -> "nofilter"
bed_label <- function(path) {
  sub("\\..*$", "", sub(".*_", "", basename(path)))
}

# Load the chromosome name mapping produced from the BAM header (chrom_name_map.sh):
# tab-separated original_name -> UCSC name, with a "# convention: X" comment line.
# Returns a named character vector (names = original, values = UCSC), or NULL.
load_chrom_map <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  lines <- readLines(path)
  convention <- sub("^# convention: ", "", grep("^# convention:", lines, value = TRUE)[1])
  message("Chromosome naming convention detected from BAM header: ", convention)
  if (convention %in% c("MIXED", "UNKNOWN")) {
    warning("Chromosome naming convention is ", convention, "; check the reference/BAM header.")
  }
  fields <- strsplit(grep("^#", lines, value = TRUE, invert = TRUE), "\t")
  fields <- Filter(function(x) length(x) == 2, fields)
  if (length(fields) == 0) {
    return(NULL)
  }
  setNames(vapply(fields, `[`, "", 2), vapply(fields, `[`, "", 1))
}

# Rename chromosomes using the mapping; names absent from the mapping are kept as-is
normalize_chr <- function(chr) {
  chr <- as.character(chr)
  if (is.null(chrom_map)) {
    return(chr)
  }
  unname(ifelse(chr %in% names(chrom_map), chrom_map[chr], chr))
}

chrom_map <- load_chrom_map(opt$chrmap)

# "0001_chr1_GNB1" / "0001_1_GNB1" / "0001_NC_060925.1_GNB1" -> "GNB1"
# Strips "<index>_<chrom>_". NCBI accessions (NC_xxxxxx.x) contain an underscore,
# so they are matched explicitly; other chromosome names must not contain one.
# The chromosome in the name is not necessarily the row's chromosome (e.g. PAR
# genes on chrY labelled with chrX), so it is matched by pattern, not by value.
strip_gene_prefix <- function(gene) {
  sub("^\\d+_(NC_\\d+(\\.\\d+)?|[^_]+)_", "", as.character(gene))
}

# Read, rename, and de-duplicate gene names in a bed file
process_bed <- function(input_bed) {
  bed <- read.delim(input_bed, header = FALSE) %>%
    dplyr::rename(chr = V1, start = V2, end = V3, gene = V4, coverage = V5) %>%
    mutate(chr = normalize_chr(chr)) %>%
    mutate(gene = strip_gene_prefix(gene)) %>%
    arrange(chr, start, end) %>%
    rename_with(~ bed_label(input_bed), coverage) %>%
    mutate(gene = ifelse(duplicated(gene), paste(gene, chr, sep = "_"), gene))
  return(bed)
}

# Single source of truth for the y-axis / capping limit
compute_limit <- function(maximum) {
  if (!is.finite(maximum) || maximum <= 0) {
    return(0.05)
  }

  limit <- if (maximum >= 10) {
    ceiling(maximum / 10) * 10
  } else if (maximum >= 1) {
    ceiling(maximum)
  } else {
    ceiling(maximum * 10) / 10
  }

  if (limit == 0) {
    limit <- 0.05
  }
  return(limit)
}

detect_outliers_high <- function(bed) {
  bed %>%
    filter(!gene %in% genes_low_fidelity) %>%
    filter(!chr == "chrX" & !chr == "chrY") %>%
    mutate(zscore = (mapq60 - mean(mapq60)) / sd(mapq60)) %>%
    filter(zscore > 2.75) %>%
    pull(gene)
}

detect_outliers_low <- function(bed) {
  bed %>%
    filter(!gene %in% genes_low_fidelity) %>%
    filter(!chr == "chrX" & !chr == "chrY") %>%
    mutate(zscore = (mapq60 - mean(mapq60)) / sd(mapq60)) %>%
    filter(zscore < -2.75) %>%
    pull(gene)
}

normalize_bed <- function(bed, high_coverage_limit) {
  bed <- bed %>%
    mutate(
      nofilter = case_when(
        nofilter > high_coverage_limit ~ high_coverage_limit,
        TRUE ~ nofilter
      ),
      primary = case_when(
        primary > high_coverage_limit ~ high_coverage_limit,
        TRUE ~ primary
      ),
      mapq60 = case_when(
        mapq60 > high_coverage_limit ~ high_coverage_limit,
        TRUE ~ mapq60
      )
    )

  bed <- bed %>%
    mutate(fidelity = ifelse(gene %in% genes_low_fidelity, "Low fidelity",
                       ifelse(gene %in% outliers_high, "Possible increased copy number",
                       ifelse(gene %in% outliers_low, "Possible decreased copy number",
                              "Normal (-2.75 < zscore < 2.75)"))))
  return(bed)
}

# Stacked segments, guarded against negative heights
df_long <- function(bed) {
  bed <- bed %>%
    mutate(primary = pmax(primary - mapq60, 0)) %>%
    mutate(nofilter = pmax(nofilter - mapq60 - primary, 0)) %>%
    pivot_longer(c(nofilter:mapq60), names_to = "set", values_to = "coverage")

  bed$set <- gsub("nofilter", "no filter", bed$set)
  bed$set <- gsub("primary", "primary only", bed$set)
  bed$set <- gsub("mapq60", "unique", bed$set)

  return(bed)
}

# Label a gene only if capping actually clipped one of its values
generate_ann_out <- function(bed_long, bed_norm, bed_all) {
  genes_capped <- bed_all %>%
    inner_join(
      bed_norm %>% select(gene, nofilter_capped = nofilter, primary_capped = primary, mapq60_capped = mapq60),
      by = "gene"
    ) %>%
    filter(nofilter > nofilter_capped | primary > primary_capped | mapq60 > mapq60_capped) %>%
    pull(gene)

  ann <- bed_all %>%
    filter(gene %in% genes_capped) %>%
    pivot_longer(c(nofilter:mapq60), names_to = "set", values_to = "coverage") %>%
    filter(set == "nofilter") %>%
    mutate(ann = paste0(as.character(round(coverage)), "X")) %>%
    select(-coverage) %>%
    left_join(select(bed_long, gene, coverage), by = "gene") %>%
    group_by_at(vars(chr:ann)) %>%
    summarise(coverage = max(coverage), .groups = "drop")

  ann$set <- gsub("nofilter", "no filter", ann$set)
  return(ann)
}

# Annotations labelling median and background coverage
general_ann <- function(bed) {
  chrs <- c("chr5", "chr10", "chr15", "chrY")

  bed1 <- bed %>%
    filter(chr %in% chrs) %>%
    group_by(chr) %>%
    slice(which.max(end)) %>%
    mutate(coverage = round(median), ann = paste0("Median (", round(median), "X)"))
  bed2 <- bed %>%
    filter(chr %in% chrs) %>%
    group_by(chr) %>%
    slice(which.max(end)) %>%
    mutate(coverage = round(bg_cov), ann = paste0("Background (", round(bg_cov), "X)"))

  bed <- rbind(bed1, bed2) %>%
    mutate(set = factor(set, levels = c("unique", "primary only", "no filter"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  return(bed)
}

# `ymax` is the shared limit used to cap the data
generate_plot <- function(bed, ymax, ann_out, ann_facet, output_pdf) {

  axis_ticks <- seq(0, ymax, length.out = 5)
  label_offset <- 0.15 * ymax

  bed <- bed %>%
    mutate(set = factor(set, levels = c("no filter", "primary only", "unique"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  ann_out <- ann_out %>%
    mutate(set = factor(set, levels = c("no filter", "primary only", "unique"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  p <- ggplot() +
    geom_bar(data = bed, aes(x = gene, y = coverage, fill = fidelity, alpha = set), stat = "identity") +
    geom_text(data = ann_out, aes(x = gene, y = coverage - label_offset, label = ann), size = 3, hjust = "inward") +
    geom_text(data = ann_facet, aes(x = gene, y = coverage, label = ann), size = 4, vjust = 0.5, hjust = "outward", nudge_x = 0.5) +
    geom_hline(yintercept = ceiling(median), linewidth = 1, linetype = "dashed") +
    geom_hline(yintercept = ceiling(bg_cov), linewidth = 1, linetype = "dashed") +
    facet_wrap(~ chr, nrow = 5, scales = "free_x") +
    theme(
      plot.margin = unit(c(0.5, 4, 0, 0), "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18),
      strip.text.x = element_text(size = 18),
      strip.background = element_rect(fill = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      panel.background = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 18),
      panel.grid.major = element_line(colour = "grey")
    ) +
    scale_y_continuous(breaks = axis_ticks) +
    scale_fill_manual(values = c("lavenderblush4", "orangered3", "turquoise4", "#DDAA33FF"),
                      breaks = c("Normal (-2.75 < zscore < 2.75)", "Possible increased copy number", "Possible decreased copy number", "Low fidelity")) +
    scale_alpha_manual(values = c(0.33, 0.66, 1), breaks = c("no filter", "primary only", "unique")) +
    coord_cartesian(ylim = c(0, max(axis_ticks)), expand = FALSE, clip = "off") +
    labs(y = "Mean coverage", alpha = "Alignement type filter", fill = "",
         title = paste(sub("^(.*)_.*$", "\\1", full_bed_file), " (mean coverage:", round(mean), "X)"))

  pdf(output_pdf, width = 22, height = 14)
  print(p)
  dev.off()

  svg(gsub("\\.pdf$", ".svg", output_pdf), width = 22, height = 14)
  print(p)
  dev.off()
}

# ---- Usage ----

tryCatch({
  input <- c(full_bed_file, prim_bed_file, mapq60_bed_file)
  bed_list <- lapply(input, process_bed)
  names(bed_list) <- bed_label(input)
  stopifnot(setequal(names(bed_list), c("nofilter", "primary", "mapq60")))

  # Safety check: downstream logic expects UCSC names (chrX/chrY exclusion, facet annotations)
  bed_chrs <- unique(bed_list[["nofilter"]]$chr)
  if (!any(grepl("^chr", bed_chrs))) {
    warning("No UCSC-style chromosome names after mapping (found: ",
            paste(head(bed_chrs, 5), collapse = ", "),
            "). Sex-chromosome exclusion and facet annotations will not work.")
  }

  # Median / mean for primary alignments only
  median <- median(bed_list[["primary"]]$primary)
  mean <- bed_list[["primary"]] %>%
    filter(!gene %in% genes_low_fidelity) %>%
    filter(!chr %in% c("chrY", "chrX")) %>%
    pull(primary) %>%
    mean()

  # Join the three bed files into one data frame
  bed_all <- bed_list[["nofilter"]] %>%
    left_join(select(bed_list[["primary"]], gene, primary), by = "gene") %>%
    left_join(select(bed_list[["mapq60"]], gene, mapq60), by = "gene")

  outliers_high <- detect_outliers_high(bed_all)
  outliers_low  <- detect_outliers_low(bed_all)

  max_normal_coverage <- suppressWarnings(
    max(bed_all$nofilter[bed_all$nofilter < 1.5 * median])
  )
  if (!is.finite(max_normal_coverage)) {
    max_normal_coverage <- 0
  }

  plot_limit <- compute_limit(max_normal_coverage)

  bed_norm  <- normalize_bed(bed_all, plot_limit)
  bed_long  <- df_long(bed_norm)
  ann_df    <- generate_ann_out(bed_long, bed_norm, bed_all)
  ann_facet <- general_ann(bed_long)

  generate_plot(bed_long, plot_limit, ann_df, ann_facet, opt$output)

  means_df <- data.frame(sample = sub("_mapq.pdf", "", opt$output),
                         panel_cov = mean,
                         background_cov = bg_cov)

  write.csv(means_df, paste0(means_df$sample, "_mean_coverage.csv"), quote = FALSE, row.names = FALSE)

}, error = function(e) {
  warning("[WARNING] No plot generated: ", e$message)
})
