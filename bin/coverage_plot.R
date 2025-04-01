#!/usr/bin/env Rscript

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
    make_option(c("-b", "--bcov"), type = "integer", default = NULL,
                            help = "Background coverage [default: %default]", metavar = "NUMBER"),
    make_option(c("-l", "--lowgenes"), type = "character", default = NULL,
                            help = "Path to text file containing list of low fidelity genes separated by a newline [default: %default]", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = "output.pdf",
                            help = "Output pdf name [default: %default]", metavar = "FILE")
)

# ---- Parse options ----
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Read the bed files ----

full_bed_file <- opt$nofilt
prim_bed_file <- opt$primary
mapq60_bed_file <- opt$unique
bg_cov <- as.numeric(opt$bcov)
genes_low_fidelity <- readLines(opt$lowgenes)

# Functions ####

# Helper function to read, process and de-duplicate gene names in bed files

process_bed <- function(input_bed) {
  bed <- read.delim(input_bed, header = FALSE) %>%
    dplyr::rename(chr = V1, start = V2, end = V3, gene = V4, coverage = V5) %>%
    arrange(chr,start,end) %>%
    rename_with(~ gsub(".*_(.*)\\.regions\\.bed$", "\\1", input_bed), coverage) %>%
    mutate(gene = gsub("^\\d{3}_.+?_", "", gene)) %>%
    mutate(gene = ifelse(duplicated(gene), paste(gene, chr, sep = "_"), gene))

  return(bed)
}

detect_outliers_high <- function(bed) {

  outliers <- bed_all %>%
    filter(!gene %in% genes_low_fidelity) %>%
    filter(!chr == "chrX" & !chr == "chrY") %>%
    mutate(zscore = (mapq60-mean(mapq60))/sd(mapq60)) %>%
    filter(zscore > 2.75) %>%
    pull(gene)

  return(outliers)

}

detect_outliers_low <- function(bed) {

  outliers <- bed_all %>%
    filter(!gene %in% genes_low_fidelity) %>%
    filter(!chr == "chrX" & !chr == "chrY") %>%
    mutate(zscore = (mapq60-mean(mapq60))/sd(mapq60)) %>%
    filter(zscore < -2.75) %>%
    pull(gene)

  return(outliers)

}

normalize_bed <- function(bed, maximum) {

  # Normalize outliers in coverage
  high_coverage_limit <- (ceiling(maximum / 10) * 10)

  bed <- bed %>%
    mutate(nofilter = case_when(
      nofilter > high_coverage_limit & primary > high_coverage_limit ~ (ceiling(maximum / 10) * 10),
      nofilter > high_coverage_limit ~ (ceiling(maximum / 10) * 10),
      TRUE ~ nofilter
    ),
    primary = case_when(
      primary > high_coverage_limit & mapq60 > high_coverage_limit ~ (ceiling(maximum / 10) * 10),
      primary > high_coverage_limit ~ (ceiling(maximum / 10) * 10),
      TRUE ~ primary
    ),
    mapq60 = case_when(
      mapq60 > high_coverage_limit ~ (ceiling(maximum / 10) * 10),
      TRUE ~ mapq60
    ))

  # Add fidelity variable for coloring
  bed <- bed %>%
    mutate(fidelity = ifelse(gene %in% genes_low_fidelity, "Low fidelity", ifelse(gene %in% outliers_high, "Possible increased copy number", ifelse(gene %in% outliers_low, "Possible decreased copy number", "Normal (-2.75 < zscore < 2.75)"))))
  return(bed)
}

df_long <- function(bed) {
  bed <- bed %>%
    mutate(primary = primary - mapq60) %>%
    mutate(nofilter = nofilter - mapq60 - primary) %>%
    pivot_longer(c(nofilter:mapq60), names_to = "set", values_to = "coverage")

  #Rename the sets with more informative names
  bed$set <- gsub("nofilter", "no filter", bed$set)
  bed$set <- gsub("primary", "primary only", bed$set)

  return(bed)
}

# Function to dynamically identify and annotate genes with coverage > 2× median
generate_ann_out <- function(bed_long, bed_all) {

  ann <- bed_all %>%
    filter(gene %in% genes_high) %>%
    pivot_longer(c(nofilter:mapq60), names_to = "set", values_to = "coverage") %>%
    filter(set == "nofilter") %>%
    mutate(ann = paste0(as.character(round(coverage)), "X")) %>%
    select(-coverage) %>%
    left_join(select(filter(bed_long),gene,coverage)) %>%
    group_by_at(vars(chr:ann)) %>%
    summarise(coverage = max(coverage)) %>%
    ungroup()

  ann$set <- gsub("nofilter", "no filter", ann$set)

  return(ann)

}

# Make annotation to label median and background coverage
general_ann <- function(bed) {
  bed1 <- bed %>%
    filter(chr == "chr8" | chr == "chr16" | chr == "chrY") %>%
    group_by(chr) %>%
    slice(which.max(end)) %>%
    mutate(coverage = round(median), ann = paste0("Median (", round(median), "X)"))
  bed2 <- bed %>%
    filter(chr == "chr8" | chr == "chr16" | chr == "chrY") %>%
    group_by(chr) %>%
    slice(which.max(end)) %>%
    mutate(coverage = round(bg_cov), ann = paste0("Background (", round(bg_cov), "X)"))

  bed <- rbind(bed1,bed2) %>%
    mutate(set = factor(set, levels = c("mapq60", "primary only", "no filter"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  return(bed)

}

# Function to generate a coverage plot
generate_plot <- function(bed, maximum, ann_out, ann_facet, output_pdf) {
  # Calculate axis parameters
  axis_ticks <- seq(0, (ceiling(maximum / 10) * 10), length.out = 5)

  # Reorder chromosomes for plotting

  bed <- bed %>%
    mutate(set = factor(set, levels = c("no filter", "primary only", "mapq60"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  ann_out <- ann_out %>%
    mutate(set = factor(set, levels = c("no filter", "primary only", "mapq60"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  # Plot and save as PDF
  pdf(output_pdf, width = 22, height = 14)
  print( ggplot() +
           geom_bar(data = bed, aes(x = gene, y = coverage, fill = fidelity, alpha = set), stat = "identity") +
           geom_text(data=ann_out, aes(x = gene, y = coverage-50, label = ann), size =4, hjust = "inward") +
           geom_text(data=ann_facet, aes(x = gene, y = coverage, label = ann), size =4, vjust = 0.5, hjust = "outward", nudge_x = 0.5) +
           geom_hline(yintercept = ceiling(median), linewidth = 1, linetype = 'dashed') +
           geom_hline(yintercept = ceiling(bg_cov), linewidth = 1, linetype = 'dashed') +
           facet_wrap(~ chr, nrow = 3, scales = "free_x") +
           theme(
             plot.margin = unit(c(0.5,4,0,0), "cm"),
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
           scale_y_continuous(limits = c(0, max(axis_ticks)), breaks = axis_ticks) +
           scale_fill_manual(values = c("lavenderblush4","orangered3" ,"turquoise4", "#DDAA33FF"), breaks = c("Normal (-2.75 < zscore < 2.75)", "Possible increased copy number", "Possible decreased copy number", "Low fidelity")) +
           scale_alpha_manual(values = c(0.33,0.66,1), breaks = c("no filter", "primary only", "mapq60")) +
           coord_cartesian(expand = FALSE, clip = "off") +
           labs(y = "Mean coverage", alpha = "Alignement type filter", fill = "", title = sub("^(.*)_.*$", "\\1", full_bed_file)) )
  dev.off()
}
# Usage ####

## Join input bed files together in a list
input <- c(full_bed_file, prim_bed_file, mapq60_bed_file)
bed_list <- lapply(input, process_bed)
names(bed_list) <- gsub(".*_(.*)\\.regions\\.bed$", "\\1", input)

# Store median for pirmary alignment only as a variable for future usage
median <- median(bed_list[["primary"]]$primary)

# Join bed_files into one dataframe:
bed_all <- bed_list[[1]] %>%
  left_join(select(bed_list[[2]], c(4,5)), by = "gene") %>%
  left_join(select(bed_list[[3]], c(4,5)), by = "gene")

outliers_high <- detect_outliers_high(bed_all)
outliers_low <- detect_outliers_low(bed_all)

genes_high <- bed_all %>% filter(mapq60 > 1.5*median | primary > 1.5*median | nofilter > 1.5*median) %>% pull(gene)

# Execute functions to make data ready for plotting:
max_normal_coverage <- max(bed_all$nofilter[bed_all$nofilter < 1.5 * median])
bed_norm <- normalize_bed(bed_all,max_normal_coverage)
bed_long <- df_long(bed_norm)
ann_df <- generate_ann_out(bed_long, bed_all)
ann_facet <- general_ann(bed_long)

# Generate the plot
generate_plot(bed_long, max_normal_coverage, ann_df, ann_facet, opt$output)
