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
    make_option(c("-b", "--bcov"), type = "double", default = NULL,
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
    rename_with(~ gsub(".*_(.*)\\.bed", "\\1", input_bed), coverage) %>%
    mutate(gene = gsub("^\\d{4}_.+?_", "", gene)) %>%
    mutate(gene = ifelse(duplicated(gene), paste(gene, chr, sep = "_"), gene))

  return(bed)
}

# FIX 1: single source of truth for the y-axis / capping limit.
# This used to be computed twice, with different rounding rules, which meant
# genes capped at one value could fall outside the axis range and disappear.
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

# Takes the shared limit directly instead of recomputing its own.
normalize_bed <- function(bed, high_coverage_limit) {

  bed <- bed %>%
    mutate(nofilter = case_when(
      nofilter > high_coverage_limit & primary > high_coverage_limit ~ high_coverage_limit,
      nofilter > high_coverage_limit ~ high_coverage_limit,
      TRUE ~ nofilter
    ),
    primary = case_when(
      primary > high_coverage_limit & mapq60 > high_coverage_limit ~ high_coverage_limit,
      primary > high_coverage_limit ~ high_coverage_limit,
      TRUE ~ primary
    ),
    mapq60 = case_when(
      mapq60 > high_coverage_limit ~ high_coverage_limit,
      TRUE ~ mapq60
    ))

  bed <- bed %>%
    mutate(fidelity = ifelse(gene %in% genes_low_fidelity, "Low fidelity",
                       ifelse(gene %in% outliers_high, "Possible increased copy number",
                       ifelse(gene %in% outliers_low, "Possible decreased copy number",
                              "Normal (-2.75 < zscore < 2.75)"))))
  return(bed)
}

# FIX 4: guard the stacked segments against negative heights. After capping (or
# with non-monotonic input beds where nofilter < primary) these differences can
# go below zero, which stacks downward and looks like a missing or inverted bar.
df_long <- function(bed) {
  bed <- bed %>%
    mutate(primary = pmax(primary - mapq60, 0)) %>%
    mutate(nofilter = pmax(nofilter - mapq60 - primary, 0)) %>%
    pivot_longer(c(nofilter:mapq60), names_to = "set", values_to = "coverage")

  #Rename the sets with more informative names
  bed$set <- gsub("nofilter", "no filter", bed$set)
  bed$set <- gsub("primary", "primary only", bed$set)
  bed$set <- gsub("mapq60", "unique", bed$set)

  return(bed)
}

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
    make_option(c("-b", "--bcov"), type = "double", default = NULL,
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
    rename_with(~ gsub(".*_(.*)\\.bed", "\\1", input_bed), coverage) %>%
    mutate(gene = gsub("^\\d{4}_.+?_", "", gene)) %>%
    mutate(gene = ifelse(duplicated(gene), paste(gene, chr, sep = "_"), gene))

  return(bed)
}

compute_limit <- function(maximum) {

  if (!is.finite(maximum) || maximum <= 0) {
    return(0.05)
  }

  limit <- if (maximum >= 10) {
    ceiling(maximum / 10) * 10
  } else if (maximum >= 1) {
    ceiling(maximum / 5) * 5
  } else {
    ceiling(maximum * 10) / 10
  }

  if (limit == 0) {
    limit <- 0.05
  }

  return(limit)
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

# Takes the shared limit directly instead of recomputing its own.
normalize_bed <- function(bed, high_coverage_limit) {

  bed <- bed %>%
    mutate(nofilter = case_when(
      nofilter > high_coverage_limit & primary > high_coverage_limit ~ high_coverage_limit,
      nofilter > high_coverage_limit ~ high_coverage_limit,
      TRUE ~ nofilter
    ),
    primary = case_when(
      primary > high_coverage_limit & mapq60 > high_coverage_limit ~ high_coverage_limit,
      primary > high_coverage_limit ~ high_coverage_limit,
      TRUE ~ primary
    ),
    mapq60 = case_when(
      mapq60 > high_coverage_limit ~ high_coverage_limit,
      TRUE ~ mapq60
    ))

  bed <- bed %>%
    mutate(fidelity = ifelse(gene %in% genes_low_fidelity, "Low fidelity",
                       ifelse(gene %in% outliers_high, "Possible increased copy number",
                       ifelse(gene %in% outliers_low, "Possible decreased copy number",
                              "Normal (-2.75 < zscore < 2.75)"))))
  return(bed)
}

# FIX 4: guard the stacked segments against negative heights. After capping (or
# with non-monotonic input beds where nofilter < primary) these differences can
# go below zero, which stacks downward and looks like a missing or inverted bar.
df_long <- function(bed) {
  bed <- bed %>%
    mutate(primary = pmax(primary - mapq60, 0)) %>%
    mutate(nofilter = pmax(nofilter - mapq60 - primary, 0)) %>%
    pivot_longer(c(nofilter:mapq60), names_to = "set", values_to = "coverage")

  #Rename the sets with more informative names
  bed$set <- gsub("nofilter", "no filter", bed$set)
  bed$set <- gsub("primary", "primary only", bed$set)
  bed$set <- gsub("mapq60", "unique", bed$set)

  return(bed)
}

# Function to dynamically identify and annotate genes with coverage > 2× median
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
    filter(chr == "chr5" | chr == "chr10" | chr == "chr15" | chr == "chrY") %>%
    group_by(chr) %>%
    slice(which.max(end)) %>%
    mutate(coverage = round(median), ann = paste0("Median (", round(median), "X)"))
  bed2 <- bed %>%
    filter(chr == "chr5" | chr == "chr10" | chr == "chr15" | chr == "chrY") %>%
    group_by(chr) %>%
    slice(which.max(end)) %>%
    mutate(coverage = round(bg_cov), ann = paste0("Background (", round(bg_cov), "X)"))

  bed <- rbind(bed1,bed2) %>%
    mutate(set = factor(set, levels = c("unique", "primary only", "no filter"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  return(bed)

}

# Function to generate a coverage plot
# `ymax` is the shared limit used to cap the data, so the axis can never cut off
# a bar that was deliberately capped.
generate_plot <- function(bed, ymax, ann_out, ann_facet, output_pdf) {

  axis_ticks <- seq(0, ymax, length.out = 5)

  # FIX 3: place the out-of-range labels relative to the axis rather than at a
  # hardcoded -50, so they stay on canvas at 2X as well as at 500X.
  label_offset <- 0.15 * ymax

  # Reorder chromosomes for plotting

  bed <- bed %>%
    mutate(set = factor(set, levels = c("no filter", "primary only", "unique"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  ann_out <- ann_out %>%
    mutate(set = factor(set, levels = c("no filter", "primary only", "unique"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  # Build the plot once, then render it to both devices.
  p <- ggplot() +
    geom_bar(data = bed, aes(x = gene, y = coverage, fill = fidelity, alpha = set), stat = "identity") +
    geom_text(data = ann_out, aes(x = gene, y = coverage - label_offset, label = ann), size = 3, hjust = "inward") +
    geom_text(data = ann_facet, aes(x = gene, y = coverage, label = ann), size = 4, vjust = 0.5, hjust = "outward", nudge_x = 0.5) +
    geom_hline(yintercept = ceiling(median), linewidth = 1, linetype = 'dashed') +
    geom_hline(yintercept = ceiling(bg_cov), linewidth = 1, linetype = 'dashed') +
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
    # FIX 2: breaks on the scale, range on the coord. Setting `limits` on the
    # scale silently turns out-of-range values into NA before the bars are
    # stacked, which is what made the capped gene vanish entirely.
    scale_y_continuous(breaks = axis_ticks) +
    scale_fill_manual(values = c("lavenderblush4", "orangered3", "turquoise4", "#DDAA33FF"),
                      breaks = c("Normal (-2.75 < zscore < 2.75)", "Possible increased copy number", "Possible decreased copy number", "Low fidelity")) +
    scale_alpha_manual(values = c(0.33, 0.66, 1), breaks = c("no filter", "primary only", "unique")) +
    coord_cartesian(ylim = c(0, max(axis_ticks)), expand = FALSE, clip = "off") +
    labs(y = "Mean coverage", alpha = "Alignement type filter", fill = "",
         title = paste(sub("^(.*)_.*$", "\\1", full_bed_file), " (mean coverage:", round(mean), "X)"))

  # Plot and save as PDF
  pdf(output_pdf, width = 22, height = 14)
  print(p)
  dev.off()

  # Generate a svg plot alongside the pdf plot
  svg(gsub("\\.pdf$", ".svg", output_pdf), width = 22, height = 14)
  print(p)
  dev.off()

}

# Usage ####

tryCatch({
  ## Join input bed files together in a list
  input <- c(full_bed_file, prim_bed_file, mapq60_bed_file)
  bed_list <- lapply(input, process_bed)
  names(bed_list) <- gsub(".*_(.*)\\.bed", "\\1", input)

  # Store median for pirmary alignment only as a variable for future usage
  median <- median(bed_list[["primary"]]$primary)
  mean <- bed_list[["primary"]] %>%
      filter(!gene %in% genes_low_fidelity) %>%
      filter(!chr %in% c("chrY", "chrX")) %>%
      pull(primary) %>%
      mean()

  # Join bed_files into one dataframe:
  bed_all <- bed_list[[1]] %>%
    left_join(select(bed_list[[2]], c(4,5)), by = "gene") %>%
    left_join(select(bed_list[[3]], c(4,5)), by = "gene")

  outliers_high <- detect_outliers_high(bed_all)
  outliers_low <- detect_outliers_low(bed_all)

  # Execute functions to make data ready for plotting:
  max_normal_coverage <- suppressWarnings(
    max(bed_all$nofilter[bed_all$nofilter < 1.5 * median])
  )
  if (!is.finite(max_normal_coverage)) {
    max_normal_coverage <- 0
  }

  plot_limit <- compute_limit(max_normal_coverage)

  bed_norm <- normalize_bed(bed_all, plot_limit)
  bed_long <- df_long(bed_norm)
  ann_df <- generate_ann_out(bed_long, bed_norm, bed_all)
  ann_facet <- general_ann(bed_long)

  # Generate the plot
  generate_plot(bed_long, plot_limit, ann_df, ann_facet, opt$output)

  means_df <- data.frame(sample = sub("_mapq.pdf", "", opt$output),
                        panel_cov = mean,
                        background_cov = bg_cov)

  write.csv(means_df, paste0(means_df$sample,"_mean_coverage.csv"), quote = FALSE, row.names = FALSE)

}, error = function(e) {
  warning("[WARNING] No plot generated: ", e$message)
}
)

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
    mutate(set = factor(set, levels = c("unique", "primary only", "no filter"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  return(bed)

}

# Function to generate a coverage plot
# `ymax` is the shared limit used to cap the data, so the axis can never cut off
# a bar that was deliberately capped.
generate_plot <- function(bed, ymax, ann_out, ann_facet, output_pdf) {

  axis_ticks <- seq(0, ymax, length.out = 5)

  # FIX 3: place the out-of-range labels relative to the axis rather than at a
  # hardcoded -50, so they stay on canvas at 2X as well as at 500X.
  label_offset <- 0.15 * ymax

  # Reorder chromosomes for plotting

  bed <- bed %>%
    mutate(set = factor(set, levels = c("no filter", "primary only", "unique"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  ann_out <- ann_out %>%
    mutate(set = factor(set, levels = c("no filter", "primary only", "unique"))) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE))) %>%
    mutate(gene = factor(gene, levels = unique(gene)))

  # Build the plot once, then render it to both devices.
  p <- ggplot() +
    geom_bar(data = bed, aes(x = gene, y = coverage, fill = fidelity, alpha = set), stat = "identity") +
    geom_text(data = ann_out, aes(x = gene, y = coverage - label_offset, label = ann), size = 4, hjust = "inward") +
    geom_text(data = ann_facet, aes(x = gene, y = coverage, label = ann), size = 4, vjust = 0.5, hjust = "outward", nudge_x = 0.5) +
    geom_hline(yintercept = ceiling(median), linewidth = 1, linetype = 'dashed') +
    geom_hline(yintercept = ceiling(bg_cov), linewidth = 1, linetype = 'dashed') +
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
    # FIX 2: breaks on the scale, range on the coord. Setting `limits` on the
    # scale silently turns out-of-range values into NA before the bars are
    # stacked, which is what made the capped gene vanish entirely.
    scale_y_continuous(breaks = axis_ticks) +
    scale_fill_manual(values = c("lavenderblush4", "orangered3", "turquoise4", "#DDAA33FF"),
                      breaks = c("Normal (-2.75 < zscore < 2.75)", "Possible increased copy number", "Possible decreased copy number", "Low fidelity")) +
    scale_alpha_manual(values = c(0.33, 0.66, 1), breaks = c("no filter", "primary only", "unique")) +
    coord_cartesian(ylim = c(0, max(axis_ticks)), expand = FALSE, clip = "off") +
    labs(y = "Mean coverage", alpha = "Alignement type filter", fill = "",
         title = paste(sub("^(.*)_.*$", "\\1", full_bed_file), " (mean coverage:", round(mean), "X)"))

  # Plot and save as PDF
  pdf(output_pdf, width = 22, height = 14)
  print(p)
  dev.off()

  # Generate a svg plot alongside the pdf plot
  svg(gsub("\\.pdf$", ".svg", output_pdf), width = 22, height = 14)
  print(p)
  dev.off()

}

# Usage ####

tryCatch({
  ## Join input bed files together in a list
  input <- c(full_bed_file, prim_bed_file, mapq60_bed_file)
  bed_list <- lapply(input, process_bed)
  names(bed_list) <- gsub(".*_(.*)\\.bed", "\\1", input)

  # Store median for pirmary alignment only as a variable for future usage
  median <- median(bed_list[["primary"]]$primary)
  mean <- bed_list[["primary"]] %>%
      filter(!gene %in% genes_low_fidelity) %>%
      filter(!chr %in% c("chrY", "chrX")) %>%
      pull(primary) %>%
      mean()

  # Join bed_files into one dataframe:
  bed_all <- bed_list[[1]] %>%
    left_join(select(bed_list[[2]], c(4,5)), by = "gene") %>%
    left_join(select(bed_list[[3]], c(4,5)), by = "gene")

  outliers_high <- detect_outliers_high(bed_all)
  outliers_low <- detect_outliers_low(bed_all)

  genes_high <- bed_all %>% filter(mapq60 > 1.5*median | primary > 1.5*median | nofilter > 1.5*median) %>% pull(gene)

  # Execute functions to make data ready for plotting:
  max_normal_coverage <- suppressWarnings(
    max(bed_all$nofilter[bed_all$nofilter < 1.5 * median])
  )
  if (!is.finite(max_normal_coverage)) {
    max_normal_coverage <- 0
  }

  # Computed once, used both for capping and for the axis.
  plot_limit <- compute_limit(max_normal_coverage)

  bed_norm <- normalize_bed(bed_all, plot_limit)
  bed_long <- df_long(bed_norm)
  ann_df <- generate_ann_out(bed_long, bed_all)
  ann_facet <- general_ann(bed_long)

  # Generate the plot
  generate_plot(bed_long, plot_limit, ann_df, ann_facet, opt$output)

  means_df <- data.frame(sample = sub("_mapq.pdf", "", opt$output),
                        panel_cov = mean,
                        background_cov = bg_cov)

  write.csv(means_df, paste0(means_df$sample,"_mean_coverage.csv"), quote = FALSE, row.names = FALSE)

}, error = function(e) {
  warning("[WARNING] No plot generated: ", e$message)
}
)