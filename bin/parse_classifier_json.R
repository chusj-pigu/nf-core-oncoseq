#!/usr/bin/env Rscript
# Copyright (c) the nf-core community under an open-source MIT license.
#
# Usage: Rscript classy_tables.R --sample sample1 --json sample1_combined_classification.json

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(purrr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
})

TOP_N <- 5

# =============================================================
# 1. Inputs
# =============================================================

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-s", "--sample"), type = "character", help = "Sample ID used to name outputs"),
  make_option(c("-j", "--json"),   type = "character", help = "Combined classy JSON", metavar = "file")
)))

sample_id <- opt$sample
parsed    <- fromJSON(opt$json, simplifyVector = FALSE)
scores    <- parsed$combined$classification_scores   # one entry per (classifier, class)

classifier_ids <- scores %>% map_chr("classifier_run_id") %>% unique()


# =============================================================
# 2. Small helpers
# =============================================================

pass_fail <- function(passed) if_else(passed, "PASS", "FAIL")

# Scores for one classifier (a list of entries)
scores_for <- function(id) keep(scores, ~ .x$classifier_run_id == id)

# Hierarchy of one entry as a named character vector,
# e.g. c(Family = "LYMPHOMA", Type = "B-cell lymphoblastic lymphoma")
path_levels <- function(x) {
  kinds <- map_chr(x$native_path, "level_kind")
  names <- map_chr(x$native_path, "name")
  setNames(names, str_to_title(kinds))
}

# Name of one level in an entry's hierarchy, or NA
level_value <- function(x, kind) {
  hit <- keep(x$native_path, ~ .x$level_kind == kind)
  if (length(hit) == 0) NA_character_ else hit[[1]]$name
}

# Turn a list of records into a named numeric lookup
# e.g. lookup(parsed$inference$lineage_probabilities, "label", "probability")
lookup <- function(records, name_field, value_field) {
  setNames(map_dbl(records, ~ as.numeric(.x[[value_field]])),
           map_chr(records, name_field))
}


# =============================================================
# 3. Table builders (one per kind of classifier)
# =============================================================

# ---- Generic: any classifier whose hierarchy is in native_path ----
# Columns: <one per hierarchy level>, Score, Cutoff
make_generic_table <- function(id) {
  scores_for(id) %>%
    map_dfr(function(x) {
      as_tibble(as.list(path_levels(x))) %>%
        mutate(Score  = x$score,
               Cutoff = pass_fail(x$score >= x$cutoff))
    })
}


# ---- MARLIN: adds lineage and group scores ----
# Columns: Lineage, Lineage_score, Group, Group_score, Class, Class_score, Cutoff
make_marlin_table <- function() {
  lineage_scores <- lookup(parsed$inference$lineage_probabilities, "label", "probability")
  group_scores   <- lookup(parsed$inference$group_probabilities,   "label", "probability")

  scores_for("marlin") %>%
    map_dfr(function(x) {
      lineage <- level_value(x, "lineage")
      group   <- level_value(x, "mcf")
      tibble(
        Lineage       = lineage,
        Lineage_score = unname(lineage_scores[lineage]),
        Group         = group,
        Group_score   = unname(group_scores[group]),
        Class         = x$leaf_name,
        Class_score   = x$score,
        Cutoff        = pass_fail(x$score >= x$cutoff)
      )
    }) %>%
    arrange(desc(Class_score)) %>%
    head(TOP_N)
}


# ---- Tucan: family and class scores come from the sampling block ----
# Columns: Family, Family_score, Class, Class_score, Cutoff
make_tucan_table <- function() {
  tucan         <- parsed$tucan
  sampling      <- tucan$samplings[[1]]
  family_scores <- lookup(sampling$family_scores, "family", "score")

  sampling$scores %>%
    map_dfr(~ tibble(
      Family      = .x$family,
      Family_score = unname(family_scores[.x$family]),
      Class       = .x$name,
      Class_score = as.numeric(.x$score)
    )) %>%
    mutate(Cutoff = pass_fail(Class_score >= tucan$cutoff)) %>%
    arrange(desc(Class_score)) %>%
    head(TOP_N)
}


# ---- Optic: PASS/FAIL comes from accepted_call, not score >= cutoff ----
# Optic only stores leaf probabilities, so the score shown at each level is the
# sum of the leaf probabilities underneath it (the "unconditional mass").
# Columns: Level1, Level1_score, Level2, Level2_score, ..., Score, Cutoff
make_optic_table <- function() {
  rows   <- scores_for("optic")
  paths  <- map(rows, ~ map_chr(.x$native_path, "name"))
  leaf_p <- map_dbl(rows, "score")
 
  node_score <- function(prefix) {
    n      <- length(prefix)
    inside <- map_lgl(paths, ~ length(.x) >= n && identical(.x[seq_len(n)], prefix))
    sum(leaf_p[inside])
  }
 
  top   <- order(-leaf_p)[seq_len(min(TOP_N, length(rows)))]
  depth <- max(lengths(paths[top]))
 
  map_dfr(top, function(i) {
    path <- paths[[i]]
    row  <- list()
    for (d in seq_len(depth)) {
      has_level <- d <= length(path)
      row[[paste0("Level", d)]]            <- if (has_level) path[d] else NA_character_
      row[[paste0("Level", d, "_score")]]  <- if (has_level) round(node_score(path[seq_len(d)]), 4) else NA_real_
    }
    row$Cutoff <- pass_fail(isTRUE(rows[[i]]$accepted_call))
    as_tibble(row)
  })
}


# =============================================================
# 4. Build and write everything
# =============================================================

build_table <- function(id) {
  switch(id,
    marlin = make_marlin_table(),
    tucan  = make_tucan_table(),
    optic  = make_optic_table(),
    make_generic_table(id)          # everything else
  )
}

write_table <- function(table, name) {
  file <- paste0(sample_id, "_", str_replace_all(name, ":", "_"), ".tsv")
  write_tsv(table, file)
}

for (id in classifier_ids) {
  message("Processing: ", id)
  write_table(build_table(id), id)
}