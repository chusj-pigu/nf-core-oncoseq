#!/usr/bin/env Rscript

library(readr)
library(jsonlite)
library(optparse)
library(purrr)
library(stringr)
library(dplyr)

# -----------------------------
# CLI Options
# -----------------------------
option_list <- list(
  make_option(c("-s", "--sample"), type = "character",
              help = "Sample ID to name outputs with"),
  make_option(c("-j", "--json"), type = "character",
              help = "Path to json containing all classifier results from classy", metavar = "file")
)

opt <- parse_args(OptionParser(option_list = option_list))

input <- opt$json
sample_id <- opt$sample

# ============================================================
# Helper functions
# ============================================================

# Extract the value associated with a particular level in native_path
get_path_value <- function(x, level_kind) {

  hit <- purrr::keep(
    x$native_path,
    ~ .x$level_kind == level_kind
  )

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[[1]]$name
}


# Convert a level_kind into a nice column name
# e.g. "supergroup" -> "Supergroup"
level_to_name <- function(x) {
  stringr::str_to_title(x)
}


# ============================================================
# Get classifier list
# ============================================================

parsed_json <- fromJSON(input, simplifyVector = FALSE)

classification_scores <- parsed_json$combined$classification_scores

classifier_list <- classification_scores %>%
  purrr::map_chr(~ .x$classifier_run_id) %>%
  unique()


# ============================================================
# Generic classifier
# ============================================================

make_generic_table <- function(scores, classifier_id) {

  results <- scores %>%
    purrr::keep(
      ~ .x$classifier_run_id == classifier_id
    )

  if (length(results) == 0) {
    stop("No results found for: ", classifier_id)
  }

  # --------------------------------------------------------
  # Determine which hierarchy levels are present
  # --------------------------------------------------------

  level_kinds <- results %>%
    purrr::map(
      ~ .x$native_path
    ) %>%
    purrr::flatten() %>%
    purrr::map_chr(~ .x$level_kind) %>%
    unique()


  # --------------------------------------------------------
  # Build one row for each classification
  # --------------------------------------------------------

  output <- results %>%
    purrr::map_dfr(function(x) {

      row <- list()

      # -----------------------------
      # Hierarchy
      # -----------------------------

      for (level in x$native_path) {

        level_name <- level_to_name(
          level$level_kind
        )

        row[[level_name]] <- level$name
      }


      # -----------------------------
      # Score
      #
      # For generic classifiers the
      # score corresponds to the leaf.
      # -----------------------------

      last_level <- x$native_path[[length(x$native_path)]]

      last_level_name <- level_to_name(
        last_level$level_kind
      )

      row[[paste0(last_level_name, "_score")]] <- x$score


      # -----------------------------
      # Cutoff
      # -----------------------------

      row$Cutoff <- ifelse(
        x$score >= x$cutoff,
        "PASS",
        "FAIL"
      )

      tibble::as_tibble(row)
    })


  # --------------------------------------------------------
  # Put columns in hierarchy order
  # --------------------------------------------------------

  hierarchy_columns <- level_kinds %>%
    purrr::map(level_to_name)

  score_columns <- paste0(
    hierarchy_columns,
    "_score"
  )

  # Keep only score for the leaf for generic classifiers
  # and arrange hierarchy before score.
  existing_hierarchy <- hierarchy_columns[
    hierarchy_columns %in% names(output)
  ]

  existing_scores <- score_columns[
    score_columns %in% names(output)
  ]

  output %>%
    select(
      all_of(unlist(hierarchy_columns)),
      all_of(unlist(existing_scores)),
      Cutoff
    )
}


# ============================================================
# MARLIN
# ============================================================

make_marlin_table <- function(df) {

  class_results <- df$combined$classification_scores %>%
    purrr::keep(
      ~ .x$classifier_run_id == "marlin"
    )


  # --------------------------------------------------------
  # Extract lineage probabilities
  # --------------------------------------------------------

  lineage_scores <- df$inference$lineage_probabilities %>%
    purrr::map_dfr(function(x) {

      tibble(
        Lineage = x$label,
        Lineage_score = as.numeric(x$probability)
      )
    })


  # --------------------------------------------------------
  # Extract group probabilities
  # --------------------------------------------------------

  group_scores <- df$inference$group_probabilities %>%
    purrr::map_dfr(function(x) {

      tibble(
        Group = x$label,
        Group_score = as.numeric(x$probability)
      )
    })


  # --------------------------------------------------------
  # Extract class/top probabilities
  # --------------------------------------------------------

  output <- class_results %>%
    purrr::map_dfr(function(x) {

      tibble(
        Lineage = get_path_value(
          x,
          "lineage"
        ),

        Group = get_path_value(
          x,
          "mcf"
        ),

        Class = x$leaf_name,

        Class_score = as.numeric(
          x$score
        ),

        Class_cutoff = as.numeric(
          x$cutoff
        )
      )
    })


  # --------------------------------------------------------
  # Add independent lineage/group scores
  # --------------------------------------------------------

  output <- output %>%
    left_join(
      lineage_scores,
      by = "Lineage"
    ) %>%
    left_join(
      group_scores,
      by = "Group"
    )


  # --------------------------------------------------------
  # Final cutoff
  #
  # Here PASS/FAIL refers to the final class score.
  # --------------------------------------------------------

  output %>%
    mutate(
      Cutoff = if_else(
        Class_score >= Class_cutoff,
        "PASS",
        "FAIL"
      )
    ) %>%
    select(
      Lineage,
      Lineage_score,
      Group,
      Group_score,
      Class,
      Class_score,
      Cutoff
    ) %>%
    arrange(desc(Class_score)) %>%
    head(5)
}


# ============================================================
# TUCAN
# ============================================================

make_tucan_table <- function(df) {

  Class_cutoff <- 0.7

  sampling <- df$tucan$samplings[[1]]

  family_scores <- sampling[["family_scores"]] %>%
    purrr::map_dfr(function(x) {

      tibble(
        Family = x$family,
        Family_score = as.numeric(x$score))
    })


  # --------------------------------------------------------
  # Connect class to family
  #
  # If each Tucan class result has a family/path field,
  # use that to create the relationship.
  # --------------------------------------------------------

  class_results <- sampling[["scores"]] %>%
    purrr::map_dfr(function(x) {

      family <- NA_character_

      if (!is.null(x$native_path)) {
        family <- get_path_value(
          x,
          "family"
        )
      }

      tibble(
        Family = x$family,
        Class = x$name,
        Class_score = as.numeric(x$score)
      )
    })


  output <- class_results %>%
    left_join(
      family_scores,
      by = "Family"
    )


  output %>%
    mutate(
      Cutoff = if_else(
        Class_score >= Class_cutoff,
        "PASS",
        "FAIL"
      )
    ) %>%
    select(
      Family,
      Family_score,
      Class,
      Class_score,
      Cutoff
    ) %>%
    arrange(desc(Class_score)) %>%
    head(5)
}


# ============================================================
# Generate all classifier tables
# ============================================================

classifier_tables <- list()

for (classifier_id in classifier_list) {

  message("Processing: ", classifier_id)

  if (classifier_id == "marlin") {

    classifier_tables[[classifier_id]] <-
      make_marlin_table(
        parsed_json
      )

  } else if (classifier_id == "tucan") {

    classifier_tables[[classifier_id]] <-
      make_tucan_table(
        parsed_json
      )

  } else {

    classifier_tables[[classifier_id]] <-
      make_generic_table(
        classification_scores,
        classifier_id
      )
  }
}


# ============================================================
# Write TSVs
# ============================================================

purrr::iwalk(
  classifier_tables,
  function(table, classifier_id) {

    filename <- classifier_id %>%
      stringr::str_replace_all(":", "_") %>%
      paste0(".tsv")

    filename_withsample <- paste(sample_id, filename, sep = "_")

    readr::write_tsv(
      table,
      filename_withsample
    )
  }
)
