library(SAHA)
library(tidyverse)
library(clustree)

# ==============================================================================
# 3_ SEA-AD ExN / InN / NonNeuronal consensus annotation + resolution QC
#
# Purpose:
#   - Process ExN, InN, and NonNeuronal together.
#   - Use the EC SAHA AutoAnnotate best_match when marker-based and marker-free
#     agree.
#   - ExN / InN: use the EC SAHA AutoAnnotate best_match when marker-based
#     and marker-free agree; otherwise use the EC marker-free call.
#   - NonNeuronal: use the independent SNACCS marker-free-only AutoAnnotate
#     result from script 2. No EC AutoAnnotate files are required for NN.
#   - Record ExN/InN marker-free fallbacks separately for each resolution.
#   - Record weak marker-free associations separately for each resolution.
#   - Add the selected annotation to the original SEA-AD resolution-assignment
#     taxonomy as leiden_<resolution>_anno_coarse.
#   - Save clustrees and resolution-level QC summaries for manual resolution
#     selection.
#
# This script does NOT choose the final resolution and does NOT create the final
# "anno" column. That is reserved for script 4 after anno_coarse is reviewed.
# ==============================================================================

# ==============================================================================
# 0. Paths and settings
# ==============================================================================

input_dir <- "/data/ADRD/brain_aging/phase2/public/seaad/resolution_selection"

resolution_root <- "/data/ADRD/brain_aging/exploration/res/Phase2/SEAAD_EC/resolution_selection"
ec_saha_dir <- file.path(resolution_root, "SAHA")
snaccs_saha_dir <- file.path(resolution_root, "SNACCS_NN_markerfree", "SAHA")

output_root <- file.path(resolution_root, "consensus_resolution_selection")
table_root <- file.path(output_root, "tables")
fallback_root <- file.path(table_root, "markerfree_fallbacks_by_resolution")
weak_root <- file.path(table_root, "weak_associations_by_resolution")

plot_root <- "/data/ADRD/brain_aging/exploration/plots/Phase2/SEAAD_EC/resolution_selection/consensus_resolution_selection"

cell_types <- c("ExN", "InN", "NonNeuronal")
resolution_values <- seq(0.1, 0.9, by = 0.1)
resolution_strings <- sprintf("%.1f", resolution_values)

# Diagnostic thresholds only.
# These thresholds DO NOT determine whether EC AutoAnnotate falls back to a
# marker-free call. They only flag weak selected associations for review.
weak_correlation_threshold <- 0.50
weak_margin_threshold <- 0.05

# NonNeuronal uses the independent SNACCS marker-free-only analysis from
# script 2. Its SAHA outputs are stored separately from the EC SAHA outputs.

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fallback_root, recursive = TRUE, showWarnings = FALSE)
dir.create(weak_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)

# These nested lists remain available when the script is sourced interactively.
tax_by_celltype <- list()
decision_by_celltype <- list()
markerfree_fallback_by_celltype <- list()
weak_association_by_celltype <- list()
resolution_summary_by_celltype <- list()

all_decisions_list <- list()
all_fallbacks_list <- list()
all_weak_list <- list()
all_resolution_summary_list <- list()

# ==============================================================================
# 1. Process each SEA-AD broad compartment
# ==============================================================================

for (cell_type in cell_types) {
  
  message("\n============================================================")
  message("Processing: ", cell_type)
  message("============================================================")
  
  tax_file <- file.path(
    input_dir,
    paste0("seaad_ec_", cell_type, "_resolution_assignments_obs.csv")
  )
  
  if (!file.exists(tax_file)) {
    stop("Missing SEA-AD resolution assignment table: ", tax_file)
  }
  
  tax <- read.csv(tax_file, check.names = FALSE)
  
  # read.csv(..., check.names = FALSE) can preserve blank header fields.
  # dplyr verbs reject data frames with NA or empty column names, so repair
  # only invalid/duplicate names while leaving all meaningful names unchanged.
  tax_names <- colnames(tax)
  bad_tax_names <- which(is.na(tax_names) | trimws(tax_names) == "")
  
  if (length(bad_tax_names) > 0) {
    tax_names[bad_tax_names] <- paste0("unnamed_", bad_tax_names)
  }
  
  colnames(tax) <- make.unique(tax_names, sep = "_")
  
  if (nrow(tax) == 0) {
    stop("SEA-AD resolution assignment table is empty: ", tax_file)
  }
  
  celltype_table_dir <- file.path(table_root, cell_type)
  celltype_fallback_dir <- file.path(fallback_root, cell_type)
  celltype_weak_dir <- file.path(weak_root, cell_type)
  celltype_plot_dir <- file.path(plot_root, cell_type)
  
  dir.create(celltype_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(celltype_fallback_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(celltype_weak_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(celltype_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # Resolve the original Leiden assignment column for every resolution.
  # ---------------------------------------------------------------------------
  
  leiden_cols <- setNames(rep(NA_character_, length(resolution_strings)), resolution_strings)
  
  for (res_i in seq_along(resolution_strings)) {
    res_str <- resolution_strings[res_i]
    
    candidates <- c(
      paste0("leiden_", res_str),
      paste0("leiden_scvi_", res_str),
      paste0("leiden_res", res_str),
      paste0("leiden_res_", res_str)
    )
    
    hits <- candidates[candidates %in% colnames(tax)]
    
    if (length(hits) == 0) {
      res_pattern <- gsub("\\.", "\\\\.", res_str)
      hits <- grep(
        paste0("^leiden.*", res_pattern, "$"),
        colnames(tax),
        value = TRUE
      )
    }
    
    if (length(hits) != 1) {
      stop(
        "Could not uniquely resolve the Leiden column for ",
        cell_type,
        " res ",
        res_str,
        ". Matches: ",
        paste(hits, collapse = ", ")
      )
    }
    
    leiden_cols[[res_str]] <- hits[1]
  }
  
  # ---------------------------------------------------------------------------
  # Resolve a donor column when present. Used only for resolution diagnostics.
  # ---------------------------------------------------------------------------
  
  donor_candidates <- c(
    "donor_id",
    "donor",
    "Donor ID",
    "donor_label",
    "individual_id",
    "individualID",
    "Individual ID",
    "individual"
  )
  
  donor_hits <- donor_candidates[donor_candidates %in% colnames(tax)]
  
  if (length(donor_hits) > 0) {
    donor_col <- donor_hits[1]
    message("Donor column: ", donor_col)
  } else {
    donor_col <- NA_character_
    message("No donor column detected for ", cell_type, ".")
  }
  
  # ---------------------------------------------------------------------------
  # Per-resolution objects for this compartment.
  # ---------------------------------------------------------------------------
  
  decision_by_resolution <- vector("list", length(resolution_values))
  markerfree_fallback_by_resolution <- vector("list", length(resolution_values))
  weak_association_by_resolution <- vector("list", length(resolution_values))
  names(decision_by_resolution) <- resolution_strings
  names(markerfree_fallback_by_resolution) <- resolution_strings
  names(weak_association_by_resolution) <- resolution_strings
  
  # =============================================================================
  # 2. Build annotation decision and QC tables for every resolution
  # =============================================================================
  
  for (res_i in seq_along(resolution_values)) {
    res <- resolution_values[res_i]
    res_str <- resolution_strings[res_i]
    leiden_col <- leiden_cols[[res_str]]
    
    message("\n---------- ", cell_type, " res ", res_str, " ----------")
    
    # -------------------------------------------------------------------------
    # Read the annotation source appropriate for this broad compartment.
    #
    # ExN / InN:
    #   EC atlas SAHA run with data_type = "Both".
    #
    # NonNeuronal:
    #   SNACCS marker-free-only SAHA run from script 2. These files live in:
    #   .../resolution_selection/SNACCS_NN_markerfree/SAHA
    # -------------------------------------------------------------------------
    
    ec_auto <- NULL
    ec_corr_summary <- NULL
    snaccs_auto <- NULL
    snaccs_corr_summary <- NULL
    
    if (cell_type %in% c("ExN", "InN")) {
      ec_auto_file <- file.path(
        ec_saha_dir,
        paste0("seaad_ec_", cell_type, "_res", res_str, "_SAHAvsECatlas_auto.csv")
      )
      
      ec_ann_file <- file.path(
        ec_saha_dir,
        paste0("seaad_ec_", cell_type, "_res", res_str, "_SAHAvsECatlas_ANN.RData")
      )
      
      if (!file.exists(ec_auto_file)) {
        stop("Missing EC AutoAnnotate file: ", ec_auto_file)
      }
      
      if (!file.exists(ec_ann_file)) {
        stop("Missing EC SAHA ANN file: ", ec_ann_file)
      }
      
      ec_auto <- read.csv(ec_auto_file, check.names = FALSE)
      
      ec_auto_names <- colnames(ec_auto)
      bad_ec_auto_names <- which(is.na(ec_auto_names) | trimws(ec_auto_names) == "")
      
      if (length(bad_ec_auto_names) > 0) {
        ec_auto_names[bad_ec_auto_names] <- paste0("unnamed_", bad_ec_auto_names)
      }
      
      colnames(ec_auto) <- make.unique(ec_auto_names, sep = "_")
      
      required_ec_auto_cols <- c(
        "cluster",
        "marker_based",
        "marker_free",
        "consensus",
        "best_match"
      )
      
      missing_ec_auto_cols <- setdiff(required_ec_auto_cols, colnames(ec_auto))
      
      if (length(missing_ec_auto_cols) > 0) {
        stop(
          "Unexpected EC AutoAnnotate format for ",
          cell_type,
          " res ",
          res_str,
          ". Missing: ",
          paste(missing_ec_auto_cols, collapse = ", ")
        )
      }
      
      ec_auto <- ec_auto %>%
        transmute(
          cluster = as.character(cluster),
          cluster = str_remove(cluster, "^query\\."),
          cluster = str_remove(cluster, "^RNA\\.g"),
          cluster = str_remove(cluster, "^RNA\\."),
          cluster = str_remove(cluster, "^g(?=[0-9]+$)"),
          ec_marker_based = as.character(marker_based),
          ec_marker_free = as.character(marker_free),
          ec_consensus = as.character(consensus),
          ec_best_match = as.character(best_match)
        )
      
      if (anyDuplicated(ec_auto$cluster) > 0) {
        duplicate_clusters <- unique(ec_auto$cluster[duplicated(ec_auto$cluster)])
        stop(
          "Duplicate EC AutoAnnotate clusters for ",
          cell_type,
          " res ",
          res_str,
          ": ",
          paste(duplicate_clusters, collapse = ", ")
        )
      }
      
      ec_env <- new.env()
      ec_loaded_objects <- load(ec_ann_file, envir = ec_env)
      
      if ("ann" %in% ec_loaded_objects) {
        ec_ann <- ec_env$ann
      } else if (length(ec_loaded_objects) == 1) {
        ec_ann <- ec_env[[ec_loaded_objects[1]]]
      } else {
        stop(
          "Could not uniquely identify the EC SAHA object in ",
          ec_ann_file,
          ". Objects: ",
          paste(ec_loaded_objects, collapse = ", ")
        )
      }
      
      ec_corr <- as.matrix(ec_ann@results$marker_free$corr)
      
      ec_corr_names <- colnames(ec_corr)
      bad_ec_corr_names <- which(is.na(ec_corr_names) | trimws(ec_corr_names) == "")
      
      if (length(bad_ec_corr_names) > 0) {
        ec_corr_names[bad_ec_corr_names] <- paste0("unnamed_reference_", bad_ec_corr_names)
      }
      
      colnames(ec_corr) <- make.unique(ec_corr_names, sep = "_")
      
      if (nrow(ec_corr) == 0 || ncol(ec_corr) == 0) {
        stop("Empty EC marker-free correlation matrix for ", cell_type, " res ", res_str)
      }
      
      ec_corr_long <- as.data.frame(ec_corr, check.names = FALSE) %>%
        rownames_to_column("cluster") %>%
        pivot_longer(
          cols = -cluster,
          names_to = "reference_raw",
          values_to = "correlation"
        ) %>%
        mutate(
          cluster = as.character(cluster),
          cluster = str_remove(cluster, "^query\\."),
          cluster = str_remove(cluster, "^RNA\\.g"),
          cluster = str_remove(cluster, "^RNA\\."),
          cluster = str_remove(cluster, "^g(?=[0-9]+$)"),
          reference = as.character(reference_raw),
          reference = str_remove(reference, "^db\\."),
          reference = str_replace_all(reference, "\\.", " "),
          reference = str_squish(reference),
          correlation = as.numeric(correlation)
        ) %>%
        filter(!is.na(correlation)) %>%
        arrange(cluster, desc(correlation), reference)
      
      ec_corr_summary <- ec_corr_long %>%
        group_by(cluster) %>%
        summarize(
          ec_markerfree_best_from_corr = first(reference),
          ec_markerfree_correlation = first(correlation),
          ec_markerfree_second_match = nth(reference, 2, default = NA_character_),
          ec_markerfree_second_correlation = nth(correlation, 2, default = NA_real_),
          ec_markerfree_margin = ec_markerfree_correlation - ec_markerfree_second_correlation,
          .groups = "drop"
        )
    }
    
    if (cell_type == "NonNeuronal") {
      snaccs_auto_file <- file.path(
        snaccs_saha_dir,
        paste0("seaad_ec_NonNeuronal_res", res_str, "_SAHAvsSNACCS_markerfree_auto.csv")
      )
      
      snaccs_ann_file <- file.path(
        snaccs_saha_dir,
        paste0("seaad_ec_NonNeuronal_res", res_str, "_SAHAvsSNACCS_markerfree_ANN.RData")
      )
      
      if (!file.exists(snaccs_auto_file)) {
        stop("Missing SNACCS marker-free AutoAnnotate file: ", snaccs_auto_file)
      }
      
      if (!file.exists(snaccs_ann_file)) {
        stop("Missing SNACCS marker-free ANN file: ", snaccs_ann_file)
      }
      
      snaccs_auto <- read.csv(snaccs_auto_file, check.names = FALSE)
      
      snaccs_auto_names <- colnames(snaccs_auto)
      bad_snaccs_auto_names <- which(is.na(snaccs_auto_names) | trimws(snaccs_auto_names) == "")
      
      if (length(bad_snaccs_auto_names) > 0) {
        snaccs_auto_names[bad_snaccs_auto_names] <- paste0("unnamed_", bad_snaccs_auto_names)
      }
      
      colnames(snaccs_auto) <- make.unique(snaccs_auto_names, sep = "_")
      
      required_snaccs_auto_cols <- c("cluster", "best_match", "correlation")
      missing_snaccs_auto_cols <- setdiff(required_snaccs_auto_cols, colnames(snaccs_auto))
      
      if (length(missing_snaccs_auto_cols) > 0) {
        stop(
          "Unexpected SNACCS marker-free format for NonNeuronal res ",
          res_str,
          ". Missing: ",
          paste(missing_snaccs_auto_cols, collapse = ", ")
        )
      }
      
      snaccs_auto <- snaccs_auto %>%
        transmute(
          cluster = as.character(cluster),
          cluster = str_remove(cluster, "^query\\."),
          cluster = str_remove(cluster, "^RNA\\.g"),
          cluster = str_remove(cluster, "^RNA\\."),
          cluster = str_remove(cluster, "^g(?=[0-9]+$)"),
          snaccs_marker_free = as.character(best_match),
          snaccs_auto_correlation = as.numeric(correlation)
        )
      
      if (anyDuplicated(snaccs_auto$cluster) > 0) {
        duplicate_clusters <- unique(snaccs_auto$cluster[duplicated(snaccs_auto$cluster)])
        stop(
          "Duplicate SNACCS marker-free clusters for NonNeuronal res ",
          res_str,
          ": ",
          paste(duplicate_clusters, collapse = ", ")
        )
      }
      
      snaccs_env <- new.env()
      snaccs_loaded_objects <- load(snaccs_ann_file, envir = snaccs_env)
      
      if ("ann" %in% snaccs_loaded_objects) {
        snaccs_ann <- snaccs_env$ann
      } else if (length(snaccs_loaded_objects) == 1) {
        snaccs_ann <- snaccs_env[[snaccs_loaded_objects[1]]]
      } else {
        stop(
          "Could not uniquely identify the SNACCS SAHA object in ",
          snaccs_ann_file,
          ". Objects: ",
          paste(snaccs_loaded_objects, collapse = ", ")
        )
      }
      
      snaccs_corr <- as.matrix(snaccs_ann@results$marker_free$corr)
      
      snaccs_corr_names <- colnames(snaccs_corr)
      bad_snaccs_corr_names <- which(is.na(snaccs_corr_names) | trimws(snaccs_corr_names) == "")
      
      if (length(bad_snaccs_corr_names) > 0) {
        snaccs_corr_names[bad_snaccs_corr_names] <- paste0("unnamed_reference_", bad_snaccs_corr_names)
      }
      
      colnames(snaccs_corr) <- make.unique(snaccs_corr_names, sep = "_")
      
      if (nrow(snaccs_corr) == 0 || ncol(snaccs_corr) == 0) {
        stop("Empty SNACCS marker-free correlation matrix for NonNeuronal res ", res_str)
      }
      
      snaccs_corr_long <- as.data.frame(snaccs_corr, check.names = FALSE) %>%
        rownames_to_column("cluster") %>%
        pivot_longer(
          cols = -cluster,
          names_to = "reference_raw",
          values_to = "correlation"
        ) %>%
        mutate(
          cluster = as.character(cluster),
          cluster = str_remove(cluster, "^query\\."),
          cluster = str_remove(cluster, "^RNA\\.g"),
          cluster = str_remove(cluster, "^RNA\\."),
          cluster = str_remove(cluster, "^g(?=[0-9]+$)"),
          reference = as.character(reference_raw),
          reference = str_remove(reference, "^db\\."),
          reference = str_replace_all(reference, "\\.", " "),
          reference = str_squish(reference),
          correlation = as.numeric(correlation)
        ) %>%
        filter(!is.na(correlation)) %>%
        arrange(cluster, desc(correlation), reference)
      
      snaccs_corr_summary <- snaccs_corr_long %>%
        group_by(cluster) %>%
        summarize(
          snaccs_markerfree_best_from_corr = first(reference),
          snaccs_markerfree_correlation = first(correlation),
          snaccs_markerfree_second_match = nth(reference, 2, default = NA_character_),
          snaccs_markerfree_second_correlation = nth(correlation, 2, default = NA_real_),
          snaccs_markerfree_margin = snaccs_markerfree_correlation - snaccs_markerfree_second_correlation,
          .groups = "drop"
        )
    }
    
    # -------------------------------------------------------------------------
    # Cell and donor representation of each query cluster.
    # -------------------------------------------------------------------------
    
    cluster_counts <- tax %>%
      transmute(
        cluster = as.character(.data[[leiden_col]]),
        cluster = str_remove(cluster, "^query\\."),
        cluster = str_remove(cluster, "^RNA\\.g"),
        cluster = str_remove(cluster, "^RNA\\."),
        cluster = str_remove(cluster, "^g(?=[0-9]+$)")
      ) %>%
      filter(!is.na(cluster)) %>%
      count(cluster, name = "n_cells")
    
    if (!is.na(donor_col)) {
      donor_counts <- tax %>%
        transmute(
          cluster = as.character(.data[[leiden_col]]),
          cluster = str_remove(cluster, "^query\\."),
          cluster = str_remove(cluster, "^RNA\\.g"),
          cluster = str_remove(cluster, "^RNA\\."),
          cluster = str_remove(cluster, "^g(?=[0-9]+$)"),
          donor = .data[[donor_col]]
        ) %>%
        filter(!is.na(cluster)) %>%
        group_by(cluster) %>%
        summarize(
          n_donors = n_distinct(donor, na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      donor_counts <- cluster_counts %>%
        transmute(
          cluster,
          n_donors = NA_integer_
        )
    }
    
    # -------------------------------------------------------------------------
    # Verify every metadata cluster can be annotated.
    # -------------------------------------------------------------------------
    
    metadata_clusters <- cluster_counts$cluster
    
    if (cell_type %in% c("ExN", "InN")) {
      missing_ec_auto_clusters <- setdiff(metadata_clusters, ec_auto$cluster)
      missing_ec_corr_clusters <- setdiff(metadata_clusters, ec_corr_summary$cluster)
      
      if (length(missing_ec_auto_clusters) > 0) {
        stop(
          "Missing EC AutoAnnotate results for ",
          cell_type,
          " res ",
          res_str,
          ": ",
          paste(missing_ec_auto_clusters, collapse = ", ")
        )
      }
      
      if (length(missing_ec_corr_clusters) > 0) {
        stop(
          "Missing EC marker-free correlations for ",
          cell_type,
          " res ",
          res_str,
          ": ",
          paste(missing_ec_corr_clusters, collapse = ", ")
        )
      }
    }
    
    if (cell_type == "NonNeuronal") {
      missing_snaccs_auto_clusters <- setdiff(metadata_clusters, snaccs_auto$cluster)
      missing_snaccs_corr_clusters <- setdiff(metadata_clusters, snaccs_corr_summary$cluster)
      
      if (length(missing_snaccs_auto_clusters) > 0) {
        stop(
          "Missing SNACCS marker-free results for NonNeuronal res ",
          res_str,
          ": ",
          paste(missing_snaccs_auto_clusters, collapse = ", ")
        )
      }
      
      if (length(missing_snaccs_corr_clusters) > 0) {
        stop(
          "Missing SNACCS marker-free correlations for NonNeuronal res ",
          res_str,
          ": ",
          paste(missing_snaccs_corr_clusters, collapse = ", ")
        )
      }
    }
    
    # -------------------------------------------------------------------------
    # Join the relevant result tables.
    # -------------------------------------------------------------------------
    
    if (cell_type %in% c("ExN", "InN")) {
      cluster_decision <- cluster_counts %>%
        left_join(donor_counts, by = "cluster") %>%
        left_join(ec_auto, by = "cluster") %>%
        left_join(ec_corr_summary, by = "cluster") %>%
        mutate(
          snaccs_marker_free = NA_character_,
          snaccs_auto_correlation = NA_real_,
          snaccs_markerfree_best_from_corr = NA_character_,
          snaccs_markerfree_correlation = NA_real_,
          snaccs_markerfree_second_match = NA_character_,
          snaccs_markerfree_second_correlation = NA_real_,
          snaccs_markerfree_margin = NA_real_,
          snaccs_markerfree_matches_corr = NA
        )
      
      cluster_decision <- cluster_decision %>%
        mutate(
          cell_type = cell_type,
          resolution = res,
          ec_consensus_upper = str_to_upper(str_squish(ec_consensus)),
          ec_best_match_upper = str_to_upper(str_squish(ec_best_match)),
          use_marker_free = is.na(ec_consensus_upper) |
            ec_consensus_upper != "MATCH" |
            is.na(ec_best_match_upper) |
            ec_best_match_upper %in% c("INCONCLUSIVE", "UNKNOWN", ""),
          markerfree_fallback = use_marker_free,
          fallback_reason = case_when(
            !markerfree_fallback ~ NA_character_,
            is.na(ec_consensus_upper) ~ "EC_CONSENSUS_MISSING",
            ec_consensus_upper == "DISAGREEMENT" ~ "EC_DISAGREEMENT",
            is.na(ec_best_match_upper) ~ "EC_BEST_MATCH_MISSING",
            ec_best_match_upper %in% c("INCONCLUSIVE", "UNKNOWN", "") ~ "EC_INCONCLUSIVE",
            TRUE ~ paste0("EC_", ec_consensus_upper)
          ),
          anno_coarse = if_else(markerfree_fallback, ec_marker_free, ec_best_match),
          anno_source = if_else(
            markerfree_fallback,
            "EC_marker_free_fallback",
            "EC_auto_match"
          ),
          selected_correlation = ec_markerfree_correlation,
          selected_second_match = ec_markerfree_second_match,
          selected_second_correlation = ec_markerfree_second_correlation,
          selected_margin = ec_markerfree_margin,
          ec_markerfree_matches_corr = str_remove_all(str_to_lower(ec_marker_free), "[^a-z0-9]") ==
            str_remove_all(str_to_lower(ec_markerfree_best_from_corr), "[^a-z0-9]")
        )
    }
    
    if (cell_type == "NonNeuronal") {
      cluster_decision <- cluster_counts %>%
        left_join(donor_counts, by = "cluster") %>%
        left_join(snaccs_auto, by = "cluster") %>%
        left_join(snaccs_corr_summary, by = "cluster") %>%
        mutate(
          cell_type = cell_type,
          resolution = res,
          ec_marker_based = NA_character_,
          ec_marker_free = NA_character_,
          ec_consensus = NA_character_,
          ec_best_match = NA_character_,
          ec_markerfree_best_from_corr = NA_character_,
          ec_markerfree_correlation = NA_real_,
          ec_markerfree_second_match = NA_character_,
          ec_markerfree_second_correlation = NA_real_,
          ec_markerfree_margin = NA_real_,
          ec_markerfree_matches_corr = NA,
          use_marker_free = TRUE,
          markerfree_fallback = FALSE,
          fallback_reason = NA_character_,
          anno_coarse = snaccs_marker_free,
          anno_source = "SNACCS_marker_free_only",
          selected_correlation = snaccs_markerfree_correlation,
          selected_second_match = snaccs_markerfree_second_match,
          selected_second_correlation = snaccs_markerfree_second_correlation,
          selected_margin = snaccs_markerfree_margin,
          snaccs_markerfree_matches_corr = str_remove_all(str_to_lower(snaccs_marker_free), "[^a-z0-9]") ==
            str_remove_all(str_to_lower(snaccs_markerfree_best_from_corr), "[^a-z0-9]")
        )
    }
    
    cluster_decision <- cluster_decision %>%
      mutate(
        weak_correlation = is.na(selected_correlation) |
          selected_correlation < weak_correlation_threshold,
        weak_margin = is.na(selected_margin) |
          selected_margin < weak_margin_threshold,
        weak_association = weak_correlation | weak_margin,
        weak_reason = case_when(
          is.na(selected_correlation) & is.na(selected_margin) ~ "MISSING_CORRELATION_AND_MARGIN",
          is.na(selected_correlation) ~ "MISSING_SELECTED_CORRELATION",
          is.na(selected_margin) ~ "MISSING_TOP2_MARGIN",
          selected_correlation < weak_correlation_threshold & selected_margin < weak_margin_threshold ~ "LOW_CORRELATION_AND_SMALL_MARGIN",
          selected_correlation < weak_correlation_threshold ~ "LOW_CORRELATION",
          selected_margin < weak_margin_threshold ~ "SMALL_TOP2_MARGIN",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::select(
        cell_type,
        resolution,
        cluster,
        n_cells,
        n_donors,
        ec_marker_based,
        ec_marker_free,
        ec_consensus,
        ec_best_match,
        ec_markerfree_best_from_corr,
        ec_markerfree_correlation,
        ec_markerfree_second_match,
        ec_markerfree_second_correlation,
        ec_markerfree_margin,
        ec_markerfree_matches_corr,
        snaccs_marker_free,
        snaccs_auto_correlation,
        snaccs_markerfree_best_from_corr,
        snaccs_markerfree_correlation,
        snaccs_markerfree_second_match,
        snaccs_markerfree_second_correlation,
        snaccs_markerfree_margin,
        snaccs_markerfree_matches_corr,
        use_marker_free,
        markerfree_fallback,
        fallback_reason,
        anno_coarse,
        anno_source,
        selected_correlation,
        selected_second_match,
        selected_second_correlation,
        selected_margin,
        weak_correlation,
        weak_margin,
        weak_association,
        weak_reason
      )
    
    if (any(is.na(cluster_decision$anno_coarse)) || any(cluster_decision$anno_coarse == "")) {
      bad_clusters <- cluster_decision %>%
        filter(is.na(anno_coarse) | anno_coarse == "") %>%
        pull(cluster)
      
      stop(
        "Missing final anno_coarse for ",
        cell_type,
        " res ",
        res_str,
        ": ",
        paste(bad_clusters, collapse = ", ")
      )
    }
    
    # -------------------------------------------------------------------------
    # Store and write the full decisions, marker-free fallbacks, and weak calls.
    # -------------------------------------------------------------------------
    
    markerfree_fallback_table <- cluster_decision %>%
      filter(.data$markerfree_fallback)
    
    weak_association_table <- cluster_decision %>%
      filter(.data$weak_association)
    
    decision_by_resolution[[res_str]] <- cluster_decision
    markerfree_fallback_by_resolution[[res_str]] <- markerfree_fallback_table
    weak_association_by_resolution[[res_str]] <- weak_association_table
    
    write.csv(
      cluster_decision,
      file.path(
        celltype_table_dir,
        paste0("seaad_ec_", cell_type, "_res", res_str, "_annotation_decisions.csv")
      ),
      row.names = FALSE
    )
    
    write.csv(
      markerfree_fallback_table,
      file.path(
        celltype_fallback_dir,
        paste0("seaad_ec_", cell_type, "_res", res_str, "_markerfree_fallbacks.csv")
      ),
      row.names = FALSE
    )
    
    write.csv(
      weak_association_table,
      file.path(
        celltype_weak_dir,
        paste0("seaad_ec_", cell_type, "_res", res_str, "_weak_associations.csv")
      ),
      row.names = FALSE
    )
    
    # -------------------------------------------------------------------------
    # Add selected calls and QC columns to the cell-level taxonomy.
    # -------------------------------------------------------------------------
    
    anno_col <- paste0("leiden_", res_str, "_anno_coarse")
    source_col <- paste0("leiden_", res_str, "_anno_source")
    fallback_col <- paste0("leiden_", res_str, "_markerfree_fallback")
    weak_col <- paste0("leiden_", res_str, "_weak_association")
    correlation_col <- paste0("leiden_", res_str, "_selected_correlation")
    margin_col <- paste0("leiden_", res_str, "_selected_margin")
    
    tax_lookup <- cluster_decision %>%
      transmute(
        cluster,
        .anno_coarse = anno_coarse,
        .anno_source = anno_source,
        .markerfree_fallback = .data$markerfree_fallback,
        .weak_association = .data$weak_association,
        .selected_correlation = selected_correlation,
        .selected_margin = selected_margin
      )
    
    tax <- tax %>%
      mutate(
        .cluster_join = as.character(.data[[leiden_col]]),
        .cluster_join = str_remove(.cluster_join, "^query\\."),
        .cluster_join = str_remove(.cluster_join, "^RNA\\.g"),
        .cluster_join = str_remove(.cluster_join, "^RNA\\."),
        .cluster_join = str_remove(.cluster_join, "^g(?=[0-9]+$)")
      ) %>%
      left_join(tax_lookup, by = c(".cluster_join" = "cluster")) %>%
      rename(
        !!anno_col := .anno_coarse,
        !!source_col := .anno_source,
        !!fallback_col := .markerfree_fallback,
        !!weak_col := .weak_association,
        !!correlation_col := .selected_correlation,
        !!margin_col := .selected_margin
      ) %>%
      dplyr::select(-.cluster_join)
    
    if (any(is.na(tax[[anno_col]]))) {
      stop("Missing cell-level anno_coarse after join for ", cell_type, " res ", res_str)
    }
  }
  
  # =============================================================================
  # 3. Combine per-resolution tables for this cell compartment
  # ==============================================================================
  
  decision_all <- bind_rows(decision_by_resolution)
  markerfree_fallback_all <- bind_rows(markerfree_fallback_by_resolution)
  weak_association_all <- bind_rows(weak_association_by_resolution)
  
  write.csv(
    decision_all,
    file.path(celltype_table_dir, paste0("seaad_ec_", cell_type, "_annotation_decisions_all_resolutions.csv")),
    row.names = FALSE
  )
  
  write.csv(
    markerfree_fallback_all,
    file.path(celltype_table_dir, paste0("seaad_ec_", cell_type, "_markerfree_fallbacks_all_resolutions.csv")),
    row.names = FALSE
  )
  
  write.csv(
    weak_association_all,
    file.path(celltype_table_dir, paste0("seaad_ec_", cell_type, "_weak_associations_all_resolutions.csv")),
    row.names = FALSE
  )
  
  write.csv(
    tax,
    file.path(celltype_table_dir, paste0("seaad_ec_", cell_type, "_resolution_assignments_anno_coarse.csv")),
    row.names = FALSE
  )
  
  # =============================================================================
  # 4. Resolution-selection summary
  # ==============================================================================
  
  resolution_summary <- decision_all %>%
    mutate(
      cluster_n_cells = n_cells,
      cluster_n_donors = n_donors
    ) %>%
    group_by(cell_type, resolution) %>%
    summarize(
      n_clusters = n(),
      n_cells = sum(cluster_n_cells),
      min_cluster_cells = min(cluster_n_cells),
      median_cluster_cells = median(cluster_n_cells),
      min_cluster_donors = ifelse(
        all(is.na(cluster_n_donors)),
        NA_real_,
        min(cluster_n_donors, na.rm = TRUE)
      ),
      median_cluster_donors = ifelse(
        all(is.na(cluster_n_donors)),
        NA_real_,
        median(cluster_n_donors, na.rm = TRUE)
      ),
      n_markerfree_fallback_clusters = sum(markerfree_fallback),
      fraction_markerfree_fallback_clusters = mean(markerfree_fallback),
      n_markerfree_fallback_cells = sum(cluster_n_cells[markerfree_fallback]),
      fraction_markerfree_fallback_cells = sum(cluster_n_cells[markerfree_fallback]) / sum(cluster_n_cells),
      n_weak_clusters = sum(weak_association),
      fraction_weak_clusters = mean(weak_association),
      n_weak_cells = sum(cluster_n_cells[weak_association]),
      fraction_weak_cells = sum(cluster_n_cells[weak_association]) / sum(cluster_n_cells),
      mean_selected_correlation_weighted = ifelse(
        all(is.na(selected_correlation)),
        NA_real_,
        sum(selected_correlation * cluster_n_cells, na.rm = TRUE) /
          sum(cluster_n_cells[!is.na(selected_correlation)])
      ),
      median_selected_correlation = ifelse(
        all(is.na(selected_correlation)),
        NA_real_,
        median(selected_correlation, na.rm = TRUE)
      ),
      min_selected_correlation = ifelse(
        all(is.na(selected_correlation)),
        NA_real_,
        min(selected_correlation, na.rm = TRUE)
      ),
      mean_selected_margin_weighted = ifelse(
        all(is.na(selected_margin)),
        NA_real_,
        sum(selected_margin * cluster_n_cells, na.rm = TRUE) /
          sum(cluster_n_cells[!is.na(selected_margin)])
      ),
      median_selected_margin = ifelse(
        all(is.na(selected_margin)),
        NA_real_,
        median(selected_margin, na.rm = TRUE)
      ),
      min_selected_margin = ifelse(
        all(is.na(selected_margin)),
        NA_real_,
        min(selected_margin, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  write.csv(
    resolution_summary,
    file.path(celltype_table_dir, paste0("seaad_ec_", cell_type, "_resolution_selection_summary.csv")),
    row.names = FALSE
  )
  
  # =============================================================================
  # 5. Clustree of selected anno_coarse calls
  # ==============================================================================
  
  anno_coarse_cols <- paste0("leiden_", resolution_strings, "_anno_coarse")
  
  missing_anno_coarse_cols <- setdiff(anno_coarse_cols, colnames(tax))
  
  if (length(missing_anno_coarse_cols) > 0) {
    stop(
      "Missing anno_coarse columns before clustree for ",
      cell_type,
      ": ",
      paste(missing_anno_coarse_cols, collapse = ", ")
    )
  }
  
  clustree_df <- tax %>%
    dplyr::select(all_of(anno_coarse_cols))
  
  png(
    file.path(celltype_plot_dir, paste0("seaad_ec_", cell_type, "_anno_coarse_clustree.png")),
    width = 20,
    height = 12,
    units = "in",
    res = 300
  )
  print(clustree(clustree_df, prefix = "leiden_", suffix = "_anno_coarse"))
  dev.off()
  
  pdf(
    file.path(celltype_plot_dir, paste0("seaad_ec_", cell_type, "_anno_coarse_clustree.pdf")),
    width = 20,
    height = 12
  )
  print(clustree(clustree_df, prefix = "leiden_", suffix = "_anno_coarse"))
  dev.off()
  
  # =============================================================================
  # 6. Resolution QC plots
  # ==============================================================================
  
  resolution_qc_long <- resolution_summary %>%
    dplyr::select(
      cell_type,
      resolution,
      fraction_markerfree_fallback_cells,
      fraction_weak_cells
    ) %>%
    pivot_longer(
      cols = c(fraction_markerfree_fallback_cells, fraction_weak_cells),
      names_to = "metric",
      values_to = "fraction_cells"
    ) %>%
    mutate(
      metric = recode(
        metric,
        fraction_markerfree_fallback_cells = "Marker-free fallback",
        fraction_weak_cells = "Weak selected association"
      )
    )
  
  p_resolution_qc <- ggplot(
    resolution_qc_long,
    aes(x = resolution, y = fraction_cells, group = metric, linetype = metric)
  ) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = resolution_values) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      x = "Leiden resolution",
      y = "Fraction of cells",
      linetype = NULL,
      title = paste0("SEA-AD ", cell_type, ": fallback and weak-association burden")
    ) +
    theme_bw()
  
  ggsave(
    file.path(celltype_plot_dir, paste0("seaad_ec_", cell_type, "_resolution_fallback_weak_fraction.pdf")),
    p_resolution_qc,
    width = 8,
    height = 5
  )
  
  p_correlation <- ggplot(
    decision_all,
    aes(x = factor(resolution), y = selected_correlation)
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(size = n_cells), width = 0.15, alpha = 0.6) +
    geom_hline(yintercept = weak_correlation_threshold, linetype = "dashed") +
    labs(
      x = "Leiden resolution",
      y = "Selected marker-free correlation",
      size = "Cluster cells",
      title = paste0("SEA-AD ", cell_type, ": selected annotation correlation")
    ) +
    theme_bw()
  
  ggsave(
    file.path(celltype_plot_dir, paste0("seaad_ec_", cell_type, "_selected_correlation_by_resolution.pdf")),
    p_correlation,
    width = 9,
    height = 6
  )
  
  p_margin <- ggplot(
    decision_all,
    aes(x = factor(resolution), y = selected_margin)
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(size = n_cells), width = 0.15, alpha = 0.6) +
    geom_hline(yintercept = weak_margin_threshold, linetype = "dashed") +
    labs(
      x = "Leiden resolution",
      y = "Best minus second-best marker-free correlation",
      size = "Cluster cells",
      title = paste0("SEA-AD ", cell_type, ": marker-free top-two margin")
    ) +
    theme_bw()
  
  ggsave(
    file.path(celltype_plot_dir, paste0("seaad_ec_", cell_type, "_selected_margin_by_resolution.pdf")),
    p_margin,
    width = 9,
    height = 6
  )
  
  # =============================================================================
  # 7. Retain this compartment in memory and in combined all-cell-type outputs
  # ==============================================================================
  
  tax_by_celltype[[cell_type]] <- tax
  decision_by_celltype[[cell_type]] <- decision_by_resolution
  markerfree_fallback_by_celltype[[cell_type]] <- markerfree_fallback_by_resolution
  weak_association_by_celltype[[cell_type]] <- weak_association_by_resolution
  resolution_summary_by_celltype[[cell_type]] <- resolution_summary
  
  all_decisions_list[[cell_type]] <- decision_all
  all_fallbacks_list[[cell_type]] <- markerfree_fallback_all
  all_weak_list[[cell_type]] <- weak_association_all
  all_resolution_summary_list[[cell_type]] <- resolution_summary
  
  message("Completed: ", cell_type)
}

# ==============================================================================
# 8. Combined ExN + InN + NonNeuronal outputs
# ==============================================================================

all_decisions <- bind_rows(all_decisions_list)
all_markerfree_fallbacks <- bind_rows(all_fallbacks_list)
all_weak_associations <- bind_rows(all_weak_list)
all_resolution_summary <- bind_rows(all_resolution_summary_list)

write.csv(
  all_decisions,
  file.path(table_root, "seaad_ec_all_celltypes_annotation_decisions_all_resolutions.csv"),
  row.names = FALSE
)

write.csv(
  all_markerfree_fallbacks,
  file.path(table_root, "seaad_ec_all_celltypes_markerfree_fallbacks_all_resolutions.csv"),
  row.names = FALSE
)

write.csv(
  all_weak_associations,
  file.path(table_root, "seaad_ec_all_celltypes_weak_associations_all_resolutions.csv"),
  row.names = FALSE
)

write.csv(
  all_resolution_summary,
  file.path(table_root, "seaad_ec_all_celltypes_resolution_selection_summary.csv"),
  row.names = FALSE
)

# Combined overview plot across all three broad compartments.
p_all_resolution_qc <- all_resolution_summary %>%
  dplyr::select(
    cell_type,
    resolution,
    fraction_markerfree_fallback_cells,
    fraction_weak_cells
  ) %>%
  pivot_longer(
    cols = c(fraction_markerfree_fallback_cells, fraction_weak_cells),
    names_to = "metric",
    values_to = "fraction_cells"
  ) %>%
  mutate(
    metric = recode(
      metric,
      fraction_markerfree_fallback_cells = "Marker-free fallback",
      fraction_weak_cells = "Weak selected association"
    )
  ) %>%
  ggplot(
    aes(x = resolution, y = fraction_cells, group = metric, linetype = metric)
  ) +
  geom_line() +
  geom_point() +
  facet_wrap(~cell_type, ncol = 1) +
  scale_x_continuous(breaks = resolution_values) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Leiden resolution",
    y = "Fraction of cells",
    linetype = NULL,
    title = "SEA-AD resolution QC across ExN, InN, and NonNeuronal"
  ) +
  theme_bw()

ggsave(
  file.path(plot_root, "seaad_ec_all_celltypes_resolution_fallback_weak_fraction.pdf"),
  p_all_resolution_qc,
  width = 9,
  height = 10
)

capture.output(
  sessionInfo(),
  file = file.path(output_root, "SEAAD_all_celltypes_consensus_resolution_selection_sessionInfo.txt")
)

message("\n============================================================")
message("Completed ExN, InN, and NonNeuronal consensus resolution selection.")
message("Results: ", output_root)
message("Plots:   ", plot_root)
message("============================================================")
