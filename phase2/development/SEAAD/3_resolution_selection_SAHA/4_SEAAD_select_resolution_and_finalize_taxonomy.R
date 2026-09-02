library(tidyverse)
library(clustree)

# ==============================================================================
# 4_ SEA-AD manual resolution selection + final taxonomy annotation
#
# Purpose:
#   1. Review the ExN / InN / NonNeuronal clustrees and resolution QC from
#      script 3.
#   2. Manually choose one final Leiden resolution for each broad compartment.
#   3. Review every cluster at the selected resolution, including the automatic
#      annotation from script 3 and its QC.
#   4. Manually assign the final cluster name in a new "anno" column.
#   5. Manually group those final names into "anno_coarse" using case_when().
#   6. Save final per-compartment and combined taxonomy tables.
#
# IMPORTANT:
#   - This script intentionally does NOT automatically choose a resolution.
#   - This script intentionally does NOT automatically choose the final "anno"
#     names.
#   - This script intentionally leaves "anno_coarse" for a manual case_when()
#     written from the final "anno" names.
#
# Recommended workflow:
#   A. First run with selected_resolution left as NA_character_. Review the
#      clustrees and QC summary produced/printed below.
#   B. Enter the selected resolution for ExN, InN, and NonNeuronal and rerun.
#      The script writes cluster-review and annotation-template files.
#   C. Fill ExN_anno_map, InN_anno_map, and NonNeuronal_anno_map below and
#      rerun. The script writes taxonomy tables containing the final "anno".
#   D. Write the three anno_coarse case_when() blocks near the bottom, set
#      finalize_anno_coarse <- TRUE, and rerun to write final taxonomy files.
# ==============================================================================

# ==============================================================================
# 0. Paths
# ==============================================================================

resolution_root <- "/data/ADRD/brain_aging/exploration/res/Phase2/SEAAD_EC/resolution_selection"
script3_root <- file.path(resolution_root, "consensus_resolution_selection")
script3_table_root <- file.path(script3_root, "tables")

script3_plot_root <- "/data/ADRD/brain_aging/exploration/plots/Phase2/SEAAD_EC/resolution_selection/consensus_resolution_selection"

output_root <- file.path(resolution_root, "final_annotation")
review_root <- file.path(output_root, "resolution_review")
final_table_root <- file.path(output_root, "tables")

plot_root <- "/data/ADRD/brain_aging/exploration/plots/Phase2/SEAAD_EC/resolution_selection/final_annotation"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(review_root, recursive = TRUE, showWarnings = FALSE)
dir.create(final_table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)

cell_types <- c("ExN", "InN", "NonNeuronal")
resolution_values <- seq(0.1, 0.9, by = 0.1)
resolution_strings <- sprintf("%.1f", resolution_values)

# ==============================================================================
# 1. USER EDIT: choose one final resolution for each broad compartment
# ==============================================================================
#
# Leave these as NA_character_ on the first run.
# Review the clustrees and resolution-selection summary, then enter one of:
# "0.1", "0.2", ..., "0.9"
#
# Example only:
selected_resolution <- c(
  ExN = "0.3",
  InN = "0.4",
  NonNeuronal = "0.2"
)

#for the first time running
# selected_resolution <- c(
#   ExN = NA_character_,
#   InN = NA_character_,
#   NonNeuronal = NA_character_
# )

# ==============================================================================
# 2. USER EDIT: final cluster-number -> anno mappings
# ==============================================================================
#
# After selecting resolutions and rerunning this script, look at the files:
#
#   resolution_review/ExN_selected_resX_cluster_review.csv
#   resolution_review/InN_selected_resX_cluster_review.csv
#   resolution_review/NonNeuronal_selected_resX_cluster_review.csv
#
# The script also writes *_anno_map_template.R files containing every actual
# cluster number at the selected resolution. Copy those lines here and replace
# each empty string with your desired FINAL annotation name.
#
# Example only:
ExN_anno_map <- c(
  "0"  = "ExN SEMA3E CCN2",
  "1"  = "ExN SEMA3E SULF1",
  "2"  = "ExN SEMA3E ITPR2",
  "3"  = "ExN SEMA3E SULF1",
  "4"  = "ExN RORB TLL1",
  "5"  = "ExN CUX2 GLP2R",
  "6"  = "ExN LAMP5 SPON1",
  "7"  = "ExN THEMIS BACE2",
  "8"  = "ExN CUX2 GLP2R",
  "9"  = "ExN CUX2 RELN",
  "10" = "ExN RELN ZNF536"
)

InN_anno_map <- c(
  "0"  = "InN SST PTPRK",
  "1"  = "InN PVALB MYO5B",
  "2"  = "InN SST PTPRK",
  "3"  = "InN SST PDE1A",
  "4"  = "InN SST PDE1A",
  "5"  = "InN VIP PTH2R",
  "6"  = "InN PAX6 PDLIM5",
  "7"  = "InN VIP PCDH18",
  "8"  = "InN VIP PCDH18",
  "9"  = "InN LAMP5 COL5A2",
  "10" = "InN LAMP5 KIT",
  "11" = "InN VIP TNS3",
  "12" = "InN VIP PCDH18",
  "13" = "InN PAX6 DDR2",
  "14" = "InN LAMP5 CHST9",
  "15" = "InN LAMP5 CHST9",
  "16" = "InN PVALB MYO5B",
  "17" = "InN PVALB MYO5B"
)

NonNeuronal_anno_map <- c(
  "0" = "OPC2",
  "1" = "Astro3",
  "2" = "Astro3",
  "3" = "Micro2",
  "4" = "Oligo3",
  "5" = "BEC venous",
  "6" = "Micro1",
  "7" = "Oligo2",
  "8" = "Peri"
)

# #as place holder
# ExN_anno_map <- character(0)
# 
# InN_anno_map <- character(0)
# 
# NonNeuronal_anno_map <- character(0)

# ==============================================================================
# 3. USER EDIT LAST: switch this on only after writing anno_coarse case_when()
# ==============================================================================

finalize_anno_coarse <- TRUE

# ==============================================================================
# 4. Load script 3 taxonomies and resolution-selection diagnostics
# ==============================================================================

tax_by_celltype <- list()
selected_leiden_col_by_celltype <- list()
cluster_review_by_celltype <- list()
final_tax_by_celltype <- list()

resolution_summary_file <- file.path(
  script3_table_root,
  "seaad_ec_all_celltypes_resolution_selection_summary.csv"
)

if (!file.exists(resolution_summary_file)) {
  stop("Missing script 3 resolution-selection summary: ", resolution_summary_file)
}

resolution_summary <- read.csv(resolution_summary_file, check.names = FALSE)

resolution_summary_names <- colnames(resolution_summary)
bad_resolution_summary_names <- which(
  is.na(resolution_summary_names) | trimws(resolution_summary_names) == ""
)

if (length(bad_resolution_summary_names) > 0) {
  resolution_summary_names[bad_resolution_summary_names] <- paste0(
    "unnamed_",
    bad_resolution_summary_names
  )
}

colnames(resolution_summary) <- make.unique(resolution_summary_names, sep = "_")

message("\n============================================================")
message("RESOLUTION-SELECTION SUMMARY FROM SCRIPT 3")
message("============================================================")
print(as.data.frame(resolution_summary))

write.csv(
  resolution_summary,
  file.path(review_root, "seaad_ec_all_celltypes_resolution_selection_summary_for_review.csv"),
  row.names = FALSE
)

for (cell_type in cell_types) {
  tax_file <- file.path(
    script3_table_root,
    cell_type,
    paste0("seaad_ec_", cell_type, "_resolution_assignments_anno_coarse.csv")
  )
  
  if (!file.exists(tax_file)) {
    stop("Missing script 3 taxonomy for ", cell_type, ": ", tax_file)
  }
  
  tax <- read.csv(tax_file, check.names = FALSE)
  
  tax_names <- colnames(tax)
  bad_tax_names <- which(is.na(tax_names) | trimws(tax_names) == "")
  
  if (length(bad_tax_names) > 0) {
    tax_names[bad_tax_names] <- paste0("unnamed_", bad_tax_names)
  }
  
  colnames(tax) <- make.unique(tax_names, sep = "_")
  
  if (nrow(tax) == 0) {
    stop("Script 3 taxonomy is empty for ", cell_type, ": ", tax_file)
  }
  
  tax_by_celltype[[cell_type]] <- tax
}

# ==============================================================================
# 5. Recreate all three clustrees in one review PDF
# ==============================================================================
#
# These use the script 3 automatic anno_coarse calls across resolutions.
# The original individual PNG/PDF files from script 3 remain untouched.
# ==============================================================================

clustree_review_pdf <- file.path(
  plot_root,
  "seaad_ec_ExN_InN_NonNeuronal_anno_coarse_clustrees_for_resolution_selection.pdf"
)

pdf(
  clustree_review_pdf,
  width = 20,
  height = 12,
  onefile = TRUE
)

for (cell_type in cell_types) {
  tax <- tax_by_celltype[[cell_type]]
  anno_coarse_cols <- paste0("leiden_", resolution_strings, "_anno_coarse")
  missing_anno_coarse_cols <- setdiff(anno_coarse_cols, colnames(tax))
  
  if (length(missing_anno_coarse_cols) > 0) {
    dev.off()
    stop(
      "Missing script 3 anno_coarse columns for ",
      cell_type,
      ": ",
      paste(missing_anno_coarse_cols, collapse = ", ")
    )
  }
  
  clustree_df <- tax %>%
    dplyr::select(all_of(anno_coarse_cols))
  
  print(
    clustree(
      clustree_df,
      prefix = "leiden_",
      suffix = "_anno_coarse"
    ) +
      ggtitle(paste0("SEA-AD ", cell_type, ": automatic annotation clustree"))
  )
}

dev.off()

message("\nClustree review PDF:")
message(clustree_review_pdf)

message("\nOriginal script 3 clustrees:")
for (cell_type in cell_types) {
  message(
    "  ",
    cell_type,
    ": ",
    file.path(
      script3_plot_root,
      cell_type,
      paste0("seaad_ec_", cell_type, "_anno_coarse_clustree.png")
    )
  )
}

# ==============================================================================
# 6. Validate whether final resolutions have been selected
# ==============================================================================

if (!setequal(names(selected_resolution), cell_types)) {
  stop(
    "selected_resolution must be a named vector containing exactly: ",
    paste(cell_types, collapse = ", ")
  )
}

selected_resolution <- selected_resolution[cell_types]

resolution_selection_complete <- all(
  !is.na(selected_resolution) &
    selected_resolution %in% resolution_strings
)

if (!resolution_selection_complete) {
  message("\n============================================================")
  message("RESOLUTION SELECTION NOT YET COMPLETE")
  message("============================================================")
  message("Review the clustree PDF and resolution-selection summary above.")
  message("Then edit selected_resolution near the top of this script and rerun.")
  message("No final resolution was chosen automatically.")
}

# ==============================================================================
# 7. After resolution selection: build cluster-review and mapping templates
# ==============================================================================

mapping_ready <- setNames(rep(FALSE, length(cell_types)), cell_types)
anno_map_by_celltype <- list(
  ExN = ExN_anno_map,
  InN = InN_anno_map,
  NonNeuronal = NonNeuronal_anno_map
)

if (resolution_selection_complete) {
  selected_resolution_table <- tibble(
    cell_type = cell_types,
    selected_resolution = unname(selected_resolution[cell_types])
  )
  
  write.csv(
    selected_resolution_table,
    file.path(review_root, "seaad_ec_selected_resolutions.csv"),
    row.names = FALSE
  )
  
  message("\n============================================================")
  message("SELECTED RESOLUTIONS")
  message("============================================================")
  print(as.data.frame(selected_resolution_table))
  
  for (cell_type in cell_types) {
    tax <- tax_by_celltype[[cell_type]]
    res_str <- selected_resolution[[cell_type]]
    
    # -------------------------------------------------------------------------
    # Resolve the original Leiden cluster column at the selected resolution.
    # -------------------------------------------------------------------------
    
    leiden_candidates <- c(
      paste0("leiden_", res_str),
      paste0("leiden_scvi_", res_str),
      paste0("leiden_res", res_str),
      paste0("leiden_res_", res_str)
    )
    
    leiden_hits <- leiden_candidates[leiden_candidates %in% colnames(tax)]
    
    if (length(leiden_hits) == 0) {
      res_pattern <- gsub("\\.", "\\\\.", res_str)
      leiden_hits <- grep(
        paste0("^leiden.*", res_pattern, "$"),
        colnames(tax),
        value = TRUE
      )
    }
    
    if (length(leiden_hits) != 1) {
      stop(
        "Could not uniquely resolve the Leiden column for ",
        cell_type,
        " selected resolution ",
        res_str,
        ". Matches: ",
        paste(leiden_hits, collapse = ", ")
      )
    }
    
    leiden_col <- leiden_hits[1]
    selected_leiden_col_by_celltype[[cell_type]] <- leiden_col
    
    auto_anno_col <- paste0("leiden_", res_str, "_anno_coarse")
    auto_source_col <- paste0("leiden_", res_str, "_anno_source")
    auto_weak_col <- paste0("leiden_", res_str, "_weak_association")
    auto_correlation_col <- paste0("leiden_", res_str, "_selected_correlation")
    auto_margin_col <- paste0("leiden_", res_str, "_selected_margin")
    
    required_selected_cols <- c(
      leiden_col,
      auto_anno_col,
      auto_source_col,
      auto_weak_col,
      auto_correlation_col,
      auto_margin_col
    )
    
    missing_selected_cols <- setdiff(required_selected_cols, colnames(tax))
    
    if (length(missing_selected_cols) > 0) {
      stop(
        "Missing selected-resolution columns for ",
        cell_type,
        " res ",
        res_str,
        ": ",
        paste(missing_selected_cols, collapse = ", ")
      )
    }
    
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
      
      cluster_review <- tax %>%
        transmute(
          cluster = as.character(.data[[leiden_col]]),
          auto_anno = as.character(.data[[auto_anno_col]]),
          auto_source = as.character(.data[[auto_source_col]]),
          weak_association = as.logical(.data[[auto_weak_col]]),
          selected_correlation = as.numeric(.data[[auto_correlation_col]]),
          selected_margin = as.numeric(.data[[auto_margin_col]]),
          donor = .data[[donor_col]]
        ) %>%
        mutate(
          cluster = str_remove(cluster, "^query\\."),
          cluster = str_remove(cluster, "^RNA\\.g"),
          cluster = str_remove(cluster, "^RNA\\."),
          cluster = str_remove(cluster, "^g(?=[0-9]+$)")
        ) %>%
        group_by(cluster) %>%
        summarize(
          auto_anno = first(auto_anno),
          auto_source = first(auto_source),
          weak_association = first(weak_association),
          selected_correlation = first(selected_correlation),
          selected_margin = first(selected_margin),
          n_cells = n(),
          n_donors = n_distinct(donor, na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      cluster_review <- tax %>%
        transmute(
          cluster = as.character(.data[[leiden_col]]),
          auto_anno = as.character(.data[[auto_anno_col]]),
          auto_source = as.character(.data[[auto_source_col]]),
          weak_association = as.logical(.data[[auto_weak_col]]),
          selected_correlation = as.numeric(.data[[auto_correlation_col]]),
          selected_margin = as.numeric(.data[[auto_margin_col]])
        ) %>%
        mutate(
          cluster = str_remove(cluster, "^query\\."),
          cluster = str_remove(cluster, "^RNA\\.g"),
          cluster = str_remove(cluster, "^RNA\\."),
          cluster = str_remove(cluster, "^g(?=[0-9]+$)")
        ) %>%
        group_by(cluster) %>%
        summarize(
          auto_anno = first(auto_anno),
          auto_source = first(auto_source),
          weak_association = first(weak_association),
          selected_correlation = first(selected_correlation),
          selected_margin = first(selected_margin),
          n_cells = n(),
          n_donors = NA_integer_,
          .groups = "drop"
        )
    }
    
    if (all(str_detect(cluster_review$cluster, "^[0-9]+$"))) {
      cluster_review <- cluster_review %>%
        mutate(.cluster_number = as.integer(cluster)) %>%
        arrange(.cluster_number) %>%
        dplyr::select(-.cluster_number)
    } else {
      cluster_review <- cluster_review %>%
        arrange(cluster)
    }
    
    cluster_review <- cluster_review %>%
      mutate(
        cell_type = cell_type,
        selected_resolution = res_str,
        .before = 1
      )
    
    cluster_review_by_celltype[[cell_type]] <- cluster_review
    
    review_file <- file.path(
      review_root,
      paste0(
        cell_type,
        "_selected_res",
        res_str,
        "_cluster_review.csv"
      )
    )
    
    write.csv(
      cluster_review,
      review_file,
      row.names = FALSE
    )
    
    mapping_template <- cluster_review %>%
      transmute(
        cell_type,
        selected_resolution,
        cluster,
        auto_anno,
        auto_source,
        weak_association,
        selected_correlation,
        selected_margin,
        n_cells,
        n_donors,
        anno = ""
      )
    
    mapping_template_file <- file.path(
      review_root,
      paste0(
        cell_type,
        "_selected_res",
        res_str,
        "_anno_mapping_template.csv"
      )
    )
    
    write.csv(
      mapping_template,
      mapping_template_file,
      row.names = FALSE
    )
    
    # Also write an R-code template that can be pasted directly into section 2.
    template_map_lines <- paste0(
      "  \"",
      cluster_review$cluster,
      "\" = \"\""
    )
    
    if (length(template_map_lines) > 1) {
      template_map_lines[seq_len(length(template_map_lines) - 1)] <- paste0(
        template_map_lines[seq_len(length(template_map_lines) - 1)],
        ","
      )
    }
    
    anno_map_template_lines <- c(
      paste0(cell_type, "_anno_map <- c("),
      template_map_lines,
      ")"
    )
    
    anno_map_template_file <- file.path(
      review_root,
      paste0(
        cell_type,
        "_selected_res",
        res_str,
        "_anno_map_template.R"
      )
    )
    
    writeLines(
      anno_map_template_lines,
      anno_map_template_file
    )
    
    message("\n------------------------------------------------------------")
    message(cell_type, " selected resolution: ", res_str)
    message("Cluster review: ", review_file)
    message("Mapping CSV template: ", mapping_template_file)
    message("R mapping template: ", anno_map_template_file)
    message("------------------------------------------------------------")
    print(as.data.frame(cluster_review))
    
    # -------------------------------------------------------------------------
    # Validate the manually entered cluster -> anno mapping.
    # -------------------------------------------------------------------------
    
    anno_map <- anno_map_by_celltype[[cell_type]]
    mapping_problems <- character(0)
    
    if (length(anno_map) == 0) {
      mapping_problems <- c(mapping_problems, "mapping is empty")
    } else {
      if (is.null(names(anno_map))) {
        mapping_problems <- c(mapping_problems, "mapping has no cluster-number names")
      } else {
        if (any(is.na(names(anno_map)) | trimws(names(anno_map)) == "")) {
          mapping_problems <- c(mapping_problems, "mapping has blank cluster-number names")
        }
        
        if (anyDuplicated(names(anno_map)) > 0) {
          duplicate_mapping_clusters <- unique(
            names(anno_map)[duplicated(names(anno_map))]
          )
          mapping_problems <- c(
            mapping_problems,
            paste0(
              "duplicate mapped clusters: ",
              paste(duplicate_mapping_clusters, collapse = ", ")
            )
          )
        }
        
        missing_mapping_clusters <- setdiff(
          cluster_review$cluster,
          names(anno_map)
        )
        
        extra_mapping_clusters <- setdiff(
          names(anno_map),
          cluster_review$cluster
        )
        
        if (length(missing_mapping_clusters) > 0) {
          mapping_problems <- c(
            mapping_problems,
            paste0(
              "missing clusters: ",
              paste(missing_mapping_clusters, collapse = ", ")
            )
          )
        }
        
        if (length(extra_mapping_clusters) > 0) {
          mapping_problems <- c(
            mapping_problems,
            paste0(
              "clusters not present at selected resolution: ",
              paste(extra_mapping_clusters, collapse = ", ")
            )
          )
        }
      }
      
      if (any(is.na(anno_map) | trimws(anno_map) == "")) {
        mapping_problems <- c(mapping_problems, "one or more final anno names are blank")
      }
    }
    
    if (length(mapping_problems) == 0) {
      mapping_ready[[cell_type]] <- TRUE
    } else {
      mapping_ready[[cell_type]] <- FALSE
      message("\n", cell_type, " final anno mapping is not complete yet:")
      for (problem_i in seq_along(mapping_problems)) {
        message("  - ", mapping_problems[problem_i])
      }
    }
  }
}

# ==============================================================================
# 8. Assign the manually curated final "anno" column
# ==============================================================================

all_anno_mappings_complete <- resolution_selection_complete && all(mapping_ready)

if (resolution_selection_complete && !all_anno_mappings_complete) {
  message("\n============================================================")
  message("FINAL anno MAPPINGS NOT YET COMPLETE")
  message("============================================================")
  message("The cluster-review and mapping-template files were written.")
  message("Fill the three *_anno_map objects in section 2 and rerun.")
  message("No cluster name will be guessed automatically.")
}

if (all_anno_mappings_complete) {
  for (cell_type in cell_types) {
    tax <- tax_by_celltype[[cell_type]]
    res_str <- selected_resolution[[cell_type]]
    leiden_col <- selected_leiden_col_by_celltype[[cell_type]]
    anno_map <- anno_map_by_celltype[[cell_type]]
    
    auto_anno_col <- paste0("leiden_", res_str, "_anno_coarse")
    auto_source_col <- paste0("leiden_", res_str, "_anno_source")
    auto_weak_col <- paste0("leiden_", res_str, "_weak_association")
    auto_correlation_col <- paste0("leiden_", res_str, "_selected_correlation")
    auto_margin_col <- paste0("leiden_", res_str, "_selected_margin")
    
    tax <- tax %>%
      mutate(
        selected_resolution = res_str,
        selected_cluster = as.character(.data[[leiden_col]]),
        selected_cluster = str_remove(selected_cluster, "^query\\."),
        selected_cluster = str_remove(selected_cluster, "^RNA\\.g"),
        selected_cluster = str_remove(selected_cluster, "^RNA\\."),
        selected_cluster = str_remove(selected_cluster, "^g(?=[0-9]+$)"),
        anno_auto = as.character(.data[[auto_anno_col]]),
        anno_auto_source = as.character(.data[[auto_source_col]]),
        anno_auto_weak_association = as.logical(.data[[auto_weak_col]]),
        anno_auto_correlation = as.numeric(.data[[auto_correlation_col]]),
        anno_auto_margin = as.numeric(.data[[auto_margin_col]]),
        anno = unname(anno_map[selected_cluster])
      )
    
    if (any(is.na(tax$anno) | trimws(tax$anno) == "")) {
      missing_anno_clusters <- tax %>%
        filter(is.na(anno) | trimws(anno) == "") %>%
        distinct(selected_cluster) %>%
        pull(selected_cluster)
      
      stop(
        "Missing final anno after mapping for ",
        cell_type,
        " clusters: ",
        paste(missing_anno_clusters, collapse = ", ")
      )
    }
    
    final_tax_by_celltype[[cell_type]] <- tax
    
    write.csv(
      tax,
      file.path(
        final_table_root,
        paste0(
          "seaad_ec_",
          cell_type,
          "_selected_res",
          res_str,
          "_taxonomy_with_anno.csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  ExN_tax <- final_tax_by_celltype[["ExN"]]
  InN_tax <- final_tax_by_celltype[["InN"]]
  NonNeuronal_tax <- final_tax_by_celltype[["NonNeuronal"]]
  
  message("\n============================================================")
  message("FINAL anno COLUMN ASSIGNED")
  message("============================================================")
  message("The selected-resolution taxonomies with final anno were written to:")
  message(final_table_root)
  message("Now write the anno_coarse case_when() blocks below.")
}

# ==============================================================================
# 9. DEFINE anno_coarse FROM FINAL anno AND SAVE FINAL TAXONOMIES
# ==============================================================================
#
# Requires:
#   - selected resolutions for all three compartments
#   - complete ExN_anno_map, InN_anno_map, NonNeuronal_anno_map
#   - finalize_anno_coarse <- TRUE
#
# Rules:
#   ExN:
#       "ExN WWW XXX" -> "ExN WWW"
#
#   InN:
#       "InN WWW XXX" -> "InN WWW"
#
#   NonNeuronal:
#       remove a trailing number
#       "Astro3" -> "Astro"
#       "Micro2" -> "Micro"
#       "Oligo3" -> "Oligo"
#       "BEC venous" -> "BEC venous"
#       "Peri" -> "Peri"
# ==============================================================================

if (all_anno_mappings_complete && finalize_anno_coarse) {
  
  # ---------------------------------------------------------------------------
  # ExN
  # ---------------------------------------------------------------------------
  
  ExN_tax <- ExN_tax %>%
    mutate(
      anno_coarse = stringr::word(anno, 1, 2)
    )
  
  # ---------------------------------------------------------------------------
  # InN
  # ---------------------------------------------------------------------------
  
  InN_tax <- InN_tax %>%
    mutate(
      anno_coarse = stringr::word(anno, 1, 2)
    )
  
  # ---------------------------------------------------------------------------
  # NonNeuronal
  # ---------------------------------------------------------------------------
  
  NonNeuronal_tax <- NonNeuronal_tax %>%
    mutate(
      anno_coarse = stringr::str_remove(anno, "\\d+$"),
      anno_coarse = stringr::str_trim(anno_coarse)
    )
  
  # ---------------------------------------------------------------------------
  # Validate anno_coarse
  # ---------------------------------------------------------------------------
  
  missing_coarse_ExN <- ExN_tax %>%
    filter(is.na(anno_coarse) | trimws(anno_coarse) == "") %>%
    distinct(anno) %>%
    pull(anno)
  
  missing_coarse_InN <- InN_tax %>%
    filter(is.na(anno_coarse) | trimws(anno_coarse) == "") %>%
    distinct(anno) %>%
    pull(anno)
  
  missing_coarse_NonNeuronal <- NonNeuronal_tax %>%
    filter(is.na(anno_coarse) | trimws(anno_coarse) == "") %>%
    distinct(anno) %>%
    pull(anno)
  
  if (length(missing_coarse_ExN) > 0) {
    stop(
      "ExN anno values missing anno_coarse: ",
      paste(missing_coarse_ExN, collapse = ", ")
    )
  }
  
  if (length(missing_coarse_InN) > 0) {
    stop(
      "InN anno values missing anno_coarse: ",
      paste(missing_coarse_InN, collapse = ", ")
    )
  }
  
  if (length(missing_coarse_NonNeuronal) > 0) {
    stop(
      "NonNeuronal anno values missing anno_coarse: ",
      paste(missing_coarse_NonNeuronal, collapse = ", ")
    )
  }
  
  # ---------------------------------------------------------------------------
  # Show final anno -> anno_coarse mappings
  # ---------------------------------------------------------------------------
  
  message("\n============================================================")
  message("FINAL ExN anno -> anno_coarse")
  message("============================================================")
  
  print(
    ExN_tax %>%
      distinct(anno, anno_coarse) %>%
      arrange(anno_coarse, anno)
  )
  
  message("\n============================================================")
  message("FINAL InN anno -> anno_coarse")
  message("============================================================")
  
  print(
    InN_tax %>%
      distinct(anno, anno_coarse) %>%
      arrange(anno_coarse, anno)
  )
  
  message("\n============================================================")
  message("FINAL NonNeuronal anno -> anno_coarse")
  message("============================================================")
  
  print(
    NonNeuronal_tax %>%
      distinct(anno, anno_coarse) %>%
      arrange(anno_coarse, anno)
  )
  
  # ---------------------------------------------------------------------------
  # Update final taxonomy list
  # ---------------------------------------------------------------------------
  
  final_tax_by_celltype[["ExN"]] <- ExN_tax
  final_tax_by_celltype[["InN"]] <- InN_tax
  final_tax_by_celltype[["NonNeuronal"]] <- NonNeuronal_tax
  
  # ---------------------------------------------------------------------------
  # Define exact final output files
  # ---------------------------------------------------------------------------
  
  ExN_final_file <- file.path(
    final_table_root,
    paste0(
      "seaad_ec_ExN_selected_res",
      selected_resolution[["ExN"]],
      "_final_taxonomy.csv"
    )
  )
  
  InN_final_file <- file.path(
    final_table_root,
    paste0(
      "seaad_ec_InN_selected_res",
      selected_resolution[["InN"]],
      "_final_taxonomy.csv"
    )
  )
  
  NonNeuronal_final_file <- file.path(
    final_table_root,
    paste0(
      "seaad_ec_NonNeuronal_selected_res",
      selected_resolution[["NonNeuronal"]],
      "_final_taxonomy.csv"
    )
  )
  
  final_lookup_file <- file.path(
    final_table_root,
    "seaad_ec_final_annotation_lookup.csv"
  )
  
  combined_final_file <- file.path(
    final_table_root,
    "seaad_ec_all_celltypes_final_taxonomy.csv"
  )
  
  # ---------------------------------------------------------------------------
  # Save final per-compartment taxonomy tables
  # ---------------------------------------------------------------------------
  
  write.csv(
    ExN_tax,
    ExN_final_file,
    row.names = FALSE
  )
  
  write.csv(
    InN_tax,
    InN_final_file,
    row.names = FALSE
  )
  
  write.csv(
    NonNeuronal_tax,
    NonNeuronal_final_file,
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Save compact cluster -> anno -> anno_coarse lookup
  # ---------------------------------------------------------------------------
  
  final_annotation_lookup <- bind_rows(
    ExN_tax %>%
      distinct(
        selected_resolution,
        selected_cluster,
        anno,
        anno_coarse
      ) %>%
      mutate(
        cell_type = "ExN",
        .before = 1
      ),
    
    InN_tax %>%
      distinct(
        selected_resolution,
        selected_cluster,
        anno,
        anno_coarse
      ) %>%
      mutate(
        cell_type = "InN",
        .before = 1
      ),
    
    NonNeuronal_tax %>%
      distinct(
        selected_resolution,
        selected_cluster,
        anno,
        anno_coarse
      ) %>%
      mutate(
        cell_type = "NonNeuronal",
        .before = 1
      )
  )
  
  if (
    all(
      stringr::str_detect(
        final_annotation_lookup$selected_cluster,
        "^[0-9]+$"
      )
    )
  ) {
    
    final_annotation_lookup <- final_annotation_lookup %>%
      mutate(
        .cluster_number = as.integer(selected_cluster)
      ) %>%
      arrange(
        factor(cell_type, levels = cell_types),
        .cluster_number
      ) %>%
      dplyr::select(-.cluster_number)
    
  } else {
    
    final_annotation_lookup <- final_annotation_lookup %>%
      arrange(
        factor(cell_type, levels = cell_types),
        selected_cluster
      )
  }
  
  write.csv(
    final_annotation_lookup,
    final_lookup_file,
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Save combined cell-level taxonomy
  # ---------------------------------------------------------------------------
  
  final_taxonomy_all <- bind_rows(
    ExN_tax %>%
      mutate(
        broad_cell_type = "ExN",
        .before = 1
      ),
    
    InN_tax %>%
      mutate(
        broad_cell_type = "InN",
        .before = 1
      ),
    
    NonNeuronal_tax %>%
      mutate(
        broad_cell_type = "NonNeuronal",
        .before = 1
      )
  )
  
  write.csv(
    final_taxonomy_all,
    combined_final_file,
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Verify that the files were actually written
  # ---------------------------------------------------------------------------
  
  final_files <- c(
    ExN_final_file,
    InN_final_file,
    NonNeuronal_final_file,
    final_lookup_file,
    combined_final_file
  )
  
  missing_final_files <- final_files[!file.exists(final_files)]
  
  if (length(missing_final_files) > 0) {
    stop(
      "One or more final taxonomy files were not written:\n",
      paste(missing_final_files, collapse = "\n")
    )
  }
  
  # ---------------------------------------------------------------------------
  # Completion message
  # ---------------------------------------------------------------------------
  
  message("\n============================================================")
  message("FINAL TAXONOMY COMPLETE")
  message("============================================================")
  
  message("\nExN:")
  message(ExN_final_file)
  
  message("\nInN:")
  message(InN_final_file)
  
  message("\nNonNeuronal:")
  message(NonNeuronal_final_file)
  
  message("\nAnnotation lookup:")
  message(final_lookup_file)
  
  message("\nCombined taxonomy:")
  message(combined_final_file)
}


if (all_anno_mappings_complete && !finalize_anno_coarse) {
  
  message("\n============================================================")
  message("anno IS READY; anno_coarse HAS NOT BEEN FINALIZED")
  message("============================================================")
  message("Set finalize_anno_coarse <- TRUE near the top of the script.")
  message("Then rerun the script.")
}

