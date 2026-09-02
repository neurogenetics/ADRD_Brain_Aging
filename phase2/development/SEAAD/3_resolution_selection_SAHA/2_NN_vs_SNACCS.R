library(SAHA)
library(tidyverse)
library(clustree)

# ==============================================================================
# SEA-AD NonNeuronal resolution validation against SNACCS
# SAHA marker-free / AvgExp only
# ==============================================================================

# Existing SEA-AD resolution-sweep outputs
input_dir <- "/data/ADRD/brain_aging/phase2/public/seaad/resolution_selection"

# SNACCS NonNeuronal average-expression reference
snaccs_avgexp_file <- "/data/ADRD/atera_consulting/data/snaccs_saha/NN_AvgExp.csv"

# Writable outputs under exploration
output_root <- "/data/ADRD/brain_aging/exploration/res/Phase2/SEAAD_EC/resolution_selection/SNACCS_NN_markerfree"
plot_dir <- "/data/ADRD/brain_aging/exploration/plots/Phase2/SEAAD_EC/resolution_selection/SNACCS_NN_markerfree"
saha_dir <- file.path(output_root, "SAHA")
table_dir <- file.path(output_root, "tables")

resolution_values <- seq(0.1, 0.9, by = 0.1)
resolution_strings <- sprintf("%.1f", resolution_values)

overwrite_saha <- FALSE

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(saha_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load and format SNACCS AvgExp reference
# ==============================================================================

if (!file.exists(snaccs_avgexp_file)) {
  stop("Missing SNACCS AvgExp reference: ", snaccs_avgexp_file)
}

SNACCS_NN_AvgExp <- read.csv(snaccs_avgexp_file, check.names = FALSE)

if (nrow(SNACCS_NN_AvgExp) == 0 || ncol(SNACCS_NN_AvgExp) < 2) {
  stop("Unexpected SNACCS AvgExp dimensions: ", nrow(SNACCS_NN_AvgExp), " rows x ", ncol(SNACCS_NN_AvgExp), " columns")
}

# SAHA AvgExp expects genes as row names and reference annotations as columns.
# Handle the common CSV format in which gene names were written as one non-numeric
# column, usually the first column. Use a column index so an empty CSV header is safe.
default_rownames <- identical(rownames(SNACCS_NN_AvgExp), as.character(seq_len(nrow(SNACCS_NN_AvgExp))))
snaccs_numeric_columns <- vapply(SNACCS_NN_AvgExp, is.numeric, logical(1))

candidate_gene_columns <- c(
  "gene",
  "Gene",
  "genes",
  "Genes",
  "gene_name",
  "gene_names",
  "feature",
  "features",
  "X"
)

candidate_gene_index <- which(colnames(SNACCS_NN_AvgExp) %in% candidate_gene_columns)
non_numeric_index <- which(!snaccs_numeric_columns)

if (length(candidate_gene_index) > 0) {
  snaccs_gene_col_index <- candidate_gene_index[1]
} else if (length(non_numeric_index) == 1) {
  snaccs_gene_col_index <- non_numeric_index[1]
} else {
  snaccs_gene_col_index <- NA_integer_
}

if (!is.na(snaccs_gene_col_index)) {
  snaccs_gene_ids <- as.character(SNACCS_NN_AvgExp[[snaccs_gene_col_index]])
  
  if (any(is.na(snaccs_gene_ids) | snaccs_gene_ids == "")) {
    stop("Missing gene identifiers in the inferred SNACCS gene column.")
  }
  
  if (anyDuplicated(snaccs_gene_ids) > 0) {
    duplicate_snaccs_genes <- unique(snaccs_gene_ids[duplicated(snaccs_gene_ids)])
    stop("Duplicate SNACCS gene identifiers: ", paste(head(duplicate_snaccs_genes, 20), collapse = ", "))
  }
  
  SNACCS_NN_AvgExp <- SNACCS_NN_AvgExp[, -snaccs_gene_col_index, drop = FALSE]
  rownames(SNACCS_NN_AvgExp) <- snaccs_gene_ids
} else if (default_rownames) {
  stop(
    "Could not identify gene names in SNACCS_NN_AvgExp. ",
    "SAHA AvgExp requires genes as row names and SNACCS annotations as columns."
  )
}

snaccs_non_numeric_cols <- colnames(SNACCS_NN_AvgExp)[!vapply(SNACCS_NN_AvgExp, is.numeric, logical(1))]
if (length(snaccs_non_numeric_cols) > 0) {
  stop("Non-numeric SNACCS AvgExp columns after gene-column removal: ", paste(snaccs_non_numeric_cols, collapse = ", "))
}

if (anyDuplicated(colnames(SNACCS_NN_AvgExp)) > 0) {
  duplicate_snaccs_annotations <- unique(colnames(SNACCS_NN_AvgExp)[duplicated(colnames(SNACCS_NN_AvgExp))])
  stop("Duplicate SNACCS annotation columns: ", paste(duplicate_snaccs_annotations, collapse = ", "))
}

message("SNACCS reference: ", nrow(SNACCS_NN_AvgExp), " genes x ", ncol(SNACCS_NN_AvgExp), " annotations")

# ==============================================================================
# 2. Load SEA-AD NonNeuronal variable features
# ==============================================================================

varfeat_file <- file.path(input_dir, "seaad_ec_NonNeuronal_variable_features.txt")

if (!file.exists(varfeat_file)) {
  stop("Missing SEA-AD NonNeuronal variable-feature file: ", varfeat_file)
}

seaad_varfeat <- unique(trimws(readLines(varfeat_file)))
seaad_varfeat <- seaad_varfeat[nzchar(seaad_varfeat)]

if (length(seaad_varfeat) == 0) {
  stop("No variable features found in: ", varfeat_file)
}

message("SEA-AD NonNeuronal variable features: ", length(seaad_varfeat))

# ==============================================================================
# 3. Run marker-free SAHA for resolutions 0.1-0.9
# ==============================================================================

for (res_i in seq_along(resolution_values)) {
  res <- resolution_values[res_i]
  res_str <- resolution_strings[res_i]
  
  message("\n========== SNACCS marker-free: NonNeuronal res ", res_str, " ==========")
  
  avgexp_file <- file.path(
    input_dir,
    paste0("seaad_ec_NonNeuronal_avgexp_res", res_str, ".csv")
  )
  
  ann_file <- file.path(
    saha_dir,
    paste0("seaad_ec_NonNeuronal_res", res_str, "_SAHAvsSNACCS_markerfree_ANN.RData")
  )
  
  auto_file <- file.path(
    saha_dir,
    paste0("seaad_ec_NonNeuronal_res", res_str, "_SAHAvsSNACCS_markerfree_auto.csv")
  )
  
  html_file <- file.path(
    saha_dir,
    paste0("seaad_ec_NonNeuronal_res", res_str, "_SAHAvsSNACCS_markerfree.html")
  )
  
  overlap_file <- file.path(
    saha_dir,
    paste0("seaad_ec_NonNeuronal_res", res_str, "_SAHAvsSNACCS_markerfree_genes.txt")
  )
  
  if (!file.exists(avgexp_file)) {
    stop("Missing SEA-AD AvgExp input: ", avgexp_file)
  }
  
  if (
    !overwrite_saha &&
    file.exists(ann_file) &&
    file.exists(auto_file) &&
    file.exists(html_file) &&
    file.exists(overlap_file)
  ) {
    message("Skipping completed resolution ", res_str)
    next
  }
  
  # ---------------------------------------------------------------------------
  # Format SEA-AD AvgExp as genes x query clusters
  # ---------------------------------------------------------------------------
  
  query_avgexp_raw <- read.csv(avgexp_file, check.names = FALSE)
  
  if (nrow(query_avgexp_raw) == 0 || ncol(query_avgexp_raw) < 2) {
    stop("Unexpected SEA-AD AvgExp dimensions in: ", avgexp_file)
  }
  
  query_cluster_ids <- as.character(query_avgexp_raw[[1]])
  
  if (any(is.na(query_cluster_ids) | query_cluster_ids == "")) {
    stop("Missing query cluster identifiers in: ", avgexp_file)
  }
  
  if (anyDuplicated(query_cluster_ids) > 0) {
    duplicate_query_clusters <- unique(query_cluster_ids[duplicated(query_cluster_ids)])
    stop("Duplicate query cluster identifiers at res ", res_str, ": ", paste(duplicate_query_clusters, collapse = ", "))
  }
  
  query_cluster_ids <- ifelse(
    grepl("^RNA\\.g", query_cluster_ids),
    query_cluster_ids,
    paste0("RNA.g", query_cluster_ids)
  )
  
  query_avgexp_raw <- query_avgexp_raw[, -1, drop = FALSE]
  
  query_non_numeric_cols <- colnames(query_avgexp_raw)[!vapply(query_avgexp_raw, is.numeric, logical(1))]
  if (length(query_non_numeric_cols) > 0) {
    stop(
      "Non-numeric SEA-AD AvgExp gene columns at res ",
      res_str,
      ": ",
      paste(head(query_non_numeric_cols, 20), collapse = ", ")
    )
  }
  
  rownames(query_avgexp_raw) <- query_cluster_ids
  query_avgexp <- as.data.frame(t(query_avgexp_raw), check.names = FALSE)
  
  if (anyDuplicated(rownames(query_avgexp)) > 0) {
    duplicate_query_genes <- unique(rownames(query_avgexp)[duplicated(rownames(query_avgexp))])
    stop("Duplicate SEA-AD AvgExp genes at res ", res_str, ": ", paste(head(duplicate_query_genes, 20), collapse = ", "))
  }
  
  # Use SEA-AD variable features that are available in both expression matrices.
  common_varfeats <- Reduce(
    intersect,
    list(
      seaad_varfeat,
      rownames(query_avgexp),
      rownames(SNACCS_NN_AvgExp)
    )
  )
  
  if (length(common_varfeats) == 0) {
    stop("No overlapping SEA-AD variable features between query and SNACCS at res ", res_str)
  }
  
  message(
    "Genes used for marker-free comparison: ",
    length(common_varfeats),
    " / ",
    length(seaad_varfeat),
    " SEA-AD variable features"
  )
  
  writeLines(common_varfeats, overlap_file)
  
  # ---------------------------------------------------------------------------
  # SAHA marker-free only
  # ---------------------------------------------------------------------------
  
  ann <- Create_SAHA_object(
    query = query_avgexp,
    db = SNACCS_NN_AvgExp,
    data_type = "AvgExp"
  )
  
  ann <- Initialize_MarkerFree(ann = ann)
  ann <- Downsample(ann, custom_ds = common_varfeats)
  ann <- NormalizeDS(ann, assay_query = "RNA")
  ann <- CorrelateDS(ann)
  ann <- Create_MarkerFree_Viz(ann)
  
  # AvgExp-only auto annotation uses the maximum marker-free correlation.
  auto_raw <- AutoAnnotate(ann, data_type = "AvgExp")
  
  if (!all(c("cluster", "best_match", "correlation") %in% colnames(auto_raw))) {
    stop("Unexpected AvgExp AutoAnnotate output at res ", res_str)
  }
  
  # Preserve the raw SAHA best-match string and also make a cleaned annotation.
  auto <- auto_raw %>%
    mutate(
      cluster = sub("^RNA\\.g", "", as.character(cluster)),
      best_match_raw = as.character(best_match),
      best_match = str_remove(best_match_raw, "^db\\."),
      best_match = str_replace_all(best_match, "\\.", " "),
      resolution = res
    )
  
  Generate_SAHA_Report(
    ann,
    auto = auto_raw,
    output_file = html_file
  )
  
  save(ann, file = ann_file)
  write.csv(auto, auto_file, row.names = FALSE)
  
  message("Completed NonNeuronal res ", res_str)
  
  rm(
    query_avgexp_raw,
    query_cluster_ids,
    query_non_numeric_cols,
    query_avgexp,
    common_varfeats,
    ann,
    auto_raw,
    auto
  )
}

# ==============================================================================
# 4. Combine marker-free SNACCS calls across resolutions
# ==============================================================================

auto_list <- vector("list", length(resolution_values))

for (res_i in seq_along(resolution_values)) {
  res <- resolution_values[res_i]
  res_str <- resolution_strings[res_i]
  
  auto_file <- file.path(
    saha_dir,
    paste0("seaad_ec_NonNeuronal_res", res_str, "_SAHAvsSNACCS_markerfree_auto.csv")
  )
  
  if (!file.exists(auto_file)) {
    stop("Missing completed SNACCS auto-annotation file: ", auto_file)
  }
  
  auto_df <- read.csv(auto_file, check.names = FALSE)
  
  if (!all(c("cluster", "best_match", "correlation") %in% colnames(auto_df))) {
    stop("Unexpected auto-annotation columns in: ", auto_file)
  }
  
  auto_list[[res_i]] <- auto_df %>%
    mutate(
      cluster = as.character(cluster),
      resolution = res
    )
}

combined_auto <- bind_rows(auto_list)

write.csv(
  combined_auto,
  file.path(table_dir, "seaad_ec_NonNeuronal_SNACCS_markerfree_calls_all_resolutions.csv"),
  row.names = FALSE
)

# ==============================================================================
# 5. Attach SNACCS marker-free calls to cell-level resolution assignments
# ==============================================================================

metadata_file <- file.path(
  input_dir,
  "seaad_ec_NonNeuronal_resolution_assignments_obs.csv"
)

if (!file.exists(metadata_file)) {
  stop("Missing SEA-AD NonNeuronal resolution assignments: ", metadata_file)
}

cell_metadata <- read.csv(metadata_file, check.names = FALSE)

# Resolve Leiden assignment columns before adding the new SNACCS columns.
leiden_cols <- setNames(rep(NA_character_, length(resolution_strings)), resolution_strings)

for (res_i in seq_along(resolution_strings)) {
  res_str <- resolution_strings[res_i]
  
  candidates <- c(
    paste0("leiden_", res_str),
    paste0("leiden_scvi_", res_str),
    paste0("leiden_res", res_str),
    paste0("leiden_res_", res_str)
  )
  
  hits <- candidates[candidates %in% colnames(cell_metadata)]
  
  if (length(hits) == 0) {
    res_pattern <- gsub("\\.", "\\\\.", res_str)
    hits <- grep(
      paste0("^leiden.*", res_pattern, "$"),
      colnames(cell_metadata),
      value = TRUE
    )
  }
  
  if (length(hits) != 1) {
    stop(
      "Could not uniquely resolve the Leiden column for res ",
      res_str,
      ". Matches: ",
      paste(hits, collapse = ", ")
    )
  }
  
  leiden_cols[[res_str]] <- hits[1]
}

# Resolve donor column if present so the summary can report donor representation.
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

donor_hits <- donor_candidates[donor_candidates %in% colnames(cell_metadata)]

if (length(donor_hits) > 0) {
  donor_col <- donor_hits[1]
} else {
  donor_col <- NA_character_
}

cluster_summary_list <- vector("list", length(resolution_values))

for (res_i in seq_along(resolution_values)) {
  res <- resolution_values[res_i]
  res_str <- resolution_strings[res_i]
  leiden_col <- leiden_cols[[res_str]]
  
  lookup_df <- combined_auto %>%
    filter(resolution == res) %>%
    mutate(cluster = as.character(cluster))
  
  if (anyDuplicated(lookup_df$cluster) > 0) {
    stop("Duplicate SNACCS marker-free calls at res ", res_str)
  }
  
  metadata_clusters <- unique(sub("^RNA\\.g", "", as.character(cell_metadata[[leiden_col]])))
  metadata_clusters <- metadata_clusters[!is.na(metadata_clusters)]
  unmatched_clusters <- setdiff(metadata_clusters, lookup_df$cluster)
  
  if (length(unmatched_clusters) > 0) {
    stop(
      "Unmatched SEA-AD clusters at res ",
      res_str,
      ": ",
      paste(unmatched_clusters, collapse = ", ")
    )
  }
  
  annotation_col <- paste0("leiden_", res_str, "_SNACCS")
  correlation_col <- paste0("leiden_", res_str, "_SNACCS_correlation")
  
  join_lookup <- lookup_df %>%
    transmute(
      .cluster_join = cluster,
      .SNACCS_annotation = best_match,
      .SNACCS_correlation = correlation
    )
  
  cell_metadata <- cell_metadata %>%
    mutate(.cluster_join = sub("^RNA\\.g", "", as.character(.data[[leiden_col]]))) %>%
    left_join(join_lookup, by = ".cluster_join") %>%
    rename(
      !!annotation_col := .SNACCS_annotation,
      !!correlation_col := .SNACCS_correlation
    ) %>%
    dplyr::select(-.cluster_join)
  
  if (any(is.na(cell_metadata[[annotation_col]]))) {
    stop("Missing SNACCS annotation after cell-level join at res ", res_str)
  }
  
  cluster_counts <- cell_metadata %>%
    mutate(cluster = sub("^RNA\\.g", "", as.character(.data[[leiden_col]]))) %>%
    count(cluster, name = "n_cells")
  
  if (!is.na(donor_col)) {
    donor_counts <- cell_metadata %>%
      mutate(cluster = sub("^RNA\\.g", "", as.character(.data[[leiden_col]]))) %>%
      group_by(cluster) %>%
      summarize(
        n_donors = n_distinct(.data[[donor_col]], na.rm = TRUE),
        .groups = "drop"
      )
    
    cluster_counts <- cluster_counts %>%
      left_join(donor_counts, by = "cluster")
  } else {
    cluster_counts <- cluster_counts %>%
      mutate(n_donors = NA_integer_)
  }
  
  cluster_summary_list[[res_i]] <- lookup_df %>%
    dplyr::select(cluster, best_match, correlation) %>%
    left_join(cluster_counts, by = "cluster") %>%
    mutate(resolution = res)
}

cluster_summary <- bind_rows(cluster_summary_list)

write.csv(
  cell_metadata,
  file.path(table_dir, "seaad_ec_NonNeuronal_resolution_assignments_SNACCS_markerfree.csv"),
  row.names = FALSE
)

write.csv(
  cluster_summary,
  file.path(table_dir, "seaad_ec_NonNeuronal_SNACCS_markerfree_cluster_summary.csv"),
  row.names = FALSE
)

# ==============================================================================
# 6. Across-resolution summaries and plots
# ==============================================================================

if (!is.na(donor_col)) {
  annotation_summary <- cluster_summary %>%
    group_by(resolution, best_match) %>%
    summarize(
      total_cells = sum(n_cells),
      max_donors = max(n_donors),
      n_clusters = n_distinct(cluster),
      mean_correlation = weighted.mean(correlation, w = n_cells),
      min_correlation = min(correlation),
      max_correlation = max(correlation),
      .groups = "drop"
    ) %>%
    mutate(plot_label = paste0("D=", max_donors, "\nK=", n_clusters))
} else {
  annotation_summary <- cluster_summary %>%
    group_by(resolution, best_match) %>%
    summarize(
      total_cells = sum(n_cells),
      n_clusters = n_distinct(cluster),
      mean_correlation = weighted.mean(correlation, w = n_cells),
      min_correlation = min(correlation),
      max_correlation = max(correlation),
      .groups = "drop"
    ) %>%
    mutate(
      max_donors = NA_integer_,
      plot_label = paste0("K=", n_clusters)
    )
}

write.csv(
  annotation_summary,
  file.path(table_dir, "seaad_ec_NonNeuronal_SNACCS_markerfree_annotation_summary.csv"),
  row.names = FALSE
)

p_cells <- ggplot(
  annotation_summary,
  aes(x = best_match, y = total_cells, fill = best_match)
) +
  geom_col() +
  geom_text(aes(label = plot_label), vjust = -0.2, size = 2.5) +
  facet_wrap(~resolution, nrow = 3, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    x = NULL,
    y = "SEA-AD NonNeuronal cells",
    title = "SEA-AD NonNeuronal marker-free matches to SNACCS"
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  )

ggsave(
  file.path(plot_dir, "seaad_ec_NonNeuronal_SNACCS_markerfree_cells_across_resolutions.pdf"),
  p_cells,
  width = 12,
  height = 9
)

p_correlation <- ggplot(
  cluster_summary,
  aes(x = best_match, y = correlation)
) +
  geom_point(aes(size = n_cells), alpha = 0.75) +
  facet_wrap(~resolution, nrow = 3, scales = "free_x") +
  labs(
    x = NULL,
    y = "SAHA marker-free correlation",
    size = "Cells",
    title = "SEA-AD NonNeuronal cluster correlations to best SNACCS match"
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top"
  )

ggsave(
  file.path(plot_dir, "seaad_ec_NonNeuronal_SNACCS_markerfree_best_match_correlations.pdf"),
  p_correlation,
  width = 12,
  height = 9
)

# ==============================================================================
# 7. Clustree using the full SNACCS best-match label at each resolution
# ==============================================================================

snaccs_annotation_cols <- paste0("leiden_", resolution_strings, "_SNACCS")

missing_snaccs_annotation_cols <- setdiff(snaccs_annotation_cols, colnames(cell_metadata))
if (length(missing_snaccs_annotation_cols) > 0) {
  stop("Missing SNACCS annotation columns for clustree: ", paste(missing_snaccs_annotation_cols, collapse = ", "))
}

snaccs_clustree_df <- cell_metadata %>%
  dplyr::select(all_of(snaccs_annotation_cols))

png(
  file.path(plot_dir, "seaad_ec_NonNeuronal_SNACCS_markerfree_clustree.png"),
  width = 18,
  height = 10,
  units = "in",
  res = 300
)
print(clustree(snaccs_clustree_df, prefix = "leiden_", suffix = "_SNACCS"))
dev.off()

# ==============================================================================
# 8. Session information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(output_root, "SEAAD_NonNeuronal_SNACCS_markerfree_sessionInfo.txt")
)

message("\nCompleted SEA-AD NonNeuronal marker-free comparison against SNACCS.")
message("Results: ", output_root)
message("Plots: ", plot_dir)
