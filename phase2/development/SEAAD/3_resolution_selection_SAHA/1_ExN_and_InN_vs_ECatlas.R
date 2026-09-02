library(SAHA)
library(tidyverse)
library(clustree)

# ==============================================================================
# SEA-AD EC re-analysis: SAHA annotation and resolution selection
# ==============================================================================

# Existing read-only SEA-AD resolution sweep
input_dir <- "/data/ADRD/brain_aging/phase2/public/seaad/resolution_selection"

# Custom EC reference used for the primary samples
reference_dir <- "/data/ADRD/DLB_VHD/building_EC_reference/markers"

# New writable outputs
output_root <- "/data/ADRD/brain_aging/exploration/res/Phase2/SEAAD_EC/resolution_selection"
plot_dir <- "/data/ADRD/brain_aging/exploration/plots/Phase2/SEAAD_EC/resolution_selection"
saha_dir <- file.path(output_root, "SAHA")
table_dir <- file.path(output_root, "tables")

cell_types <- c("ExN", "InN")
resolution_values <- seq(0.1, 0.9, by = 0.1)
resolution_strings <- sprintf("%.1f", resolution_values)

overwrite_saha <- FALSE

# The uploaded primary workflow only establishes the ExN/InN custom-reference
# filenames. Leave these NA to let the script search for a unique NonNeuronal
# match in reference_dir, or fill them manually if its filenames differ.
# nonneuronal_db_markers <- NA_character_
# nonneuronal_db_avgexp <- NA_character_
# nonneuronal_db_varfeats <- NA_character_
# 
# Fill these in only AFTER inspecting the resolution-selection diagnostics.
selected_resolution <- c(ExN = NA_real_, InN = NA_real_)
selected_level <- c(ExN = "subtype", InN = "subtype")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(saha_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

completed_cell_types <- character(0)

# ==============================================================================
# 1. Run SAHA for every available cell compartment and resolution
# ==============================================================================

for (cell_type in cell_types) {
  message("\n========== SAHA: ", cell_type, " ==========")
  
  db_markers_file <- file.path(reference_dir, paste0(cell_type, "_cleaned_NAMED_markers.rds"))
  db_avgexp_file <- file.path(reference_dir, paste0(cell_type, "_cleaned_NAMED_avgexp.rds"))
  db_varfeats_file <- file.path(reference_dir, paste0(cell_type, "_cleaned_varfeats_for_hclust.rds"))
  
  if (cell_type == "NonNeuronal") {
    if (!is.na(nonneuronal_db_markers)) db_markers_file <- nonneuronal_db_markers
    if (!is.na(nonneuronal_db_avgexp)) db_avgexp_file <- nonneuronal_db_avgexp
    if (!is.na(nonneuronal_db_varfeats)) db_varfeats_file <- nonneuronal_db_varfeats
    
    if (!file.exists(db_markers_file)) {
      x <- list.files(reference_dir, pattern = "Non.*NAMED_markers\\.rds$", full.names = TRUE, ignore.case = TRUE)
      if (length(x) == 1) db_markers_file <- x[1]
    }
    
    if (!file.exists(db_avgexp_file)) {
      x <- list.files(reference_dir, pattern = "Non.*NAMED_avgexp\\.rds$", full.names = TRUE, ignore.case = TRUE)
      if (length(x) == 1) db_avgexp_file <- x[1]
    }
    
    if (!file.exists(db_varfeats_file)) {
      x <- list.files(reference_dir, pattern = "Non.*varfeats.*\\.rds$", full.names = TRUE, ignore.case = TRUE)
      if (length(x) == 1) db_varfeats_file <- x[1]
    }
  }
  
  reference_files <- c(db_markers_file, db_avgexp_file, db_varfeats_file)
  
  if (any(!file.exists(reference_files))) {
    missing_reference_files <- reference_files[!file.exists(reference_files)]
    
    if (cell_type == "NonNeuronal") {
      message("Skipping NonNeuronal because its custom EC reference could not be uniquely resolved.")
      message("Missing/unresolved: ", paste(missing_reference_files, collapse = "; "))
      message("Set nonneuronal_db_markers/nonneuronal_db_avgexp/nonneuronal_db_varfeats at the top and rerun.")
      next
    }
    
    stop(
      "Missing custom EC reference file(s) for ",
      cell_type,
      ": ",
      paste(missing_reference_files, collapse = "; ")
    )
  }
  
  db_markers <- readRDS(db_markers_file)
  db_avgexp <- readRDS(db_avgexp_file)
  db_varfeats <- readRDS(db_varfeats_file)
  
  for (res_i in seq_along(resolution_values)) {
    res_str <- resolution_strings[res_i]
    
    marker_file <- file.path(
      input_dir,
      paste0("seaad_ec_", cell_type, "_markers_res", res_str, ".csv")
    )
    
    avgexp_file <- file.path(
      input_dir,
      paste0("seaad_ec_", cell_type, "_avgexp_res", res_str, ".csv")
    )
    
    varfeat_file <- file.path(
      input_dir,
      paste0("seaad_ec_", cell_type, "_variable_features.txt")
    )
    
    ann_file <- file.path(
      saha_dir,
      paste0("seaad_ec_", cell_type, "_res", res_str, "_SAHAvsECatlas_ANN.RData")
    )
    
    auto_file <- file.path(
      saha_dir,
      paste0("seaad_ec_", cell_type, "_res", res_str, "_SAHAvsECatlas_auto.csv")
    )
    
    html_file <- file.path(
      saha_dir,
      paste0("seaad_ec_", cell_type, "_res", res_str, "_SAHAvsECatlas.html")
    )
    
    input_files <- c(marker_file, avgexp_file, varfeat_file)
    
    if (any(!file.exists(input_files))) {
      stop(
        "Missing SEA-AD input(s) for ",
        cell_type,
        " res ",
        res_str,
        ": ",
        paste(input_files[!file.exists(input_files)], collapse = "; ")
      )
    }
    
    if (
      !overwrite_saha &&
      file.exists(ann_file) &&
      file.exists(auto_file) &&
      file.exists(html_file)
    ) {
      message("Skipping completed: ", cell_type, " res ", res_str)
      next
    }
    
    message("Running: ", cell_type, " res ", res_str)
    
    # Marker table -> SAHA marker format
    query_markers_raw <- read.csv(marker_file, check.names = FALSE)
    
    required_marker_columns <- c(
      "pvals",
      "logfoldchanges",
      "pct_nz_group",
      "pct_nz_reference",
      "pvals_adj",
      "group",
      "names"
    )
    
    missing_marker_columns <- setdiff(
      required_marker_columns,
      colnames(query_markers_raw)
    )
    
    if (length(missing_marker_columns) > 0) {
      stop(
        "Unexpected marker format in ",
        marker_file,
        ". Missing: ",
        paste(missing_marker_columns, collapse = ", ")
      )
    }
    
    query_markers <- query_markers_raw %>%
      transmute(
        p_val = pvals,
        avg_log2FC = logfoldchanges,
        pct.1 = pct_nz_group,
        pct.2 = pct_nz_reference,
        p_val_adj = pvals_adj,
        cluster = group,
        gene = names
      ) %>%
      dplyr::select(
        p_val,
        avg_log2FC,
        pct.1,
        pct.2,
        p_val_adj,
        cluster,
        gene
      )
    
    finite_fc <- query_markers$avg_log2FC[
      is.finite(query_markers$avg_log2FC)
    ]
    
    if (length(finite_fc) == 0) {
      stop("No finite log-fold changes in ", marker_file)
    }
    
    max_fc <- max(finite_fc, na.rm = TRUE)
    min_fc <- min(finite_fc, na.rm = TRUE)
    
    query_markers <- query_markers %>%
      mutate(
        avg_log2FC = case_when(
          is.infinite(avg_log2FC) & avg_log2FC > 0 ~ max_fc,
          is.infinite(avg_log2FC) & avg_log2FC < 0 ~ min_fc,
          TRUE ~ avg_log2FC
        )
      )
    
    # AvgExp table -> genes x clusters matrix. The primary workflow used the
    # first column as the cluster ID and prefixed it with RNA.g.
    query_avgexp_raw <- read.csv(
      avgexp_file,
      check.names = FALSE
    )
    
    if (ncol(query_avgexp_raw) < 2) {
      stop("Unexpected AvgExp format in ", avgexp_file)
    }
    
    avgexp_cluster_ids <- as.character(query_avgexp_raw[[1]])
    
    avgexp_cluster_ids <- ifelse(
      grepl("^RNA\\.g", avgexp_cluster_ids),
      avgexp_cluster_ids,
      paste0("RNA.g", avgexp_cluster_ids)
    )
    
    query_avgexp_raw <- query_avgexp_raw[, -1, drop = FALSE]
    rownames(query_avgexp_raw) <- avgexp_cluster_ids
    
    query_avgexp <- as.data.frame(
      t(query_avgexp_raw),
      check.names = FALSE
    )
    
    varfeat <- unique(trimws(readLines(varfeat_file)))
    varfeat <- varfeat[nzchar(varfeat)]
    
    common_varfeats <- intersect(
      varfeat,
      db_varfeats
    )
    
    if (length(common_varfeats) == 0) {
      stop(
        "No variable-feature overlap for ",
        cell_type,
        " res ",
        res_str
      )
    }
    
    ann <- Create_SAHA_object(
      query = query_markers,
      db = db_markers,
      data_type = "Markers"
    )
    
    ann <- Create_SAHA_object(
      query = query_avgexp,
      db = db_avgexp,
      data_type = "AvgExp",
      existing = ann
    )
    
    ann <- Initialize_Markers(
      ann,
      spec_thresh = 1.1
    )
    
    ann <- Tune_Markers(
      ann = ann,
      method = "absolute",
      method_value = 100,
      method_var = "avg_log2FC",
      set = "db"
    )
    
    marker_richness <- capture.output(
      Marker_Richness(ann)
    ) %>%
      read.table(
        text = .,
        header = TRUE
      )
    
    if (
      "gene" %in% colnames(marker_richness) &&
      any(!is.na(marker_richness$gene)) &&
      median(marker_richness$gene, na.rm = TRUE) > 500
    ) {
      ann <- Tune_Markers(
        ann = ann,
        method = "absolute",
        method_value = 500,
        method_var = "avg_log2FC",
        set = "query"
      )
    }
    
    ann <- Run_Marker_Based(ann)
    ann <- Create_MarkerBased_Viz(ann)
    
    ann <- Initialize_MarkerFree(
      ann = ann
    )
    
    ann <- Downsample(
      ann,
      custom_ds = common_varfeats
    )
    
    ann <- NormalizeDS(
      ann,
      assay_query = "RNA"
    )
    
    ann <- CorrelateDS(ann)
    ann <- Create_MarkerFree_Viz(ann)
    
    # Same correction as the working primary script, but match a literal
    # leading "RNA." rather than treating "." as a regex wildcard.
    colnames(ann@results$marker_free$corr) <- sub(
      "^RNA\\.",
      "",
      colnames(ann@results$marker_free$corr)
    )
    
    auto <- AutoAnnotate(
      ann,
      data_type = "Both"
    )
    
    Generate_SAHA_Report(
      ann,
      auto = auto,
      output_file = html_file
    )
    
    # Preserve save() because the working primary script explicitly used it
    # after RDS serialization of the SAHA object had caused problems.
    save(
      ann,
      file = ann_file
    )
    
    write.csv(
      auto,
      auto_file,
      row.names = FALSE
    )
    
    rm(
      query_markers_raw,
      query_markers,
      finite_fc,
      max_fc,
      min_fc
    )
    
    rm(
      query_avgexp_raw,
      query_avgexp,
      avgexp_cluster_ids,
      varfeat,
      common_varfeats,
      marker_richness,
      ann,
      auto
    )
  }
  
  rm(
    db_markers,
    db_avgexp,
    db_varfeats
  )
  
  completed_cell_types <- c(
    completed_cell_types,
    cell_type
  )
}

if (length(completed_cell_types) == 0) {
  stop("No cell compartments completed SAHA annotation.")
}

# ==============================================================================
# 2. Attach SAHA calls to cell-level resolution assignments and summarize
# ==============================================================================

for (cell_type in completed_cell_types) {
  message(
    "\n========== Resolution summaries: ",
    cell_type,
    " =========="
  )
  
  metadata_file <- file.path(
    input_dir,
    paste0(
      "seaad_ec_",
      cell_type,
      "_resolution_assignments_obs.csv"
    )
  )
  
  if (!file.exists(metadata_file)) {
    stop("Missing metadata: ", metadata_file)
  }
  
  cell_metadata <- read.csv(
    metadata_file,
    check.names = FALSE
  )
  
  # Resolve raw Leiden columns before adding standardized annotation columns.
  leiden_cols <- setNames(
    rep(
      NA_character_,
      length(resolution_strings)
    ),
    resolution_strings
  )
  
  for (res_i in seq_along(resolution_strings)) {
    res_str <- resolution_strings[res_i]
    
    candidates <- c(
      paste0("leiden_", res_str),
      paste0("leiden_scvi_", res_str),
      paste0("leiden_res", res_str),
      paste0("leiden_res_", res_str)
    )
    
    hits <- candidates[
      candidates %in% colnames(cell_metadata)
    ]
    
    if (length(hits) == 0) {
      res_pattern <- gsub(
        "\\.",
        "\\\\.",
        res_str
      )
      
      hits <- grep(
        paste0(
          "^leiden.*",
          res_pattern,
          "$"
        ),
        colnames(cell_metadata),
        value = TRUE
      )
    }
    
    if (length(hits) != 1) {
      stop(
        "Could not uniquely resolve Leiden column for ",
        cell_type,
        " res ",
        res_str,
        ". Hits: ",
        paste(hits, collapse = ", ")
      )
    }
    
    leiden_cols[[res_str]] <- hits[1]
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
  
  donor_hits <- donor_candidates[
    donor_candidates %in% colnames(cell_metadata)
  ]
  
  donor_col <- if (length(donor_hits) > 0) {
    donor_hits[1]
  } else {
    NA_character_
  }
  
  age_candidates <- c(
    "age",
    "Age",
    "age_at_death",
    "Age at Death",
    "age_death",
    "ageDeath"
  )
  
  age_hits <- age_candidates[
    age_candidates %in% colnames(cell_metadata)
  ]
  
  age_col <- if (length(age_hits) > 0) {
    age_hits[1]
  } else {
    NA_character_
  }
  
  if (!is.na(donor_col)) {
    message("Donor column: ", donor_col)
  } else {
    message(
      "No donor column detected; donor counts will be omitted."
    )
  }
  
  if (!is.na(age_col)) {
    message("Age column: ", age_col)
  } else {
    message(
      "No age column detected; age/cell-count plot will be omitted."
    )
  }
  
  lookup_list <- list()
  cluster_summary_list <- list()
  
  for (res_i in seq_along(resolution_values)) {
    res <- resolution_values[res_i]
    res_str <- resolution_strings[res_i]
    leiden_col <- leiden_cols[[res_str]]
    
    auto_file <- file.path(
      saha_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_res",
        res_str,
        "_SAHAvsECatlas_auto.csv"
      )
    )
    
    if (!file.exists(auto_file)) {
      stop("Missing SAHA auto file: ", auto_file)
    }
    
    auto_df <- read.csv(
      auto_file,
      check.names = FALSE
    )
    
    missing_auto_columns <- setdiff(
      c(
        "cluster",
        "marker_free"
      ),
      colnames(auto_df)
    )
    
    if (length(missing_auto_columns) > 0) {
      stop(
        "Unexpected SAHA auto format in ",
        auto_file
      )
    }
    
    lookup_df <- auto_df %>%
      transmute(
        cluster_char = sub(
          "^RNA\\.g",
          "",
          as.character(cluster)
        ),
        annotation_full = as.character(marker_free),
        annotation_class = word(
          annotation_full,
          1
        ),
        annotation_subtype = str_extract(
          annotation_full,
          "^\\S+(?:\\s+\\S+)?"
        ),
        annotation_marker = str_extract(
          annotation_full,
          "^\\S+(?:\\s+\\S+){0,2}"
        )
      )
    
    if (anyDuplicated(lookup_df$cluster_char) > 0) {
      stop(
        "Duplicate SAHA cluster calls in ",
        auto_file
      )
    }
    
    metadata_clusters <- unique(
      sub(
        "^RNA\\.g",
        "",
        as.character(
          cell_metadata[[leiden_col]]
        )
      )
    )
    
    metadata_clusters <- metadata_clusters[
      !is.na(metadata_clusters)
    ]
    
    unmatched_clusters <- setdiff(
      metadata_clusters,
      lookup_df$cluster_char
    )
    
    if (length(unmatched_clusters) > 0) {
      stop(
        "Unmatched clusters for ",
        cell_type,
        " res ",
        res_str,
        ": ",
        paste(
          unmatched_clusters,
          collapse = ", "
        )
      )
    }
    
    full_col <- paste0(
      "leiden_",
      res_str,
      "_annotation"
    )
    
    class_col <- paste0(
      "leiden_",
      res_str,
      "_class"
    )
    
    subtype_col <- paste0(
      "leiden_",
      res_str,
      "_subtype"
    )
    
    marker_col <- paste0(
      "leiden_",
      res_str,
      "_marker"
    )
    
    cell_metadata <- cell_metadata %>%
      mutate(
        .cluster_join = sub(
          "^RNA\\.g",
          "",
          as.character(
            .data[[leiden_col]]
          )
        )
      ) %>%
      left_join(
        lookup_df,
        by = c(
          ".cluster_join" = "cluster_char"
        )
      ) %>%
      rename(
        !!full_col := annotation_full,
        !!class_col := annotation_class,
        !!subtype_col := annotation_subtype,
        !!marker_col := annotation_marker
      ) %>%
      dplyr::select(
        -.cluster_join
      )
    
    if (any(is.na(cell_metadata[[full_col]]))) {
      stop(
        "Missing annotation after join for ",
        cell_type,
        " res ",
        res_str
      )
    }
    
    lookup_list[[res_str]] <- lookup_df %>%
      mutate(
        cell_type = cell_type,
        resolution = res
      )
    
    cluster_counts <- cell_metadata %>%
      mutate(
        cluster_char = sub(
          "^RNA\\.g",
          "",
          as.character(
            .data[[leiden_col]]
          )
        )
      ) %>%
      count(
        cluster_char,
        name = "n_cells"
      )
    
    if (!is.na(donor_col)) {
      donor_counts <- cell_metadata %>%
        mutate(
          cluster_char = sub(
            "^RNA\\.g",
            "",
            as.character(
              .data[[leiden_col]]
            )
          )
        ) %>%
        group_by(
          cluster_char
        ) %>%
        summarize(
          n_donors = n_distinct(
            .data[[donor_col]],
            na.rm = TRUE
          ),
          .groups = "drop"
        )
      
      cluster_counts <- cluster_counts %>%
        left_join(
          donor_counts,
          by = "cluster_char"
        )
    } else {
      cluster_counts <- cluster_counts %>%
        mutate(
          n_donors = NA_integer_
        )
    }
    
    if (
      !is.na(age_col) &&
      is.numeric(cell_metadata[[age_col]])
    ) {
      age_summary <- cell_metadata %>%
        mutate(
          cluster_char = sub(
            "^RNA\\.g",
            "",
            as.character(
              .data[[leiden_col]]
            )
          )
        ) %>%
        group_by(
          cluster_char
        ) %>%
        summarize(
          min_age = ifelse(
            all(
              is.na(
                .data[[age_col]]
              )
            ),
            NA_real_,
            min(
              .data[[age_col]],
              na.rm = TRUE
            )
          ),
          max_age = ifelse(
            all(
              is.na(
                .data[[age_col]]
              )
            ),
            NA_real_,
            max(
              .data[[age_col]],
              na.rm = TRUE
            )
          ),
          .groups = "drop"
        )
      
      cluster_counts <- cluster_counts %>%
        left_join(
          age_summary,
          by = "cluster_char"
        )
    } else {
      cluster_counts <- cluster_counts %>%
        mutate(
          min_age = NA_real_,
          max_age = NA_real_
        )
    }
    
    cluster_summary_list[[res_str]] <- cluster_counts %>%
      left_join(
        lookup_df,
        by = "cluster_char"
      ) %>%
      mutate(
        cell_type = cell_type,
        resolution = res,
        leiden_column = leiden_col
      )
  }
  
  lookup_all <- bind_rows(
    lookup_list
  )
  
  cluster_summary <- bind_rows(
    cluster_summary_list
  )
  
  write.csv(
    cell_metadata,
    file.path(
      table_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_resolution_assignments_annotated.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    lookup_all,
    file.path(
      table_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_SAHA_annotation_lookup_all_resolutions.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    cluster_summary,
    file.path(
      table_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_resolution_cluster_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Across-resolution summaries
  # ---------------------------------------------------------------------------
  
  if (!is.na(donor_col)) {
    marker_summary <- cluster_summary %>%
      group_by(
        resolution,
        annotation_marker,
        annotation_subtype
      ) %>%
      summarize(
        total_cells = sum(n_cells),
        max_donors = max(n_donors),
        n_clusters = n_distinct(cluster_char),
        .groups = "drop"
      ) %>%
      mutate(
        plot_label = paste0(
          "D=",
          max_donors,
          "\nK=",
          n_clusters
        )
      )
    
    subtype_summary <- cluster_summary %>%
      group_by(
        resolution,
        annotation_subtype
      ) %>%
      summarize(
        total_cells = sum(n_cells),
        max_donors = max(n_donors),
        n_clusters = n_distinct(cluster_char),
        .groups = "drop"
      ) %>%
      mutate(
        plot_label = paste0(
          "D=",
          max_donors,
          "\nK=",
          n_clusters
        )
      )
  } else {
    marker_summary <- cluster_summary %>%
      group_by(
        resolution,
        annotation_marker,
        annotation_subtype
      ) %>%
      summarize(
        total_cells = sum(n_cells),
        n_clusters = n_distinct(cluster_char),
        .groups = "drop"
      ) %>%
      mutate(
        max_donors = NA_integer_,
        plot_label = paste0(
          "K=",
          n_clusters
        )
      )
    
    subtype_summary <- cluster_summary %>%
      group_by(
        resolution,
        annotation_subtype
      ) %>%
      summarize(
        total_cells = sum(n_cells),
        n_clusters = n_distinct(cluster_char),
        .groups = "drop"
      ) %>%
      mutate(
        max_donors = NA_integer_,
        plot_label = paste0(
          "K=",
          n_clusters
        )
      )
  }
  
  write.csv(
    marker_summary,
    file.path(
      table_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_across_resolution_marker_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    subtype_summary,
    file.path(
      table_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_across_resolution_subtype_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  p_marker <- ggplot(
    marker_summary,
    aes(
      x = annotation_marker,
      y = total_cells,
      fill = annotation_subtype
    )
  ) +
    geom_col() +
    geom_text(
      aes(
        label = plot_label
      ),
      vjust = -0.2,
      size = 2.5
    ) +
    facet_wrap(
      ~resolution,
      nrow = 3,
      scales = "free_x"
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0,
          0.25
        )
      )
    ) +
    labs(
      x = NULL,
      y = "Cells",
      title = paste0(
        cell_type,
        ": cells per SAHA marker annotation"
      )
    ) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      legend.position = "none"
    )
  
  ggsave(
    file.path(
      plot_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_AcrossRes_CellsPerMarker.pdf"
      )
    ),
    p_marker,
    width = 12,
    height = 9
  )
  
  p_subtype <- ggplot(
    subtype_summary,
    aes(
      x = annotation_subtype,
      y = total_cells,
      fill = annotation_subtype
    )
  ) +
    geom_col() +
    geom_text(
      aes(
        label = plot_label
      ),
      vjust = -0.2,
      size = 2.7
    ) +
    facet_wrap(
      ~resolution,
      nrow = 3,
      scales = "free_x"
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0,
          0.25
        )
      )
    ) +
    labs(
      x = NULL,
      y = "Cells",
      title = paste0(
        cell_type,
        ": cells per SAHA subtype annotation"
      )
    ) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      legend.position = "none"
    )
  
  ggsave(
    file.path(
      plot_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_AcrossRes_CellsPerSubtype.pdf"
      )
    ),
    p_subtype,
    width = 12,
    height = 9
  )
  
  # ---------------------------------------------------------------------------
  # Clustrees on SAHA-derived subtype and marker labels
  # ---------------------------------------------------------------------------
  
  subtype_columns <- paste0(
    "leiden_",
    resolution_strings,
    "_subtype"
  )
  
  marker_columns <- paste0(
    "leiden_",
    resolution_strings,
    "_marker"
  )
  
  subtype_df <- cell_metadata %>%
    dplyr::select(
      all_of(subtype_columns)
    )
  
  marker_df <- cell_metadata %>%
    dplyr::select(
      all_of(marker_columns)
    )
  
  png(
    file.path(
      plot_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_clustree_subtype.png"
      )
    ),
    width = 15,
    height = 10,
    units = "in",
    res = 300
  )
  
  print(
    clustree(
      subtype_df,
      prefix = "leiden_",
      suffix = "_subtype"
    )
  )
  
  dev.off()
  
  png(
    file.path(
      plot_dir,
      paste0(
        "seaad_ec_",
        cell_type,
        "_clustree_marker.png"
      )
    ),
    width = 20,
    height = 10,
    units = "in",
    res = 300
  )
  
  print(
    clustree(
      marker_df,
      prefix = "leiden_",
      suffix = "_marker"
    )
  )
  
  dev.off()
  
  # ---------------------------------------------------------------------------
  # Optional final-resolution export after you choose the resolution
  # ---------------------------------------------------------------------------
  
  selected_res <- selected_resolution[[cell_type]]
  
  if (!is.na(selected_res)) {
    selected_res_str <- sprintf(
      "%.1f",
      selected_res
    )
    
    level_oi <- selected_level[[cell_type]]
    
    if (!level_oi %in% c("subtype", "marker")) {
      stop(
        "selected_level must be subtype or marker for ",
        cell_type
      )
    }
    
    selected_annotation_col <- paste0(
      "leiden_",
      selected_res_str,
      "_",
      level_oi
    )
    
    if (!selected_annotation_col %in% colnames(cell_metadata)) {
      stop(
        "Missing selected annotation column: ",
        selected_annotation_col
      )
    }
    
    final_metadata <- cell_metadata %>%
      mutate(
        final_anno = .data[[selected_annotation_col]]
      )
    
    write.csv(
      final_metadata,
      file.path(
        table_dir,
        paste0(
          "seaad_ec_",
          cell_type,
          "_final_res",
          selected_res_str,
          "_",
          level_oi,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    if (
      !is.na(donor_col) &&
      !is.na(age_col) &&
      is.numeric(final_metadata[[age_col]])
    ) {
      donor_cell_summary <- final_metadata %>%
        group_by(
          final_anno,
          donor = .data[[donor_col]]
        ) %>%
        summarize(
          n_cell = n(),
          age = mean(
            .data[[age_col]],
            na.rm = TRUE
          ),
          .groups = "drop"
        )
      
      p_power <- ggplot(
        donor_cell_summary,
        aes(
          x = age,
          y = n_cell
        )
      ) +
        geom_point() +
        scale_y_log10() +
        geom_hline(
          yintercept = 10,
          linetype = "dashed"
        ) +
        facet_wrap(
          ~final_anno
        ) +
        labs(
          x = "Age",
          y = "Cells per donor",
          title = paste0(
            cell_type,
            " res ",
            selected_res_str,
            " ",
            level_oi
          )
        ) +
        theme_linedraw()
      
      ggsave(
        file.path(
          plot_dir,
          paste0(
            "seaad_ec_",
            cell_type,
            "_nDonor_nCells_res",
            selected_res_str,
            "_",
            level_oi,
            ".pdf"
          )
        ),
        p_power,
        width = 10,
        height = 6
      )
    }
  }
}

capture.output(
  sessionInfo(),
  file = file.path(
    output_root,
    "SEAAD_EC_resolution_selection_sessionInfo.txt"
  )
)

message("\n========== COMPLETE ==========")
message(
  "Completed compartments: ",
  paste(
    completed_cell_types,
    collapse = ", "
  )
)

if (length(setdiff(cell_types, completed_cell_types)) > 0) {
  message(
    "Not completed: ",
    paste(
      setdiff(
        cell_types,
        completed_cell_types
      ),
      collapse = ", "
    )
  )
}

message(
  "Results: ",
  output_root
)

message(
  "Plots:   ",
  plot_dir
)

