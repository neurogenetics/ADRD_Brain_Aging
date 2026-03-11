library(tidyverse)
library(Seurat)

setwd("/data/ADRD/brain_aging/")
out_dir <- "./exploration/revision2_senpaths/C_PCA/"

# Load processed metadata from 01_pbulk.R
meta <- readRDS(paste0(out_dir, "metadata_with_pb_keys.rds"))

# Define Senescence Pathways
sen_paths <- readRDS("./exploration/revision2_senpaths/SenPaths_completed_final.rds")
sen_paths <- sen_paths %>%
  keep(~ !is.null(.x) && length(.x) > 0 && any(.x != ""))

# Helper function to run PCA and GLMs per pathway
run_pathway_pca_glm <- function(mat, type = c("global", "regional"), pathways) {
  type <- match.arg(type)
  pb_col <- ifelse(type == "global", "pb_global", "pb_regional")
  
  pathway_results <- list()
  
  for (p_name in names(pathways)) {
    genes <- intersect(pathways[[p_name]], rownames(mat))
    if (length(genes) < 3) next
    
    # 1. PCA on Pathway-specific genes
    pca_res <- prcomp(t(mat[genes, ]), scale. = TRUE)
    pca_df <- as.data.frame(pca_res$x[, 1:2]) 
    
    # Standardize IDs to underscores
    pca_df[[pb_col]] <- gsub("-", "_", rownames(pca_df))
    
    # 2. Join with Metadata (Ensure Age_group is selected)
    pb_meta <- meta %>%
      dplyr::select(!!sym(pb_col), donor_id, Brain_region, broad_celltype, Age_group, Sex) %>%
      distinct() %>%
      mutate(!!sym(pb_col) := gsub("-", "_", !!sym(pb_col))) %>%
      inner_join(pca_df, by = pb_col)
    
    if(nrow(pb_meta) == 0) next
    
    # CRITICAL: Relevel Age_group so "young" is the reference
    # This makes the "Age_groupold" coefficient represent Old vs Young
    pb_meta$Age_group <- factor(pb_meta$Age_group, levels = c("young", "old"))
    
    # 3. GLM Loop per Celltype/Region
    group_factor <- if(type == "global") "broad_celltype" else c("broad_celltype", "Brain_region")
    groups <- pb_meta %>% select(all_of(group_factor)) %>% distinct()
    
    for(i in 1:nrow(groups)) {
      sub_df <- inner_join(pb_meta, groups[i,,drop=F], by = group_factor)
      if(nrow(sub_df) < 10) next
      
      # Correctly catch names for labels
      test_label <- paste(as.character(unlist(groups[i,])), collapse="_")
      
      try({
        # Model PC1 against Age_group and Sex
        fit <- glm(PC1 ~ Age_group + Sex, data = sub_df)
        stats <- as.data.frame(summary(fit)$coefficients)
        
        # We look for "Age_groupold" because R appends the non-reference level name
        age_row <- "Age_groupold"
        
        if(age_row %in% rownames(stats)) {
          res <- stats[age_row, ]
          res$test <- test_label
          res$pathway <- p_name
          res$type <- type
          pathway_results[[paste(p_name, test_label, sep=":")]] <- res
        }
      }, silent = TRUE)
    }
  }
  return(bind_rows(pathway_results))
}

# --- Execution ---

message("Running Global PCA-GLMs...")
mat_glob <- readRDS(paste0(out_dir, "pbulk_matrix_global.rds"))
res_glob <- run_pathway_pca_glm(mat_glob, "global", sen_paths)

message("Running Regional PCA-GLMs...")
mat_reg <- readRDS(paste0(out_dir, "pbulk_matrix_regional.rds"))
res_reg <- run_pathway_pca_glm(mat_reg, "regional", sen_paths)

# Combine and apply FDR correction
final_results <- bind_rows(res_glob, res_reg)
final_results$fdr <- p.adjust(final_results$`Pr(>|t|)`, method = "BH")

# Save results
saveRDS(final_results, paste0(out_dir, "glm_results_summary_final.rds"))
