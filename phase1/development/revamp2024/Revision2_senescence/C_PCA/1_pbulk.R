library(Seurat)
library(tidyverse)
library(edgeR)

setwd("/data/ADRD/brain_aging/")
out_dir <- "./exploration/revision2_senpaths/C_PCA/"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

obj <- readRDS("./exploration/aging.pegasus.leiden_085.subclustered.rds")

# Function to pbulk and normalize
get_pbulk_matrix <- function(seurat_obj, group_var) {
  # Filter for samples with >= 50 cells
  keep_samples <- seurat_obj@meta.data %>%
    group_by(!!sym(group_var)) %>%
    tally() %>% filter(n >= 50) %>% pull(!!sym(group_var))
  
  subset_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_var]] %in% keep_samples, ])
  subset_obj <- subset(seurat_obj, cells = subset_cells)
  
  # Aggregate (Seurat converts _ to - here)
  counts <- AggregateExpression(subset_obj, group.by = group_var, assays = "RNA", slot = "counts")$RNA
  
  # FIX: Convert dashes back to underscores to match metadata
  colnames(counts) <- gsub("-", "_", colnames(counts))
  
  # TMM Normalization
  dge <- DGEList(counts)
  dge <- calcNormFactors(dge, method = "TMM")
  return(edgeR::cpm(dge, log = TRUE))
}

# Global: donor_id + Celltype
obj$pb_global <- paste(obj$donor_id, obj$broad_celltype, sep = "___")
mat_global <- get_pbulk_matrix(obj, "pb_global")
saveRDS(mat_global, paste0(out_dir, "pbulk_matrix_global.rds"))

# Regional: donor_id + Region + Celltype
obj$pb_regional <- paste(obj$donor_id, obj$Brain_region, obj$broad_celltype, sep = "___")
mat_regional <- get_pbulk_matrix(obj, "pb_regional")
saveRDS(mat_regional, paste0(out_dir, "pbulk_matrix_regional.rds"))

# Save the updated meta for the next step
saveRDS(obj@meta.data, paste0(out_dir, "metadata_with_pb_keys.rds"))
