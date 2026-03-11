library(fgsea)
library(tidyverse)

setwd("/data/ADRD/brain_aging/")
out_dir <- "./exploration/revision2_senpaths/B_marker/"

# 1. Load Pathways (Assuming you saved your sen_paths list to RDS)
# sen_paths <- readRDS("./0_curation_list/senescence_pathways.rds")
# For now, using your example genes
sen_paths <- readRDS("./exploration/revision2_senpaths/SenPaths_completed_final.rds")
sen_paths <- sen_paths %>%
  keep(~ !is.null(.x) && length(.x) > 0 && any(.x != ""))
names(sen_paths)

# 2. Function to run fgsea on the lists
run_fgsea_on_list <- function(marker_list) {
  map_df(names(marker_list), function(nm) {
    df <- marker_list[[nm]]
    # Rank by -log10(p) * sign(FC)
    stats <- -log10(df$p_val + 1e-300) * sign(df$avg_log2FC)
    names(stats) <- df$gene
    
    fgseares <- fgsea(pathways = sen_paths, stats = stats, minSize = 1)
    fgseares$test_name <- nm
    return(fgseares)
  })
}

# 3. Process both Global and Regional
global_markers <- readRDS(paste0(out_dir, "markers_list_perCT.rds"))
res_global <- run_fgsea_on_list(global_markers)
saveRDS(res_global, paste0(out_dir, "fgsea_res_perCT.rds"))

regional_markers <- readRDS(paste0(out_dir, "markers_list_perCTperRegion.rds"))
res_regional <- run_fgsea_on_list(regional_markers)
saveRDS(res_regional, paste0(out_dir, "fgsea_res_perCTperRegion.rds"))
