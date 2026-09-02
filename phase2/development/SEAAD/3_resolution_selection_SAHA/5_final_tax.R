
##didn't end up running this because the tax after 4_ should be all Raph needs?

# # ==============================================================================
# # final_5_SEAAD_final_annotation_summary.R
# #
# # Brief summary of final SEA-AD annotations after resolution selection
# #
# # Final resolutions:
# #   ExN         = 0.3
# #   InN         = 0.4
# #   NonNeuronal = 0.2
# #
# # Shows:
# #   1. Where the final taxonomy files are
# #   2. Cell counts / proportions by anno
# #   3. Cell counts / proportions by anno_coarse
# #   4. anno -> anno_coarse mapping
# # ==============================================================================
# 
# library(tidyverse)
# 
# 
# # ==============================================================================
# # 1. Paths
# # ==============================================================================
# 
# final_tax_dir <- paste0(
#   "/data/ADRD/brain_aging/exploration/res/Phase2/SEAAD_EC/",
#   "resolution_selection/final_annotation/tables/"
# )
# 
# ExN_file <- file.path(
#   final_tax_dir,
#   "seaad_ec_ExN_selected_res0.3_final_taxonomy.csv"
# )
# 
# InN_file <- file.path(
#   final_tax_dir,
#   "seaad_ec_InN_selected_res0.4_final_taxonomy.csv"
# )
# 
# NonNeuronal_file <- file.path(
#   final_tax_dir,
#   "seaad_ec_NonNeuronal_selected_res0.2_final_taxonomy.csv"
# )
# 
# message("\n============================================================")
# message("FINAL TAXONOMY FILES")
# message("============================================================")
# message("ExN:         ", ExN_file)
# message("InN:         ", InN_file)
# message("NonNeuronal: ", NonNeuronal_file)
# 
# 
# # ==============================================================================
# # 2. Read final taxonomies
# # ==============================================================================
# 
# ExN_tax <- read.csv(ExN_file, check.names = FALSE) %>%
#   mutate(cell_type = "ExN")
# 
# InN_tax <- read.csv(InN_file, check.names = FALSE) %>%
#   mutate(cell_type = "InN")
# 
# NonNeuronal_tax <- read.csv(NonNeuronal_file, check.names = FALSE) %>%
#   mutate(cell_type = "NonNeuronal")
# 
# final_tax <- bind_rows(
#   ExN_tax,
#   InN_tax,
#   NonNeuronal_tax
# )
# 
# 
# # ==============================================================================
# # 3. Final anno summary
# # ==============================================================================
# 
# anno_summary <- final_tax %>%
#   count(
#     cell_type,
#     anno,
#     name = "n_cells"
#   ) %>%
#   group_by(cell_type) %>%
#   mutate(
#     pct_cells = 100 * n_cells / sum(n_cells)
#   ) %>%
#   ungroup() %>%
#   arrange(
#     cell_type,
#     desc(n_cells)
#   )
# 
# message("\n============================================================")
# message("FINAL anno SUMMARY")
# message("============================================================")
# 
# print(anno_summary)
# 
# 
# # ==============================================================================
# # 4. Final anno_coarse summary
# # ==============================================================================
# 
# anno_coarse_summary <- final_tax %>%
#   count(
#     cell_type,
#     anno_coarse,
#     name = "n_cells"
#   ) %>%
#   group_by(cell_type) %>%
#   mutate(
#     pct_cells = 100 * n_cells / sum(n_cells)
#   ) %>%
#   ungroup() %>%
#   arrange(
#     cell_type,
#     desc(n_cells)
#   )
# 
# message("\n============================================================")
# message("FINAL anno_coarse SUMMARY")
# message("============================================================")
# 
# print(anno_coarse_summary)
# 
# 
# # ==============================================================================
# # 5. Show anno -> anno_coarse hierarchy
# # ==============================================================================
# 
# anno_lookup <- final_tax %>%
#   distinct(
#     cell_type,
#     anno,
#     anno_coarse
#   ) %>%
#   arrange(
#     cell_type,
#     anno_coarse,
#     anno
#   )
# 
# message("\n============================================================")
# message("FINAL anno -> anno_coarse LOOKUP")
# message("============================================================")
# 
# print(anno_lookup)
# 
# 
# # ==============================================================================
# # 6. Save summaries
# # ==============================================================================
# 
# write.csv(
#   anno_summary,
#   file.path(final_tax_dir, "seaad_ec_final_anno_summary.csv"),
#   row.names = FALSE
# )
# 
# write.csv(
#   anno_coarse_summary,
#   file.path(final_tax_dir, "seaad_ec_final_anno_coarse_summary.csv"),
#   row.names = FALSE
# )
# 
# write.csv(
#   anno_lookup,
#   file.path(final_tax_dir, "seaad_ec_final_anno_to_anno_coarse_lookup.csv"),
#   row.names = FALSE
# )
# 
# 
# message("\n============================================================")
# message("FINAL ANNOTATION SUMMARY COMPLETE")
# message("============================================================")
# message("Summary files written to:")
# message(final_tax_dir)