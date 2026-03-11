library(ComplexHeatmap)
library(circlize)
library(tidyverse)

setwd("/data/ADRD/brain_aging/")
out_dir <- "./exploration/revision2_senpaths/C_PCA/"
res <- readRDS(paste0(out_dir, "glm_results_summary_final.rds"))

# Prepare Matrices
# Rows = test (CT or Region_CT), Cols = pathways
beta_mat <- res %>% 
  select(test, pathway, Estimate) %>% 
  pivot_wider(names_from = pathway, values_from = Estimate) %>% 
  column_to_rownames("test")

fdr_mat <- res %>% 
  select(test, pathway, fdr) %>% 
  pivot_wider(names_from = pathway, values_from = fdr) %>% 
  column_to_rownames("test")

# Ensure they match order
fdr_mat <- fdr_mat[rownames(beta_mat), colnames(beta_mat)]

# Significance Stars
star_mat <- matrix("", nrow = nrow(fdr_mat), ncol = ncol(fdr_mat))
star_mat[fdr_mat < 0.05] <- "﹡"
star_mat[fdr_mat < 0.01] <- "﹡﹡"
star_mat[fdr_mat < 0.001] <- "﹡﹡﹡"

# Metadata for row splitting (Global vs Regional)
rownames(res)<-NULL
row_info <- res %>% dplyr::select(test, type) %>% distinct() %>% column_to_rownames("test")
row_split_vec <- row_info[rownames(beta_mat), "type"]

# Define color scale (centered at 0)
max_val <- max(abs(beta_mat), na.rm = TRUE)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Plot
h <- Heatmap(as.matrix(beta_mat), 
             name = "Age Beta (PC1)",
             col = col_fun,
             row_split = row_split_vec,
             column_title = "Senescence Pathway PCA (Age Association)",
             cluster_columns = TRUE,
             show_row_dend = FALSE,
             rect_gp = gpar(col = "white", lwd = 1),
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(star_mat[i, j], x, y, gp = gpar(fontsize = 12))
             })

# http://localhost:35241/graphics/plot_zoom?width=627&height=808&scale=1
cairo_pdf(paste0(out_dir, "senescence_pathway_PCA_heatmap_final.pdf"), width = 6.27, height = 8.08)
draw(h, row_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

png(paste0(out_dir, "senescence_pathway_PCA_heatmap_final.png"), width = 6.27, height = 8.08,units = "in",res = 600)
draw(h, row_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

