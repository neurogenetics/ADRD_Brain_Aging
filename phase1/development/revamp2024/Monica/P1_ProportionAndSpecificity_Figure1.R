#Phase 1: Cell type proportion and regional specificity analysis and heatmaps
#Code written and conceptualized by Dom. J. Acri with support from LLM.

#Load Necessary Packages
library(Seurat)
library(speckle)
library(circlize)
library(tidyverse)
library(dplyr)
library(purrr)
library(ComplexHeatmap)

#Set Working Directory
setwd("/data/ADRD/brain_aging/exploration/")

# Read in obj with clustered cells
obj=readRDS("./aging.pegasus.leiden_085.subclustered.rds")


# cells = data frame of each cell with associated cell type, sample ID, and meta data
cells = obj@meta.data

cells$broad_celltype = factor(cells$broad_celltype,levels = unique(cells$broad_celltype))

cells$Brain_region

x=propeller(clusters = cells$broad_celltype, sample = cells$Sample_id, group = cells$Age_group )

# Get unique brain regions
unique_regions <- unique(cells$Brain_region)

# Function to run propeller for a specific region
run_propeller_for_region <- function(region) {
  # Subset the data for the current region
  region_cells <- cells[cells$Brain_region == region, ]
  
  # Run propeller
  result <- propeller(clusters = region_cells$broad_celltype, 
                      sample = region_cells$Sample_id, 
                      group = region_cells$Age_group)
  
  # Add region column to the result
  result$Region <- region
  
  return(result)
}

# Apply the function to each region and combine results
x <- map_dfr(unique_regions, run_propeller_for_region)

# If needed, you can reorder the columns to have Region as the first column
x <- x %>% select(Region, everything())

#since Prop here is not relative within each Region * CT, need to normalize; 1 = all old, -1 = all young
x$age_effect <- (x$PropMean.old - x$PropMean.young)/(x$PropMean.old + x$PropMean.young)

x$FDR = p.adjust(x$P.Value, method = "BH")

## add in total cell count so I can do specificty calc in same df:
x <- x %>%
  rowwise() %>%
  mutate(total_cell = sum(cells$Brain_region == Region & 
                            cells$broad_celltype == BaselineProp.clusters, 
                          na.rm = TRUE)) %>%
  ungroup()%>%
  data.frame()

#Calculate region specificity
x <- x %>%
  group_by(BaselineProp.clusters) %>%
  mutate(total_cells_per_type = sum(total_cell)) %>%
  ungroup() %>%
  mutate(region_spec = total_cell / total_cells_per_type) %>%
  ungroup()%>%
  data.frame()

hm_input <- x %>%
  mutate(ct=BaselineProp.clusters, age_FDR=FDR)%>%  
  select(Region,ct,age_FDR,age_effect,region_spec)

# Assuming hm_input_complete is the dataframe with all combinations

# Prepare the matrices for heatmaps
mat1 <- hm_input %>%
  select(-age_FDR,-age_effect)%>%
  pivot_wider(names_from = ct, values_from = region_spec) %>%
  column_to_rownames("Region") %>%
  as.matrix()

mat2 <- hm_input %>%
  select(-age_FDR,-region_spec)%>%
  pivot_wider(names_from = ct, values_from = age_effect) %>%
  column_to_rownames("Region") %>%
  as.matrix()

# Create the heatmaps
ht1 <- Heatmap(mat1,
               col = colorRamp2(c(0, 1), c("white", "green")),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               show_column_names = FALSE,
               heatmap_legend_param = list(title = "region_spec"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA))
               })

ht2 <- Heatmap(mat2,
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               show_column_names = TRUE,
               heatmap_legend_param = list(title = "age_effect"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA))
                 grid.text(ifelse(hm_input$age_FDR[i] < 0.05, "*", "ns"),
                           x = x, y = y, gp = gpar(fontsize = 8))
               })

ht_list = ht1 %v% ht2
draw(ht_list)

# # Display the combined heatmap
# draw(ht_list, heatmap_legend_side = "bottom")

