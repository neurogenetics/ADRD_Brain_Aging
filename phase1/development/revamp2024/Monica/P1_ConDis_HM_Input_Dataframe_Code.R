#Function to Create Input Dataframes for Concordance/Discordance Heatmaps
#Written/conceptualized by Monica Mesecar and Dom J. Acri with assistance from Perplexity AI for efficiency
install.packages("openxlsx")

# Libraries
library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(openxlsx)


# List of input files and corresponding cell type labels
file_list <- c("Inhib_ByRegion_aDEG_fixed.csv", "ExN_ByRegion_aDEG_fixed.csv","OD_ByRegion_aDEG_fixed.csv", "OPC_ByRegion_aDEG_fixed.csv","Mural_ByRegion_aDEG_fixed.csv", "Astro_ByRegion_aDEG_fixed.csv","Microglia_ByRegion_aDEG_fixed.csv", "Epend_ByRegion_aDEG_fixed.csv", "Endo_ByRegion_aDEG_fixed.csv") # Add more as needed
cell_types <- c("InN", "ExN", "Oligo", "OPC","Mural", "Astro", "Micro", "Epend", "Endo") # Should match file_list order

# Function to process each file
process_celltype <- function(file, celltype) {
  df <- read.csv(file)
  regions <- colnames(df)[2:5]
  combinations <- combn(regions, 2, simplify = FALSE)
  
  classify_combination <- function(val1, val2) {
    if (val1 == "Up" && val2 == "Up") {
      return("uu")
    } else if (val1 == "Down" && val2 == "Down") {
      return("dd")
    } else if (val1 == "Up" && val2 == "Down") {
      return("ud")
    } else if (val1 == "Down" && val2 == "Up") {
      return("du")
    } else {
      return("duxx")
    }
  }
  
  df_combined <- df %>%
    rowwise() %>%
    mutate(across(all_of(regions), as.character))
  for (comb in combinations) {
    col1 <- comb[1]
    col2 <- comb[2]
    new_col <- paste0(col1,"_", col2)
    df_combined <- df_combined %>%
      mutate(!!new_col := classify_combination(!!sym(col1), !!sym(col2)))
  }
  # Return the processed dataframe
  df_combined
}


processed_dfs <- mapply(process_celltype, file_list, cell_types, SIMPLIFY = FALSE)
names(processed_dfs) <- cell_types

wb <- createWorkbook()
for (ct in cell_types) {
  addWorksheet(wb, ct)
  writeData(wb, ct, processed_dfs[[ct]])
}
saveWorkbook(wb, "pairwise_classification_dataframes.xlsx", overwrite = TRUE)
