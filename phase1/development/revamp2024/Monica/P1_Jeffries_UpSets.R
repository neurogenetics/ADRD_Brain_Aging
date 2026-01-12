# Set the new working directory to the specified path
setwd("/data/ADRD/brain_aging/exploration")
library(dplyr)
library(tidyverse)
library(readxl)
library(UpSetR)
library(openxlsx)

#Allow ample time to download file (27 GB for 51 GB Seurat and others)
#options(timeout = 3600)  # 1 hour

#Specify Path
#url <- "http://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/pfc.clean.cpm.tar.gz"

#Download
#download.file(
  #url,
  #destfile = "./pfc.clean.cpm.tar.gz",
  #mode = "wb",
  #method = "libcurl"
#)

#Read in Seurat Obj. 
#pfc_clean_seurat<-readRDS("Jefferies_2025_data/pfc.clean.rds")

#Read in Excel Files
sample_info <- read_excel("Jefferies_2025_data/sample_information.xlsx")
Jeffries_DEGs_EA<-read_excel("Jefferies_2025_data/S6_DEGs_elderly_vs_adult.xlsx")
#################################################################################
#Clean Jeffries Excel Sheets to be Workable 
# 1. Assign the third row as the column names
# The index [3, ] selects the third row and all columns
names(Jeffries_DEGs_EA) <- Jeffries_DEGs_EA[3, ]

# 2. Remove the first three rows from the dataframe
# The index -c(1:3) removes rows 1, 2, and 3
Jeffries_DEGs_EA <- Jeffries_DEGs_EA[-c(1:3), ]

# (Optional) Re-index the row names
# This makes the row names sequential again starting from 1
rownames(Jeffries_DEGs_EA) <- 1:nrow(Jeffries_DEGs_EA)
################################################################################
#Filter down their DEG df to just broad cell types that are reflected in our dataset
filtered_Jeffries <- filter(Jeffries_DEGs_EA,`cell type` %in% c("oli", "opc", "ext", "inb", "micro", "ast", "endo"))

#Read in our aDEG file, filter for dataset and cell types to be analyzed
age <- read_csv("/data/ADRD/brain_aging/phase1/results/aging.glmmtmb_age_diffs_fdr.csv")
age_filtered <- age %>%
  filter(type == "region_broad_celltype") %>%
  filter(!tissue %in% c(
    "Entorhinal cortex SPN",
    "Middle temporal gyrus SPN",
    "Subventricular zone SPN"
  ))

#Adjust first column name 
names(age_filtered)[1] <- "gene" # Change column 1 name from feature to gene

#Add columns for cell type and region
age_filtered <- age_filtered %>%
  mutate(
    # celltype is the final word
    celltype = word(tissue, -1),
    
    # region is everything except the final word
    region   = str_trim(str_remove(tissue, paste0("\\s+", celltype, "$")))
  )

#Rename Jeffries CTs to match our naming scheme
filtered_Jeffries <- filtered_Jeffries%>%
  mutate(`cell type` = recode(`cell type`,
                            "oli" = "Oligodendrocyte",
                            "opc" = "OPC",
                            "ext" = "ExN",
                            "inb" = "InN",
                            "micro" = "Microglia",
                            "ast" = "Astrocyte",
                            "endo" = "Endothelial"))

#Select overlapping broad types from Jeffries and ours 
keep_types <- c("Oligodendrocyte", "OPC", "ExN", "InN",
                    "Microglia", "Astrocyte", "Endothelial")

#Filter our dataset down to just those broad types
age_filtered_Jtypes <- age_filtered %>%
  filter(celltype %in% keep_types)

#Examine to make sure they match
unique(filtered_Jeffries$`cell type`)
unique(age_filtered_Jtypes$celltype)

#Rename Jeffries columns for clarity/readability 
filtered_Jeffries <- filtered_Jeffries %>%
  rename(celltype = `cell type`)

filtered_Jeffries <- filtered_Jeffries %>%
  rename(gene = `gene name`)
###########################################################
#Upset Dataframe Prep

## Assume:
## df1: gene, celltype
## df2: gene, tissue   (e.g. "RegionA Excitatory")

# Helper functions for study 2: split "RegionA Excitatory"
get_ct     <- function(x) sub(".*\\s+", "", x)        # last token
get_region <- function(x) sub("\\s+[^ ]+$", "", x)    # everything but last

# Study 1: Jeffries, parse to specify no region and study flag
Jeffries_parsed <- filtered_Jeffries %>%
  rename(ct = celltype) %>%
  mutate(
    region = NA_character_,           # no region info in study1
    set    = paste0(ct, "_Jeffries")    # e.g. "Excitatory_Study1"
  )

# Study 2: from 'tissue' = "Region CellType", parse to get region_celltype columns (called set)
age_parsed <- age_filtered_Jtypes %>%
  rename(tissue_raw = tissue) %>%
  mutate(
    ct     = get_ct(tissue_raw),                     # "Excitatory"
    region = get_region(tissue_raw),                 # "RegionA"
    set    = paste(region, ct, sep = "_")            # "RegionA_Excitatory"
  )

#List of overlapping cell types 
celltypes_all <- intersect(
  unique(Jeffries_parsed$ct),
  unique(age_parsed$ct)
)
#Examine 
celltypes_all

################### UpSet Create Binary Matrix Function #####################
make_upset_for_celltype <- function(ct_of_interest, df1_parsed, df2_parsed) {
  
  # Keep only this cell type
  df1_ct <- df1_parsed %>%
    filter(ct == ct_of_interest) %>%
    select(gene, set)
  
  df2_ct <- df2_parsed %>%
    filter(ct == ct_of_interest) %>%
    select(gene, set)
  
  # Combine rows and create presence indicators
  df_long <- bind_rows(df1_ct, df2_ct) %>%
    distinct() %>%
    mutate(present = 1L)
  
  df_wide <- df_long %>%
    pivot_wider(
      id_cols      = gene,
      names_from   = set,
      values_from  = present,
      values_fill  = 0L
    )
  
  mat <- as.data.frame(df_wide[, -1, drop = FALSE])
  rownames(mat) <- df_wide$gene
  
  mat
}
###################################################
#Save each CT matrix into list, create matrix
upset_list <- lapply(
  celltypes_all,
  make_upset_for_celltype,
  df1_parsed = Jeffries_parsed,
  df2_parsed = age_parsed
)

# Name each list element
names(upset_list) <- celltypes_all

#################################################
#Examine for correctness
View(upset_list$Oligodendrocyte)

#################################################
# Define the filename for the PDF
pdf_filename <- "Endothelial2_UpSet_Plot.pdf"

# Open the PDF graphics device
# The 'width' and 'height' arguments (in inches) are optional but recommended
# to ensure the plot is correctly proportioned.
pdf(file = pdf_filename, width = 2.36, height = 1.65)

# Generate the UpSet plot (The output is now directed to the PDF file)
upset(
  upset_list[[ "Endothelial" ]], 
  sets = colnames(upset_list[[ "Endothelial" ]])
)

# Close the graphics device
# This finalizes the PDF file and saves it to your working directory
dev.off()

# Print a message to confirm the file has been saved
print(paste("UpSet plot saved to:", pdf_filename))
######################################################################
#Save UpSet Matricies into Excel workbook where each CT is a sheet

# Create a workbook
#wb <- createWorkbook()

#Grab all CTs 
#for (ct in names(upset_list)) {
  
  #df <- upset_list[[ct]]
  
  # Move rownames into a real column
  #df_out <- cbind(gene = rownames(df), df) #Keeps gene names
  
  #addWorksheet(wb, sheetName = ct) #CTs as sheets 
  
  #writeData(wb, sheet = ct, x = df_out, rowNames = FALSE)
#}

#Save workbook to working directory
#saveWorkbook(wb, file = "celltype_upset_sets.xlsx", overwrite = TRUE)

############## Curious about aDEG distribution####################

#We have 11K+ aDEGs (all 4 regions, total not unique)--4.6K EC, 3.5K MTG, Putamen 1.2K, SVZ 1.9K
# Astro: 1.9K, Endo: 49, 1.8K ExN, 3K InN, 750 Micro, 2.6K Oligo, 1.1K OPC

#Jeffries has just under 3K....
#200 Astro, 200 Endo, 1300 ExN, 600 InN, 224 Micro, 150 Oligo, 250 OPC

print("Jeffries Values")
table(filtered_Jeffries$celltype)

print("Our Values")
table(age_filtered_Jtypes$celltype)
table(age_filtered_Jtypes$tissue)
table(age_parsed$region)
table(age_parsed$ct)
