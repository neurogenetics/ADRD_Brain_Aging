#libraries
library(Seurat)
library(tidyverse)

#Load in senescence pathways as named list of genes
sen_paths <- readRDS(file="/data/ADRD/brain_aging/exploration/revision2_senpaths/SenPaths_completed_final.rds")

#identical(SenNet_2024_CNS_RNA, SenNet_2024_CNS_sc)--> TRUE
#Genes not in RNA: "H2AX"   "GLB1"   "CDKN2B"

#load object
obj = readRDS("/data/ADRD/brain_aging/exploration/aging.pegasus.leiden_085.subclustered.rds")

#OG HVFs
og_hvf <- rownames(obj@assays$RNA@meta.features)[obj@assays$RNA@meta.features$highly.variable_features]
#og_hvf_list <- as.list(og_hvf)

#get variable features
obj <- FindVariableFeatures(obj)
varfeat <- VariableFeatures(obj)

#Get All Expressed Gene list
exp_gene_df=data.frame(rownames(obj@assays$RNA))
exp_genes=as.character(exp_gene_df[[1]])

################upset plot1: c(varfeat, sen_paths)####################################
##Create Binary Incidence DF
#0.5 Unique SenPath Genes
senpath_genes<-unique(unlist(sen_paths))

#1.5 Euler Plot of VarFeat x Sen_Paths
library(eulerr)

# Prepare a named vector listing the sizes for each region or give a list of sets.
# Option A: give sets directly (eulerr will compute regions)
fit <- euler(list("SenPath Genes" = unique(senpath_genes), "Phase1 Variable Features (Seurat)" = unique(varfeat)))

# Inspect fit
print(fit)

# Basic plot
plot(fit,
     fills = list(fill = c("#66c2a5","#fc8d62"), alpha = 0.6),
     labels = list(font = 4),        # font style for labels
     edges = TRUE,
     quantities = TRUE)              # shows counts in the regions

# Prepare a named vector listing the sizes for each region or give a list of sets.
# Option A: give sets directly (eulerr will compute regions)
fit2 <- euler(list("SenPath Genes" = unique(senpath_genes), "Phase1 Variable Features (OG)" = unique(og_hvf)))

# Inspect fit
print(fit2)

# Basic plot
plot(fit2,
     fills = list(fill = c("#66c2a5","gold"), alpha = 0.6),
     labels = list(font = 4),        # font style for labels
     edges = TRUE,
     quantities = TRUE)              # shows counts in the regions


##############################
# 2. Convert the list of gene vectors into a binary presence/absence data frame
# This identifies every unique gene across all lists and checks its presence in each
all_genes <- unique(unlist(sen_paths))
upset_data <- data.frame(gene = all_genes)

for (path_name in names(sen_paths)) {
  upset_data[[path_name]] <- as.integer(all_genes %in% sen_paths[[path_name]])
}

# 3. Define the sets to plot
pathway_names <- names(sen_paths)

# 4. Create the UpSet Plot
# We use 'n_intersections' to limit to the most frequent overlaps for clarity
p_upset <- upset(
  upset_data,
  pathway_names,
  name = "Senescence Pathway Overlaps",
  width_ratio = 0.15,  # Adjusts the width of the set size bars on the left
  stripes = upset_stripes(
    geom = geom_segment(size = 5),
    colors = c('grey95', 'white')
  ),
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = 'bars_color')
    ) + scale_fill_manual(values = c('bars_color' = '#3182bd'), guide = 'none')
  ),
  matrix = intersection_matrix(
    geom = geom_point(size = 3)
  ),
  set_sizes = (
    upset_set_size() + 
      theme(axis.text.x = element_text(angle = 90)) +
      geom_text(aes(label = ..count..), stat = 'count', hjust = -0.1, size = 3)
  )
) +
  labs(
    title = "Overlap Analysis of Curated Senescence Gene Sets",
    caption = "PMID: 37933854 (CellAge v5), Dehkordi 2021, and SenNet 2024 consensus"
  ) +
  theme_minimal(base_family = "Arial")

p_upset

#dev.off()

############################################
library(ComplexUpset)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Your existing inputs (you already have these) ---
# sen_paths : named list of 9 pathways
# varfeat   : character vector of varfeat genes

# --- 1) Build upset_data with gene x pathway binary columns (your code) ---
all_genes <- unique(unlist(sen_paths))
upset_data <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

for (path_name in names(sen_paths)) {
  upset_data[[path_name]] <- as.integer(all_genes %in% sen_paths[[path_name]])
}

# --- 1b) Add column to indicate membership in varfeat (1/0 or TRUE/FALSE) ---
# (You asked to add this column to the input dataframe.)
upset_data$in_varfeat <- as.integer(upset_data$gene %in% varfeat)
# If you prefer logical:
#upset_data$in_varfeat <- upset_data$gene %in% varfeat

#Save CSV
#output_file_path <- "SenPath_VarFeat_BinaryIncidence.csv"
#write.csv(upset_data, file= output_file_path, row.names = FALSE)

# Define pathway names (sets)
pathway_names <- names(sen_paths)

##########WIP############################

# --- 2) Use upset_data() to compute intersections and prepare stacked counts ---
# Note: upset_data() (ComplexUpset) produces an object with intersection info
up <- upset_data(upset_data, intersect = pathway_names, intersections = "observed")

ints_df <- up$with_sizes

# Try to find a list-column containing the intersection members (common names)
members_col_name <- names(ints_df)[sapply(ints_df, function(col) inherits(col, "list") && all(sapply(col, function(x) is.character(x) || length(x) == 0))))][1]

if (!is.na(members_col_name) && !is.null(members_col_name)) {
  # If we have members as a list-column, compute counts directly
  ints_df <- ints_df %>%
    rowwise() %>%
    mutate(
      intersection_members = .data[[members_col_name]],
      n_in_intersection = length(intersection_members),
      n_in_varfeat = sum(intersection_members %in% varfeat),
      n_not_in_varfeat = n_in_intersection - n_in_varfeat
    ) %>%
    ungroup()
} else {
  # Fall-back: reconstruct intersection members from the binary columns in upset_data
  # Identify the column in ints_df that labels the intersection; common names include 'intersection' etc.
  inter_id_col <- names(ints_df)[grepl("intersect|intersection", names(ints_df), ignore.case = TRUE)][1]
  if (is.null(inter_id_col) || is.na(inter_id_col)) inter_id_col <- names(ints_df)[1]
  
  parse_sets_from_label <- function(lbl) {
    s <- gsub("[\\{\\}\\s]", "", as.character(lbl))
    parts <- unlist(strsplit(s, split = "[&:,;|\\+]+"))
    parts[parts != ""]
  }
  
  ints_df <- ints_df %>%
    mutate(.intersection_label = as.character(.data[[inter_id_col]])) %>%
    rowwise() %>%
    mutate(
      active_sets = list(parse_sets_from_label(.intersection_label)),
      intersection_members = list({
        active <- active_sets
        if (length(active) == 0) {
          character(0)
        } else {
          # Build logical vector: gene must be in all active sets and NOT in inactive sets
          keep <- rep(TRUE, nrow(upset_data))
          for (s in pathway_names) {
            if (s %in% active) keep <- keep & (upset_data[[s]] == 1)
            else keep <- keep & (upset_data[[s]] == 0)
          }
          upset_data$gene[keep]
        }
      }),
      n_in_intersection = length(intersection_members),
      n_in_varfeat = sum(intersection_members %in% varfeat),
      n_not_in_varfeat = n_in_intersection - n_in_varfeat
    ) %>%
    ungroup()
}

# --- 3) Build a stacked dataframe for plotting (one row per intersection × status) ---
# Determine which column to use as the intersection ID label (so mapping aligns)
id_col_candidates <- c("intersection", "intersect", "intersections", "intersect_id", "intersection_label", ".intersection_label")
id_col <- id_col_candidates[id_col_candidates %in% names(ints_df)][1]
if (is.null(id_col) || is.na(id_col)) id_col <- names(ints_df)[1]

# Create a column exactly named 'intersection' because ComplexUpset commonly expects that label
ints_df$intersection <- as.character(ints_df[[id_col]])

stack_df <- ints_df %>%
  select(intersection, n_in_varfeat, n_not_in_varfeat) %>%
  pivot_longer(cols = c(n_in_varfeat, n_not_in_varfeat),
               names_to = "status", values_to = "count") %>%
  mutate(
    status = factor(status, levels = c("n_in_varfeat", "n_not_in_varfeat"),
                    labels = c("in_varfeat", "not_in_varfeat"))
  )

# --- 4) Compose the UpSet plot with custom stacked top bars ---
# Colors: lower segment (in_varfeat) = yellow, remainder = black
my_colors <- c("in_varfeat" = "yellow", "not_in_varfeat" = "black")

# Use a custom top annotation: geom_col with data=stack_df (position='stack')
intersection_annotation <- list(
  'Intersection (stacked in_varfeat)' = list(
    aes = aes(x = intersection, y = count, fill = status),
    geom = list(
      geom_col(data = stack_df, stat = "identity", position = "stack", width = 0.9)
    )
  )
)

# Build the plot (replace your previous upset() call)
p_upset <- upset(
  upset_data,
  pathway_names,
  name = "Senescence Pathway Overlaps",
  width_ratio = 0.15,
  stripes = upset_stripes(
    geom = geom_segment(size = 5),
    colors = c('grey95', 'white')
  ),
  base_annotations = intersection_annotation,
  matrix = intersection_matrix(
    geom = geom_point(size = 3)
  ),
  set_sizes = (
    upset_set_size() +
      theme(axis.text.x = element_text(angle = 90)) +
      geom_text(aes(label = ..count..), stat = 'count', hjust = -0.1, size = 3)
  )
) +
  scale_fill_manual(name = "", values = my_colors, guide = "legend") +
  labs(
    title = "Overlap Analysis of Curated Senescence Gene Sets",
    caption = "PMID: 37933854 (CellAge v5), Dehkordi 2021, and SenNet 2024 consensus"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "right")

# Print it
print(p_upset)