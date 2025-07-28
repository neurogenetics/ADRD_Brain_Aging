#Code conceptualized and created by Monica E. Mesecar and Dom J. Acri, with assistance from Perplexity AI
#Specifically on plot aesthetics and looping through multiple files for efficiency.

# Libraries
library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggplot2)

# List of input files and corresponding cell type labels
#Shape: Gene Name on Y axis (rows), Region on X (columns), fill with Up, Dowm, NE (not expressed)
file_list <- c("Inhib_ByRegion_aDEG_fixed.csv", "ExN_ByRegion_aDEG_fixed.csv","OD_ByRegion_aDEG_fixed.csv", "OPC_ByRegion_aDEG_fixed.csv","Mural_ByRegion_aDEG_fixed.csv", "Astro_ByRegion_aDEG_fixed.csv","Microglia_ByRegion_aDEG_fixed.csv", "Epend_ByRegion_aDEG_fixed.csv", "Endo_ByRegion_aDEG_fixed.csv") # Add more as needed
cell_types <- c("InN", "ExN", "Oligo", "OPC","Mural", "Astro", "Micro", "Epend", "Endo") # Should match file_list order

# Function to process each file, loop through
process_celltype <- function(file, celltype) { #Input is input file name + label)
  df <- read.csv(file)
  regions <- colnames(df)[2:5] #Region list
  combinations <- combn(regions, 2, simplify = FALSE) #Pairwise combinations
  
  #Shorthand for combinations of expression patterns 
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
  
  # Prepare df_combined
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
  
  # Concordance
  #Count number of combinations and fill in HM 
  pairwise_counts <- data.frame(matrix(0, nrow = length(regions), ncol = length(regions)))
  rownames(pairwise_counts) <- colnames(pairwise_counts) <- regions
  for (comb in combn(regions, 2, simplify = FALSE)) {
    col1 <- comb[1]
    col2 <- comb[2]
    uu_count <- sum(df_combined[[paste0(col1,"_", col2)]] == "uu", na.rm = TRUE)
    dd_count <- sum(df_combined[[paste0(col1, "_",col2)]] == "dd", na.rm = TRUE)
    pairwise_counts[col1, col2] <- uu_count
    pairwise_counts[col2, col1] <- dd_count
  }
  pairwise_counts[is.na(pairwise_counts)] <- 0
  pairwise_counts <- as.matrix(pairwise_counts)
  final_counts <- matrix(0, nrow = nrow(pairwise_counts), ncol = ncol(pairwise_counts), 
                         dimnames = dimnames(pairwise_counts))
  final_counts[upper.tri(final_counts)] <- pairwise_counts[upper.tri(pairwise_counts)]
  final_counts[lower.tri(final_counts)] <- -pairwise_counts[lower.tri(pairwise_counts)]
  diag(final_counts) <- NA
  heatmap_data <- melt(final_counts, na.rm = FALSE)
  heatmap_data$border_color <- ifelse(
    as.integer(heatmap_data$Var1) < as.integer(heatmap_data$Var2), "red",
    ifelse(as.integer(heatmap_data$Var1) > as.integer(heatmap_data$Var2), "blue", "grey")
  )
  heatmap_data$fill_color <- ifelse(
    as.integer(heatmap_data$Var1) == as.integer(heatmap_data$Var2), "grey", "white"
  )
  p1 <- ggplot(data = heatmap_data, aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = fill_color, color = border_color), size = 1.2, show.legend = FALSE, na.rm = FALSE) +
    geom_text(aes(label = ifelse(is.na(value), "", abs(value))), vjust = 1, hjust = 1) +
    scale_fill_identity() +
    scale_color_identity() +
    theme_classic() +
    labs(title = paste(celltype, "Pairwise Concordant aDEGs"), 
         x = "Region 2", 
         y = "Region 1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  # Discordance
  pairwise_counts <- data.frame(matrix(0, nrow = length(regions), ncol = length(regions)))
  rownames(pairwise_counts) <- colnames(pairwise_counts) <- regions
  for (comb in combn(regions, 2, simplify = FALSE)) {
    col1 <- comb[1]
    col2 <- comb[2]
    du_count <- sum(df_combined[[paste0(col1,"_", col2)]] == "du", na.rm = TRUE)
    ud_count <- sum(df_combined[[paste0(col1, "_",col2)]] == "ud", na.rm = TRUE)
    pairwise_counts[col1, col2] <- du_count
    pairwise_counts[col2, col1] <- ud_count
  }
  pairwise_counts[is.na(pairwise_counts)] <- 0
  pairwise_counts <- as.matrix(pairwise_counts)
  final_counts <- matrix(0, nrow = nrow(pairwise_counts), ncol = ncol(pairwise_counts), 
                         dimnames = dimnames(pairwise_counts))
  final_counts[upper.tri(final_counts)] <- pairwise_counts[upper.tri(pairwise_counts)]
  final_counts[lower.tri(final_counts)] <- -pairwise_counts[lower.tri(pairwise_counts)]
  diag(final_counts) <- NA
  heatmap_data <- melt(final_counts, na.rm = FALSE)
  heatmap_data$border_color <- ifelse(
    as.integer(heatmap_data$Var1) < as.integer(heatmap_data$Var2), "magenta4",
    ifelse(as.integer(heatmap_data$Var1) > as.integer(heatmap_data$Var2), "magenta4", "grey")
  )
  heatmap_data$fill_color <- ifelse(
    as.integer(heatmap_data$Var1) == as.integer(heatmap_data$Var2), "grey", "white"
  )
  p2 <- ggplot(data = heatmap_data, aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = fill_color, color = border_color), size = 1.2, show.legend = FALSE, na.rm = FALSE) +
    geom_text(aes(label = ifelse(is.na(value), "", abs(value))), vjust = 1, hjust = 1) +
    scale_fill_identity() +
    scale_color_identity() +
    theme_classic() +
    labs(title = paste(celltype, "Pairwise Discordant aDEGs"), 
         x = "Down aDEGs", 
         y = "Up aDEGs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  # Return plots 
  list(concordant_plot = p1, discordant_plot = p2)
}

# Loop through all files/cell types
plot_results <- mapply(process_celltype, file_list, cell_types, SIMPLIFY = FALSE)

# Assuming plot_results is a list of lists
all_plots <- unlist(plot_results, recursive = FALSE)  # Flattens to a single list of plots

# Save to PDF
ggsave(
  filename = "all_celltype_heatmaps.pdf",
  plot = marrangeGrob(all_plots, nrow = 1, ncol = 1),
  width = 8, height = 6
)

