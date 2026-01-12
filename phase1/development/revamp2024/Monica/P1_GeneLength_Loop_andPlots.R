# Gene Length Interrogation - Tissue-specific Analysis (Wilcoxon Test & Faceted Plot)

# Load Packages (assuming they are already installed as per original script)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(biomaRt)
library(readr)
library(stringr)
library(ggplot2)
library(ggridges)

# --- 1. Gene Length Preparation (No changes from previous script) ---

# Load Ref GTF and get gene lengths
gtf_path <- "/gpfs/gsfs12/users/mesecarme/Phase1_Revision/Gene_Length/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

gtf <- import(gtf_path)
txdb <- makeTxDbFromGRanges(gtf)

gene_ranges <- genes(txdb)

if ("gene_id" %in% colnames(mcols(gene_ranges))) {
  gene_ids <- as.character(mcols(gene_ranges)$gene_id)
} else {
  gene_ids <- names(gene_ranges)
}

gene_span_df <- tibble(
  gene_id       = gene_ids,
  gene_length = as.numeric(width(gene_ranges))
)

# Get mapping: Ensembl gene ID -> HGNC symbol
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids_to_map <- sub("\\.\\d+$", "", gene_span_df$gene_id)

mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters     = "ensembl_gene_id",
  values      = unique(gene_ids_to_map),
  mart        = mart
)

# Add symbols; create final canon_lengths dataframe
canon_lengths <- gene_span_df %>%
  mutate(gene_id = sub("\\.\\d+$", "", gene_id)) %>%       # remove version suffix
  left_join(mapping, by = c("gene_id" = "ensembl_gene_id")) %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%      # drop unmapped
  dplyr::select(hgnc_symbol, gene_length) %>%
  rename(gene_symbol = hgnc_symbol) %>%
  mutate(gene_symbol = str_trim(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE) # ensure 1 row per symbol

# --- 2. Initial Data Loading and Filtering (No changes from previous script) ---

# Read in aDEG file
age <- read_csv("/data/ADRD/brain_aging/phase1/results/aging.glmmtmb_age_diffs_fdr.csv")
age_filtered <- age %>%
  filter(type == "region_broad_celltype") %>%
  filter(!tissue %in% c(
    "Entorhinal cortex SPN",
    "Middle temporal gyrus SPN",
    "Subventricular zone SPN"
  ))
names(age_filtered)[1] <- "gene_symbol" # Change column 1 name from feature to gene_symbol

# Read in all gene file
all_gene <- read_csv("/data/ADRD/brain_aging/phase1/results/aging.diffxpy_age_diffs.csv")
all_filtered <- all_gene %>%
  filter(type == "region_broad_celltype") %>%
  filter(!tissue %in% c(
    "Entorhinal cortex SPN",
    "Middle temporal gyrus SPN",
    "Subventricular zone SPN"
  ))
names(all_filtered)[1] <- "gene_symbol" # Change column 1 name from feature to gene_symbol

# --- 3. Identify Unique Tissues (No changes from previous script) ---

# Combine unique tissues from both filtered dataframes
unique_tissues <- unique(c(age_filtered$tissue, all_filtered$tissue))
unique_tissues <- unique_tissues[!is.na(unique_tissues)] # Remove NA if present

# --- 4. Tissue Analysis Function (UPDATED: Using wilcox.test) ---

#' Perform gene length analysis for a specific tissue.
#'
#' @param current_tissue The tissue value to filter by.
#' @param age_data The main DEG dataframe (age_filtered).
#' @param all_data The main ALL dataframe (all_filtered).
#' @param lengths_data The canonical gene length dataframe.
#' @return A list containing the df_plot and the wilcox.test result.
analyze_tissue_gene_length <- function(current_tissue, age_data, all_data, lengths_data) {
  
  # 1. Filter DEG and ALL data for the current tissue
  tissue_age_filtered <- age_data %>%
    filter(tissue == current_tissue)
  
  tissue_all_filtered <- all_data %>%
    filter(tissue == current_tissue)
  
  # 2. Map gene lengths for DEGs
  df_deg_unique <- tissue_age_filtered %>%
    distinct(gene_symbol)
  df_deg_lengths <- df_deg_unique %>%
    left_join(lengths_data, by = "gene_symbol") %>%
    filter(!is.na(gene_length)) # Exclude genes without length for analysis
  
  # 3. Map gene lengths for ALL genes and identify NON-DEGs
  df_all_unique <- tissue_all_filtered %>%
    distinct(gene_symbol)
  df_all_lengths <- df_all_unique %>%
    left_join(lengths_data, by = "gene_symbol")
  
  df_non_deg_lengths <- df_all_lengths %>%
    filter(!(gene_symbol %in% df_deg_lengths$gene_symbol)) %>% # Exclude DEGs
    filter(!is.na(gene_length)) # Exclude genes without length
  
  # 4. Combine DEGs and Non-DEGs for Plotting/Testing
  df_deg_lengths$group <- "DEG"
  df_non_deg_lengths$group <- "Non-DEG"
  
  df_plot <- rbind(
    df_deg_lengths[, c("gene_symbol", "gene_length", "group")],
    df_non_deg_lengths[, c("gene_symbol", "gene_length", "group")]
  ) %>%
    # Add tissue column
    mutate(tissue = current_tissue)
  
  # 5. Wilcoxon Test (using log10 transformed length)
  # Check if both groups exist and have enough data
  if (sum(df_plot$group == "DEG") >= 3 & sum(df_plot$group == "Non-DEG") >= 3) {
    wilcox_result <- tryCatch({
      wilcox.test(log10(gene_length) ~ group, data = df_plot)
    }, error = function(e) {
      list(
        statistic = c(W = NA),
        p.value = NA,
        error_msg = paste("Wilcox test failed:", as.character(e))
      )
    })
  } else {
    wilcox_result <- list(
      statistic = c(W = NA),
      p.value = NA,
      error_msg = "Insufficient data points (less than 3 in a group) for Wilcox test."
    )
  }
  
  return(list(df_plot = df_plot, wilcox_test = wilcox_result))
}

# --- 5. Run Analysis Loop and Store Results ---

# Initialize list to store results for each tissue
results_list <- list()
wilcox_summary_list <- list()
tissue_name=unique_tissues[1]

cat("Starting analysis for", length(unique_tissues), "unique tissues...\n")

for (tissue_name in unique_tissues) {
  cat("Analyzing tissue:", tissue_name, "...\n")
  tissue_results <- analyze_tissue_gene_length(
    current_tissue = tissue_name,
    age_data = age_filtered,
    all_data = all_filtered,
    lengths_data = canon_lengths
  )
  
  # Store the df_plot for the current tissue
  results_list[[tissue_name]] <- tissue_results$df_plot
  
  # Summarize the Wilcoxon test results
  wc_test <- tissue_results$wilcox_test
  
  # Calculate median log10 lengths for the summary table
  median_lengths <- tissue_results$df_plot %>%
    group_by(group) %>%
    summarise(median_log10_length = median(log10(gene_length)))
  
  if (is.na(wc_test$p.value)) {
    summary_row <- data.frame(
      tissue = tissue_name,
      N_DEG = sum(tissue_results$df_plot$group == "DEG"),
      N_NonDEG = sum(tissue_results$df_plot$group == "Non-DEG"),
      median_log10_length_DEG = ifelse("DEG" %in% median_lengths$group, median_lengths$median_log10_length[median_lengths$group == "DEG"], NA),
      median_log10_length_NonDEG = ifelse("Non-DEG" %in% median_lengths$group, median_lengths$median_log10_length[median_lengths$group == "Non-DEG"], NA),
      W_statistic = NA,
      p_value = NA,
      Notes = wc_test$error_msg
    )
  } else {
    summary_row <- data.frame(
      tissue = tissue_name,
      N_DEG = sum(tissue_results$df_plot$group == "DEG"),
      N_NonDEG = sum(tissue_results$df_plot$group == "Non-DEG"),
      median_log10_length_DEG = median_lengths$median_log10_length[median_lengths$group == "DEG"],
      median_log10_length_NonDEG = median_lengths$median_log10_length[median_lengths$group == "Non-DEG"],
      W_statistic = wc_test$statistic,
      p_value = wc_test$p.value,
      Notes = ""
    )
  }
  wilcox_summary_list[[tissue_name]] <- summary_row
}

# --- 6. Final Output ---

## 🧬 Wilcoxon Rank-Sum Test Summary

# Final Wilcoxon summary dataframe
df_wilcox_summary <- do.call(rbind, wilcox_summary_list)
print(df_wilcox_summary)

df_wilcox_summary$FDR_BH<-p.adjust(df_wilcox_summary$p_value, method = "BH")
write.csv(df_wilcox_summary, "wilcoxon_gene_length_summary.csv", row.names = FALSE)

df_combined_plot <- do.call(rbind, results_list)
rownames(df_combined_plot) <- NULL # Clean up row names

# Faceted ggplot
faceted_length_plot <- ggplot(df_combined_plot, aes(x = log10(gene_length), color = group)) +
  geom_density(size = 1.2) +
  facet_wrap(~ tissue, scales = "free_y", ncol = 4) + # Facet by tissue
  labs(
    x = "log10(Gene Length)",
    y = "Density",
    color = "Gene Group",
    title = "Length Distribution of DEGs vs Non-DEGs Faceted by Tissue"
  ) +
  scale_color_manual(values = c("DEG" = "red", "Non-DEG" = "blue")) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold") # Make facet labels bold
  )

print(faceted_length_plot)
# Save the Faceted Density Plot as PDF
ggsave(
  filename = "CT_faceted_gene_length_density.pdf", 
  plot = faceted_length_plot, 
  width = 12, 
  height = 8
)


df_plot_parsed <- df_combined_plot %>%
  # Separate the tissue string by the last space encountered
  mutate(
    # Capture everything before the last space (Region) and everything after (Cell Type)
    region = str_replace(tissue, " [^ ]+$", ""), # Replace last word with nothing
    cell_type = str_extract(tissue, "[^ ]+$")    # Extract last word
  ) %>%
  # Reorder cell types by median gene length for a nicer ridge plot display
  # Convert cell_type to a factor based on median log10(gene_length)
  mutate(
    log10_length = log10(gene_length)
  ) %>%
  group_by(cell_type) %>%
  mutate(
    median_length = median(log10_length, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    cell_type = factor(cell_type, levels = unique(cell_type[order(median_length, decreasing = FALSE)]))
  )
# Generate the Ridge Plot
ridge_plot <- ggplot(df_plot_parsed, aes(x = log10_length, y = cell_type, fill = group)) +
  geom_density_ridges(
    alpha = 0.7,        # Set transparency
    scale = 1.0,        # Adjust overlap
    rel_min_height = 0.01 # Filter out tiny noise densities
  ) +
  facet_wrap(~ region, scales = "free_y", ncol = 2) + # Facet by Region
  scale_fill_manual(values = c("DEG" = "#E64B35", "Non-DEG" = "#4DBBD5")) +
  labs(
    x = expression(log[10] * "(Gene Length)"),
    y = "Cell Type",
    fill = "Gene Group",
    title = "Gene Length Distribution (DEGs vs Non-DEGs) by Cell Type and Region"
  ) +
  theme_ridges(grid = TRUE) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10) # Improve facet labels
  )

print(ridge_plot)
# Save the Ridge Plot as PDF
ggsave(
  filename = "ridge_plot_gene_length_celltype.pdf", 
  plot = ridge_plot, 
  width = 15, 
  height = 11
)

# Calculate the median of log10 gene lengths for the histogram line
median_log10_length <- median(log10(gene_span_df$gene_length), na.rm = TRUE)

# Gene length (Whole Genome) - Display on screen
hist(log10(gene_span_df$gene_length), 
     main = "Whole-Genome Gene Length Distribution",
     xlab = expression(log[10] * "(Gene Length)"))
abline(v = median_log10_length, col = "red", lty = 2, lwd = 2) # Red vertical dashed line at median
# Add legend for the median line
legend("topright", 
       legend = c(paste0("Median: ", round(median_log10_length, 2))), 
       col = "red", 
       lty = 2, 
       lwd = 2, 
       bty = "n") # No box around legend

# Save the Histogram as PDF
pdf("whole_genome_gene_length_histogram.pdf")
hist(log10(gene_span_df$gene_length), 
     main = "Whole-Genome Gene Length Distribution",
     xlab = expression(log[10] * "(Gene Length)"))
abline(v = median_log10_length, col = "red", lty = 2, lwd = 2) # Red vertical dashed line at median
# Add legend for the median line
legend("topright", 
       legend = c(paste0("Median: ", round(median_log10_length, 2))), 
       col = "red", 
       lty = 2, 
       lwd = 2, 
       bty = "n") # No box around legend
dev.off()