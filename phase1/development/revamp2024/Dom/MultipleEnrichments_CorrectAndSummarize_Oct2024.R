################
#Date: 2024/10/21
#Authors: chatGPT, Gemini, Perplexity
#Code Curators: Monica, Amal, Dom
#Goals:
#    - Implement background adjustment on gprofiler2 outputs (in R)
#    - Allow multi_query (this was not possible in base gprofiler2 w/ multiple background of diff length)
#         - Did with lapply, needs to be BH corrected (?)
#         - Compared background correction vs none (saved as pdf w/ counts and ratio plots)
#    - summarized GO parent relationships with rrvgo 
#         - only did for ALL source::GO as proof of concentp
#         - need to do per celltype, cluster, region_celltype
################
#libraries
library(gprofiler2)
library(tidyverse)
library(ggpubr)
library(rrvgo) #bioconductor!!
library(plotly)


#################
#set background
dat=read.csv("/data/ADRD/brain_aging/phase1/results/aging.diffxpy_age_diffs.csv") #all observed genes
summary(factor(dat$tissue))
res=read.csv("/data/ADRD/brain_aging/phase1/results/aging.glmmtmb_age_diffs_fdr.csv") #final aDEGs
summary(factor(res$tissue))
#####Set background
# Assuming you have the following data structures:
# - res: A data frame with columns "tissue" and "feature"
# - dat: A data frame with columns "tissue" and "gene"
# Create a list of queries based on unique tissues in res$tissue
queries <- split(res$feature, res$tissue)
# Create a list of custom backgrounds based on unique tissues in dat$tissue
custom_backgrounds <- lapply(split(dat$gene, dat$tissue), identity)
# Perform multi-query without custom backgrounds
results_nobckgd <- gost(queries, organism = "hsapiens")
################
## need to rerun withevcodes = T
###################
# Access results for each tissue without custom backgrounds
source_lengths_nobkgd=summary(factor(results_nobckgd$result$query))
##cannot run multi-background with different custom_bg
# Run gost for each query with its respective background
results <- lapply(1:length(queries), function(i) {
  gost(queries[[i]], custom_bg = custom_backgrounds[[i]])
})
names(results)=names(custom_backgrounds)

# Initialize a list to store the summary of lengths
source_lengths <- sapply(names(results), function(name) {
  if (is.null(results[[name]])) {
    return(0)  # Return 0 if the result is NULL
  } else if (!is.null(results[[name]]$result)) {
    return(length(results[[name]]$result$source))  # Return the length of the 'source' column
  } else {
    return(0)  # Return 0 if 'result' itself is NULL
  }
})

# Print the summary
source_lengths
# Print the resulting data frame
print(source_lengths_df)

####################
#Make plots to show how background correction changed the amount of results
#####################
all_names <- unique(c(names(source_lengths), names(source_lengths_nobkgd)))

# Create aligned vectors, filling with NA where names are missing
source_lengths_aligned <- source_lengths[all_names]
source_lengths_nobkgd_aligned <- source_lengths_nobkgd[all_names]

# Merge into a data frame
source_lengths_df <- data.frame(
  background = source_lengths_aligned,
  no_background = source_lengths_nobkgd_aligned,
  row.names = all_names
)

# Print the resulting data frame
print(source_lengths_df)
# Convert to long format using tidyr::pivot_longer
# Convert rownames to a column before pivoting
df_long <- source_lengths_df %>%
  rownames_to_column(var = "cell_type") %>%  # Add rownames as a new column called "cell_type"
  pivot_longer(cols = c(background, no_background), 
               names_to = "category", 
               values_to = "value")

# Create the ggplot bar plot
p1=ggplot(df_long, aes(x = cell_type, y = value, fill = factor(category,levels = c("no_background","background")))) +
  geom_bar(stat = "identity", position = "dodge") +  # Use dodge to split bars by category
  theme_linedraw() +
  scale_fill_manual(values = c("gray55","green4"))+
  labs(x = "Cell Types", y = "Values", fill = "Category") +  # Axis labels and legend
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none") 


####modified to give percent remaining after bckgd adjustment
df_long <- source_lengths_df %>%
  rownames_to_column(var = "cell_type") %>%
  pivot_longer(cols = c(background, no_background), 
               names_to = "category", 
               values_to = "value")

# Pivot the dataframe back to wide to calculate background/no_background ratio
df_wide <- df_long %>%
  pivot_wider(names_from = category, values_from = value)

# Calculate the ratio of background to no_background
df_wide <- df_wide %>%
  mutate(ratio = background / no_background)

# Now pivot back to long format for plotting (only for the ratio column)
df_long_ratio <- df_wide %>%
  select(cell_type, ratio) %>%  # Keep only cell_type and the ratio
  pivot_longer(cols = ratio, names_to = "category", values_to = "value")  # Pivot the ratio column

# Create the bar plot for the ratio
p2=ggplot(df_long_ratio, aes(x = cell_type, y = value)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use dodge to keep the same format as before
  theme_linedraw() +
  ylim(c(0,1))+
  labs(x = "Cell Types", y = "Ratio (background / no_background)", fill = "Category") +  # Labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability

ggarrange(p1,p2,ncol = 1)

saveRDS(results,file = "./project/brain_aging_phase1/results/multi_query_gprofiler_wBackground.rds")
#########reorganize results into one df
combined_results <- do.call(rbind, lapply(names(results), function(name) {
  result_df <- results[[name]]$result  # Extract the results data frame
  result_df$query <- name                    # Add a new column with the query name
  return(result_df)                         # Return the modified data frame
}))

expected_columns <- colnames(result_df)

combined_results <- do.call(rbind, lapply(names(results), function(name) {
  if (is.null(results[[name]]$result)) {
    # If result is NULL, create a row of NA values
    na_row <- as.data.frame(matrix(NA, nrow = 1, ncol = length(expected_columns)))
    colnames(na_row) <- expected_columns
    na_row$query <- name  # Add the query name
    return(na_row)        # Return the NA row
  } else {
    # Extract the results data frame
    result_df <- results[[name]]$result
    result_df$query <- name  # Add a new column with the query name
    return(result_df)        # Return the modified data frame
  }
}))

#########################
# RRVGO to organize GO terms
##########################
#####Organize GO terms
# Load required libraries
library(gprofiler2)
library(rrvgo)

# Example gprofiler2 result (replace with your actual results)
# results <- gost(query = "your_query_here", organism = "hsapiens") # Uncomment and run your gprofiler query
results=combined_results[combined_results$source%in%c("GO:BP","GO:CC","GO:MF"),]
####################################################################
#WARNING: this disregardles $query, as this is a multi-query object#
####################################################################
#to do iteratively run next line of code
#####################################################################
##results=combined_results[(combined_results$source%in%c("GO:BP","GO:CC","GO:MF"))&(combined_results$query == "CELL_TYPE_OF_INTEREST),]
# For demonstration, let's assume 'results' is your gprofiler output
# Extract the relevant columns from gprofiler results
# Adjust the indexing based on your actual results structure
gene_ids <- results$term_id  # or another column with gene identifiers
p_values <- results$p_value   # P-values associated with each term

# Combine into a data frame
gprofiler_input <- data.frame(
  gene_id = gene_ids,
  p_value = p_values
)

# Optionally, you might want to filter the terms based on significance
gprofiler_input <- gprofiler_input[gprofiler_input$p_value < 0.05, ]  # Example threshold

# Prepare input for rrvgo
# rrvgo expects a data frame with two columns: 'term' and 'score'
# You can create a score based on the p-value (e.g., -log10(p_value))
gprofiler_input$score <- -log10(gprofiler_input$p_value)  # Transform p-value to score
## if you want to score by BH padj, need to specify that instead

# Rename columns to fit rrvgo expectations
rrvgo_input <- gprofiler_input[, c("gene_id", "score")]
colnames(rrvgo_input) <- c("term", "score")  # Rename for rrvgo

# Print the final input for rrvgo
print(rrvgo_input)

# Optionally, you can run rrvgo
# reduced_terms <- rrvgo(rrvgo_input, method = "mean", semdata = your_semdata)
# This line is commented out since you need to replace 'your_semdata' with your actual semantic data
simMatrix <- calculateSimMatrix(rrvgo_input$term,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
###WARNING: This will take some time and scale with size of input terms

scores <- setNames(rrvgo_input$score, rrvgo_input$term)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=2)

scatterPlot(simMatrix, reducedTerms)
#how to make interactive
x=scatterPlot(simMatrix, reducedTerms)
ggplotly(x)

treemapPlot(reducedTerms,size = "size")
