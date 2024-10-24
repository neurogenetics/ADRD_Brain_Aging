################
#Date: 2024/10/23
#Authors: chatGPT, Gemini, Perplexity
#Code Curators: Monica, Amal, Dom
#Goals:
#    - Directionality ct_xRegion_aDEGs: (*PRIORITY*)
#         - generate df for heatmap
#         - describe concordance vs discordance
#         - heatmap that depicts con/dis of aDEG across regions (per celltype)
#
#Problems for later:
#    - Turn into fxn so that this code makes con/dis heatmaps with whatever you feed in
#    - Standaradize tissue/region/ct naming
#    - Recode regions as abbreviations to make them easier to work with
#    - Directionality Region_xCellType_aDEGs:
#         - same as above - but celltype (per region)
################
#libraries
library(tidyverse)
#Dataframe that shows the direction of aDEG for genes that are differentially expressed in >1 region (per celltype)
micro_aDEGs=read.csv("./results/Microglia_ByRegion_aDEG.csv")
df=micro_aDEGs#rename for generalizability
##################
#Functions
classify_combination <- function(val1, val2) { # this defines con/dis-cordance in pairwise fashion
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
#######################
#Get colnames to do combinations of
colnames(df)
  #remove feature and store things to be combined
regions=colnames(df)[2:5]
# Define all possible combinations of columns A, B, C, D
combinations <- combn(regions, 2, simplify = FALSE)

# Summarize the rows (features) by each combination of columns
df_combined <- df %>%
  rowwise() %>%
  mutate(across(all_of(regions), as.character))  # Ensure all columns are character

# Loop through combinations to classify and create new columns
for (comb in combinations) {
  col1 <- comb[1]
  col2 <- comb[2]
  
  # Create a new column name based on the combination
  new_col <- paste0(col1,"_", col2)
  
  # Add a new column to the dataframe with classified combinations
  df_combined <- df_combined %>%
    mutate(!!new_col := classify_combination(!!sym(col1), !!sym(col2)))
}

# View the new dataframe with the classified combination columns
print(df_combined)

############################
# CONCORDENCE###############
############################
# Create a matrix to count uu and dd pairwise
pairwise_counts <- data.frame(matrix(0, nrow = length(regions), ncol = length(regions)))
rownames(pairwise_counts) <- colnames(pairwise_counts) <- regions

# Fill the matrix with counts
for (comb in combn(regions, 2, simplify = FALSE)) {
  col1 <- comb[1]
  col2 <- comb[2]
  
  # Count "uu" and "dd"
  uu_count <- sum(df_combined[[paste0(col1,"_", col2)]] == "uu", na.rm = TRUE)
  dd_count <- sum(df_combined[[paste0(col1, "_",col2)]] == "dd", na.rm = TRUE)
  
  # Assign counts to the appropriate positions
  pairwise_counts[col1, col2] <- uu_count  # Upper triangle
  pairwise_counts[col2, col1] <- dd_count  # Lower triangle
}

# Replace NA with 0 for proper conversion to long format
# Replace NA with 0 for proper conversion to long format
pairwise_counts[is.na(pairwise_counts)] <- 0

# Ensure pairwise_counts is a matrix
pairwise_counts <- as.matrix(pairwise_counts)

###################
# Create a matrix for final counts
# In order to color the triangles red for uu and blue for dd; direction is added to counts +# means uu, -# means dd
final_counts <- matrix(0, nrow = nrow(pairwise_counts), ncol = ncol(pairwise_counts), 
                       dimnames = dimnames(pairwise_counts))

# Fill the upper triangle with positive values (uu)
final_counts[upper.tri(final_counts)] <- pairwise_counts[upper.tri(pairwise_counts)]

# Fill the lower triangle with negative values (dd)
final_counts[lower.tri(final_counts)] <- -pairwise_counts[lower.tri(pairwise_counts)]

# Set diagonal to NA or keep as zero
diag(final_counts) <- NA  # Optionally set to NA or leave as 0

# Convert the matrix to long format for ggplot
suppressWarnings(heatmap_data <- melt(final_counts, na.rm = T))#suppress NA warning mesasge, we read and accept the warning

# Create a heatmap with ggplot2
ggplot(data = heatmap_data, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "white") +  # Use a white border for tiles
  geom_text(aes(label = abs(value)), vjust = 1, hjust = 1) +  # Add text labels
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, na.value = "grey50") +  # Color scale
  theme_minimal() +
  labs(title = "Heatmap of Pairwise Concordant aDEGs within Microglia", 
       x = "", 
       y = "", 
       fill = "Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

############################
# DISCORDENCE###############
############################
# Create a matrix to count du and ud pairwise
pairwise_counts <- data.frame(matrix(0, nrow = length(regions), ncol = length(regions)))
rownames(pairwise_counts) <- colnames(pairwise_counts) <- regions

# Fill the matrix with counts
for (comb in combn(regions, 2, simplify = FALSE)) {
  col1 <- comb[1]
  col2 <- comb[2]
  
  # Count "uu" and "dd"
  du_count <- sum(df_combined[[paste0(col1,"_", col2)]] == "du", na.rm = TRUE)
  ud_count <- sum(df_combined[[paste0(col1, "_",col2)]] == "ud", na.rm = TRUE)
  
  # Assign counts to the appropriate positions
  pairwise_counts[col1, col2] <- du_count  # Upper triangle
  pairwise_counts[col2, col1] <- ud_count  # Lower triangle
}

# Replace NA with 0 for proper conversion to long format
# Replace NA with 0 for proper conversion to long format
pairwise_counts[is.na(pairwise_counts)] <- 0

# Ensure pairwise_counts is a matrix
pairwise_counts <- as.matrix(pairwise_counts)

###################
# Set diagonal to NA or keep as zero
diag(pairwise_counts) <- NA  # Optionally set to NA or leave as 0

# Convert the matrix to long format for ggplot
suppressWarnings(heatmap_data <- melt(pairwise_counts, na.rm = T))

# Create a heatmap with ggplot2
ggplot(data = heatmap_data, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "white") +  # Use a white border for tiles
  geom_text(aes(label = abs(value)), vjust = 1, hjust = 1) +  # Add text labels
  scale_fill_gradient2(low = "white",high = "magenta4", 
                       midpoint = 0, na.value = "grey50") +  # Color scale
  theme_minimal() +
  labs(title = "Heatmap of Pairwise Discordant aDEGs within Microglia", 
       x = "Down aDEGs", 
       y = "Up aDEGs", 
       fill = "Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

############################
# UpSet to Con/DisCord######
############################
###### bar plot showing how many genes have 1,2,3,4 shared overlaps... (linking Upset to Con/Discord)
# Assuming your data frame is named 'df'
df <- df %>%
  mutate(
    non_ne_count = rowSums(select(., -feature) != "NE")
  )

ggplot(df)

# Assuming your data frame is named 'df'
ggplot(df, aes(x = factor(non_ne_count))) +
  geom_bar() +
  labs(
    title = "Distribution of aDEGs Shared Across Regions",
    x = "# Shared Regions",
    y = "# aDEGs"
  )+
  theme_bw()





