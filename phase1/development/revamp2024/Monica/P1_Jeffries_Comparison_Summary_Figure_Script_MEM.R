# 1. Install and Load Necessary Packages
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2"))
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set the name of your Excel file
file_name <- "Jeffries_Comp_CT_UpSet_Inputs.xlsx"

# Get the names of all sheets (which represent your cell types)
cell_types <- excel_sheets(file_name)

# Initialize an empty list to store data frames from each sheet
data_list <- list()

# 2. Read and Process Data from Each Excel Sheet
for (ct in cell_types) {
  # Read the current sheet (Cell Type)
  df <- read_excel(file_name, sheet = ct)
  
  # Ensure the first column is named 'gene'
  if (names(df)[1] != "gene") {
    names(df)[1] <- "gene"
  }
  
  # --- Data Aggregation ---
  
  # A. Process Jeffries Data: Count genes with '1' in the Jeffries column
  jeffries_col <- names(df)[grepl("_Jeffries", names(df))]
  
  jeffries_data <- df %>%
    filter(!!sym(jeffries_col) == 1) %>%
    summarise(
      CellType = ct,
      Group_Bar = "Jeffries Study", # Defines the first bar group
      Region_Stack = "Jeffries Study", # Defines the only stack in this bar
      Count = n()
    )
  
  # B. Process Our Study Data: Columns with regional flags
  regional_cols <- names(df)[!grepl("_Jeffries", names(df)) & names(df) != "gene"]
  
  our_study_data <- df %>%
    # Select only gene and regional columns
    select(gene, all_of(regional_cols)) %>%
    # Pivot longer to get a row for each gene/region combination
    pivot_longer(
      cols = all_of(regional_cols),
      names_to = "Region_Raw",
      values_to = "Flag"
    ) %>%
    # Keep only genes flagged as '1' (TRUE)
    filter(Flag == 1) %>%
    # Extract the meaningful region name by removing the cell type suffix (e.g., '_ExN')
    mutate(
      Region_Stack = gsub("_[A-Za-z0-9]+$", "", Region_Raw),
      # Clean up region names for better plotting (optional, but helpful)
      Region_Stack = gsub(" ", "\n", Region_Stack)
    ) %>%
    # Count the number of genes for each region
    group_by(Region_Stack) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    # Add CellType and Group_Bar columns
    mutate(
      CellType = ct,
      Group_Bar = "Our Study" # Defines the second bar group, which will be stacked
    )
  
  # Combine Jeffries and Our Study data for the current cell type
  cell_type_data <- bind_rows(jeffries_data, our_study_data)
  
  # Add the combined data frame to the list
  data_list[[ct]] <- cell_type_data
}

# 3. Combine All Cell Type Data into a Single Data Frame and Apply Custom Order
final_data <- bind_rows(data_list)

# Define the user-specified order for cell types
custom_cell_type_order <- c("InN", "ExN", "OPC", "Oligodendrocyte", "Astrocyte", "Microglia", "Endothelial")

final_data <- final_data %>%
  # Define the order for the Group_Bar factor to ensure Jeffries is always first
  mutate(Group_Bar = factor(Group_Bar, levels = c("Jeffries Study", "Our Study"))) %>%
  # Use the custom order for the CellType factor
  mutate(CellType = factor(CellType, levels = custom_cell_type_order))

# --- NEW STEP: Calculate Total Counts per Bar ---
# This data frame holds the total height for the label position
bar_totals <- final_data %>%
  group_by(CellType, Group_Bar) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = 'drop') %>%
  # Ensure the CellType factor levels are the same as in final_data
  mutate(CellType = factor(CellType, levels = custom_cell_type_order))


# 4. Create the Grouped and Stacked Bar Plot
plot <- ggplot(final_data, aes(x = Group_Bar, y = Count, fill = Region_Stack)) +
  
  # 1. Bar Geometry
  geom_bar(
    stat = "identity",
    position = position_stack(reverse = TRUE),
    width = 0.8
  ) +
  
  # 2. Text Geometry (Total Labels)
  geom_text(
    data = bar_totals, # Use the summary data frame
    aes(x = Group_Bar, y = Total_Count, label = Total_Count), # Map total count to label
    inherit.aes = FALSE, # Do not inherit fill (region) aesthetic
    vjust = -0.5, # Adjust position slightly above the bar
    size = 4,
    fontface = "plain"
  ) +
  
  # Crucial change: Use facet_wrap to group by CellType, which uses the new factor order
  facet_wrap(~ CellType, scales = "free_x", nrow = 1) +
  
  # --- Aesthetics and Labels ---
  
  labs(
    title = "Comparison of aDEG Counts by Cell Type: Jeffries et al. vs. Our Study",
    x = "Study", # X-axis now represents the two studies (Jeffries vs. Our Study)
    y = "Number of aDEGs",
    fill = "Study/Region"
  ) +
  
  # Customize the theme for a clean look
  theme_minimal(base_size = 14) +
  theme(
    # Ensure the x-axis text (Jeffries Study / Our Study) is visible
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Customize the facet strip (the Cell Type label)
    strip.text = element_text(face = "bold", size = 12),
    # Adjust plot title
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # Add a border to the facets to visually separate the cell type groups
    panel.border = element_rect(colour = "gray", fill = NA, linewidth = 0.5),
    # Remove the internal grid lines
    panel.grid.major.x = element_blank()
  ) +
  
  # Manually set colors
  scale_fill_manual(
    values = c(
      "Jeffries Study" = "gray50",
      "Middle\ntemporal\ngyrus" = "#2D79AA",
      "Entorhinal\ncortex" = "#832883",
      "Putamen" = "#a9254c",
      "Subventricular\nzone" = "#e4b722"
      # Include other regions if they appear in your data
    )
  )

# Display the plot
print(plot)
# 5. Save the Plot (optional)
ggsave("Gene_Count_Comparison_Plot.pdf", plot, width = 12, height = 7, dpi = 300)

# 6. Export the Final Data Frame to CSV
write.csv(final_data, "Jeffries_Comp_Final_Gene_Count_Data.csv", row.names = FALSE)
