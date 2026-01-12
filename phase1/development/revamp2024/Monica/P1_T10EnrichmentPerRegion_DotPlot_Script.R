# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(RColorBrewer)
  # ggnewscale is essential for defining multiple color scales
  library(ggnewscale) 
})

# Define the file path for your Excel workbook
excel_file <- "SuppTable7_SignifEnrich_BckgrCorr.xlsx"

# --- 2. Read and Prepare Data ---

# Get all sheet names (which represent the cell types)
sheet_names <- excel_sheets(excel_file)

# Read all sheets into a single dataframe, adding a 'cell_type' column.
master_df <- sheet_names %>%
  set_names() %>%
  map_dfr(~ read_excel(excel_file, sheet = .x) %>%
            mutate(cell_type = .x),
          .id = "sheet_name_col") %>%
  # Filter for only significant results
  filter(significant == TRUE) %>%
  
  # Select and process columns
  # Note: term_size is included here for the fraction calculation
  select(cell_type, brain_region, term_name, p_value, source, 
         intersection_size, query_size, term_size, precision) %>%
  
  # Calculate -log10(p_value) for color scale
  mutate(neg_log10_p_value = -log10(p_value),
         # !!! FIX: Changed fraction to intersection_size / term_size !!!
         path_label_with_frac = paste0(term_name, " (", 
                                       intersection_size, "/", term_size, ")"))

# Ensure factors for consistent ordering/colors
master_df$source <- as.factor(master_df$source)
master_df$cell_type <- as.factor(master_df$cell_type)

# --- 3. Define Global Variables (Sources, Regions, Colors) ---

# Get all unique sources and brain regions from the entire dataset
all_sources <- unique(master_df$source)
all_brain_regions <- unique(master_df$brain_region)
# Define the order of brain regions
all_brain_regions <- sort(all_brain_regions) 

# Create a consistent color palette for all sources
source_colors <- setNames(
  brewer.pal(length(all_sources), "Dark2")[1:length(all_sources)],
  all_sources
)

# --- 4. Process Data for Plotting (Top 10 per Region) ---

top_enrichments_df <- master_df %>%
  group_by(cell_type, brain_region) %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  ungroup()

# --- 5. Generate and Save Plots for Each Cell Type ---

top_enrichments_nested <- top_enrichments_df %>%
  group_by(cell_type) %>%
  nest(data = -cell_type)

# Loop through each cell type and create the plot
walk2(top_enrichments_nested$data, top_enrichments_nested$cell_type, function(df, ct) {
  
  # Set the Y-axis factor levels for consistent ordering
  df$path_label_with_frac <- factor(df$path_label_with_frac, 
                                    levels = rev(unique(df$path_label_with_frac)))
  
  # --- CRITICAL STEP: Get Source Color Mapping for Y-Axis Text ---
  term_to_source_map <- df %>%
    select(path_label_with_frac, source) %>%
    distinct() %>%
    arrange(path_label_with_frac) 
  
  y_axis_source_colors <- term_to_source_map %>%
    filter(path_label_with_frac %in% levels(df$path_label_with_frac)) %>%
    pull(source) %>%
    as.character()
  
  y_axis_color_map <- source_colors[y_axis_source_colors]
  
  # --- FIX 1: Dummy Data for all Brain Regions (Ensures X-axis is complete) ---
  unique_path_labels <- unique(df$path_label_with_frac)
  dummy_df_regions <- expand.grid(
    brain_region = all_brain_regions,
    path_label_with_frac = unique_path_labels,
    stringsAsFactors = FALSE
  )
  
  # --- FIX 2: Dummy Data for All Sources (Ensures Source Legend is complete) ---
  dummy_df_sources <- data.frame(
    brain_region = factor(all_brain_regions[1], levels = all_brain_regions),
    path_label_with_frac = factor(unique_path_labels[1], levels = levels(df$path_label_with_frac)),
    source = factor(all_sources, levels = all_sources)
  )
  
  # --- Create the plot ---
  p <- ggplot() +
    
    # Layer 1: Dummy points for all regions/terms to force X-axis display
    geom_point(data = dummy_df_regions, 
               aes(x = factor(brain_region, levels = all_brain_regions), 
                   y = factor(path_label_with_frac, levels = levels(df$path_label_with_frac))),
               alpha = 0, size = 0.01) + 
    
    # Layer 2: Dummy points for ALL sources (For Source Legend Keys - Uses DISCRETE color)
    geom_point(data = dummy_df_sources,
               aes(x = brain_region, y = path_label_with_frac, color = source),
               alpha = 0, size = 0) +
    
    # Scale 1: DISCRETE scale for the 'source' aesthetic (Layer 2)
    scale_color_manual(name = "Source",
                       values = source_colors,
                       breaks = all_sources,
                       drop = FALSE,
                       guide = guide_legend(override.aes = list(
                         size = 3,       
                         alpha = 1,      
                         shape = 16,     
                         color = source_colors
                       ), order = 3)) + 
    
    # Reset the 'color' aesthetic for continuous data
    new_scale_color() + 
    
    # Layer 3: Actual data points (Uses CONTINUOUS color)
    geom_point(data = df, 
               aes(x = factor(brain_region, levels = all_brain_regions), 
                   y = path_label_with_frac, 
                   color = neg_log10_p_value, # Continuous color aesthetic
                   size = precision)) +
    
    # Scale 2: CONTINUOUS scale for the 'neg_log10_p_value' aesthetic (Layer 3)
    scale_color_viridis_c(name = expression("-log"[10]*"(p_val)"), direction = 1,
                          guide = guide_colorbar(order = 1)) + 
    
    # Scale for size
    scale_size_continuous(name = "Precision Value", range = c(0.5, 6)) +
    
    # Theme Customization
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = 'black', size = 8),
      
      # Y-axis: Text color for pathway names is mapped to their source color
      axis.text.y = element_text(angle = 0, hjust = 1, size = 8, 
                                 color = y_axis_color_map),
      
      panel.grid.major = element_line(size = 0.25, color = "lightgray"),
      panel.grid.minor = element_blank(),
      
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
    ) +
    
    # Axis Titles
    xlab("Brain Region") +
    ylab("Enriched Pathway (Intersection/Pathway Set Size)") +
    ggtitle(paste(ct, " Top 10 Enrichments by Region"))
    
    # Preview Plot
    print(p)
  
  # Save the plot (CHANGE SIZE APPROPRIATELY)
  #ggsave(paste0(ct, "_Enrichment_DotPlot_T10perRegion.pdf"), p, 
         #width = 6, 
         #height = 6, 
         #dpi = 300)
})
