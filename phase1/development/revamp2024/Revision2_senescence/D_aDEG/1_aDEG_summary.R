library(fgsea)
library(tidyverse)
library(scales)

setwd("/data/ADRD/brain_aging/")
out_dir <- "./exploration/revision2_senpaths/D_aDEG/"

# 1. Load Pathways 
sen_paths <- readRDS("./exploration/revision2_senpaths/SenPaths_completed_final.rds")
sen_paths <- sen_paths %>%
  keep(~ !is.null(.x) && length(.x) > 0 && any(.x != ""))

# Create pathway labels with sizes: "Pathway Name (Size)"
path_sizes <- map_int(sen_paths, length)
pathway_label_map <- setNames(
  paste0(names(path_sizes), " (", path_sizes, ")"),
  names(path_sizes)
)

# 2. Load and process aDEG results
res <- read.csv("./phase1/results/aging.glmmtmb_age_diffs_fdr.csv.gz")

res <- res %>%
  filter(type == "region_broad_celltype") %>%
  mutate(
    region = case_when(
      str_detect(tissue, "Middle temporal gyrus") ~ "MTG",
      str_detect(tissue, "Putamen")               ~ "PUT",
      str_detect(tissue, "Entorhinal cortex")      ~ "EC",
      str_detect(tissue, "Subventricular zone")   ~ "SVZ",
      TRUE                                        ~ NA_character_
    ),
    ct = case_when(
      str_detect(tissue, "Middle temporal gyrus") ~ str_remove(tissue, "Middle temporal gyrus "),
      str_detect(tissue, "Putamen")               ~ str_remove(tissue, "Putamen "),
      str_detect(tissue, "Entorhinal cortex")      ~ str_remove(tissue, "Entorhinal cortex "),
      str_detect(tissue, "Subventricular zone")   ~ str_remove(tissue, "Subventricular zone "),
      TRUE                                        ~ tissue
    )
  ) %>%
  filter(!(ct == "SPN" & !region == "PUT"))

# 3. Calculate Background and Overlaps
bg_degs <- res %>%
  group_by(region, ct) %>%
  summarise(total_de_in_group = n(), .groups = "drop")

overlap_summary <- map_df(names(sen_paths), function(p_name) {
  path_genes <- sen_paths[[p_name]]
  
  res %>%
    filter(feature %in% path_genes) %>%
    group_by(region, ct) %>%
    summarise(
      n_path_de = n(),
      pct_of_pathway = (n() / length(path_genes)) * 100,
      .groups = "drop"
    ) %>%
    mutate(pathway = p_name)
}) %>%
  left_join(bg_degs, by = c("region", "ct")) %>%
  mutate(
    # Create labels for the axes
    pathway_label = pathway_label_map[pathway],
    ct_label = paste0(ct, " (", total_de_in_group, ")")
  )

# 4. Visualization
library(ggplot2)

p <- ggplot(overlap_summary, aes(x = ct_label, y = pathway_label)) +
  geom_tile(aes(fill = pct_of_pathway), color = "white") +
  # Shows only the number of overlapping DEGs
  geom_text(aes(label = n_path_de), size = 3, color = "black", family = "Arial") + 
  facet_grid(~region, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low = "#f7f7f7", 
    high = "#e64a02", 
    name = "% Pathway\nDE",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    oob = scales::squish
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Senescence Pathway Coverage in DEGs",
    subtitle = "Tile labels: n overlapping DEGs | Axis labels: (Total size)",
    x = NULL, y = NULL
  )

# Save to PDF
# Using cairo_pdf ensures Arial is embedded correctly

# http://localhost:35241/graphics/plot_zoom?width=1425&height=627&scale=1
png("./exploration/revision2_senpaths/D_aDEG/sen_pathway_coverage_heatmap_final.png", width = 14.25, height = 6.27,units = "in",res = 600)
print(p)
dev.off()

# 5. Statistical Enrichment (Fisher's Exact Test)
# Define the universe as all unique genes tested in the analysis
gene_universe <- unique(res$feature)
n_universe <- length(gene_universe)

stats_results <- map_df(names(sen_paths), function(p_name) {
  path_genes <- sen_paths[[p_name]]
  path_size <- length(path_genes)
  
  # Intersection of pathway genes with the universe 
  # (Important: only test genes that were actually measurable)
  path_in_universe <- intersect(path_genes, gene_universe)
  m <- length(path_in_universe) # Genes in pathway
  n <- n_universe - m          # Genes NOT in pathway
  
  bg_degs %>% 
    rowwise() %>%
    mutate(
      pathway = p_name,
      pathway_size = path_size,
      # Find overlap for this specific group
      overlap_genes = list(intersect(
        res %>% filter(region == .env$region, ct == .env$ct) %>% pull(feature),
        path_in_universe
      )),
      overlap = length(overlap_genes),
      de_size = total_de_in_group,
      
      # Fisher's Exact Test 
      # Matrix structure:
      #            In Path    Not In Path
      # DE         overlap    de_size - overlap
      # Not DE     m-overlap  n - (de_size - overlap)
      p = fisher.test(matrix(c(
        overlap, m - overlap, 
        de_size - overlap, n - (de_size - overlap)
      ), nrow = 2), alternative = "greater")$p.value
    )
}) %>%
  group_by(pathway) %>%
  mutate(fdr_pval = p.adjust(p, method = "fdr")) %>%
  ungroup()

# Final Table Formatting
final_stats_table <- stats_results %>%
  select(ct, region, de_size, pathway_size, pathway, overlap, p, fdr_pval) %>%
  arrange(fdr_pval)

# Save the stats table
write.csv(final_stats_table, 
          file = paste0(out_dir, "sen_pathway_enrichment_stats.csv"), 
          row.names = FALSE)

# Print a preview
print(head(final_stats_table))

