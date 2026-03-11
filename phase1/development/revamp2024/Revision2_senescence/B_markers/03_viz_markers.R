library(ggplot2)
library(tidyverse)

setwd("/data/ADRD/brain_aging/")
out_dir <- "./exploration/revision2_senpaths/B_marker//"

# Load results
df_glob <- readRDS(paste0(out_dir, "fgsea_res_perCT.rds")) %>% mutate(group = "Global")
df_reg <- readRDS(paste0(out_dir, "fgsea_res_perCTperRegion.rds")) %>% mutate(group = "Regional")

df_all <- rbind(df_glob,df_reg)
df_all$FDR <- p.adjust(df_all$pval,method = "BH")

df_glob<-df_all%>%  
  filter(group=="Global")
df_reg<-df_all%>%  
  filter(group=="Regional")

plot_data <- df_glob %>%
  mutate(log10padj = -log10(FDR),
         significant = FDR < 0.05)

p1 <- ggplot(plot_data, aes(x = test_name, y = pathway)) +
  geom_point(aes(size = log10padj, fill = NES, color = significant), shape = 21, stroke = 0.5) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  # facet_wrap(~group, scales = "free_x") +
  theme_bw(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Senescence Pathway Enrichment", size = "-log10(FDR)", fill = "NES")

# http://localhost:35241/graphics/plot_zoom?width=767&height=485&scale=1
cairo_pdf(paste0(out_dir, "senescence_fgsea_dotplot_global_final.pdf"), width = 7.67, height = 4.86)
print(p1)
dev.off()

plot_data2 <- df_reg %>%
  mutate(log10padj = -log10(FDR),
         significant = FDR < 0.05)%>%
  separate(test_name, into = c("region", "ct"), sep = "_", remove = FALSE) %>%
  # 2. Rename the regions to your specific abbreviations
  mutate(region = case_match(region,
                             "Entorhinal cortex"     ~ "EC",
                             "Middle temporal gyrus" ~ "MTG",
                             "Putamen"               ~ "PUT",
                             "Subventricular zone"   ~ "SVZ",
                             .default = region # Keeps the name if it doesn't match the above
  ))%>%
  filter(!(ct=="SPN"&!region=="PUT"))
  

p2 <- ggplot(plot_data2, aes(x = ct, y = pathway)) +
  geom_point(aes(size = log10padj, fill = NES, color = significant), shape = 21, stroke = 0.5) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  facet_wrap(~region, scales = "free_x") +
  theme_bw(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Senescence Pathway Enrichment", size = "-log10(FDR)", fill = "NES")

# http://localhost:35241/graphics/plot_zoom?width=906&height=695&scale=1
#7.4 x 8
cairo_pdf(paste0(out_dir, "senescence_fgsea_dotplot_region_final.pdf"), width = 9, height = 8)
print(p2)
dev.off()

