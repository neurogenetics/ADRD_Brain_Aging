library(tidyverse)
#Build Pathway Lists


#Load in longer gene lists
cellage<-read.delim("./exploration/scripts/Phase1_revision_sen/cellage3.tsv",sep = "\t")#https://genomics.senescence.info/download.html#cellage; downloaded 20260311 - CITE PMID37933854 (HAGR v 5)
cellage%>%
  group_by(Type.of.senescence)%>%
  summarize(n=n())
dekhordi<-read.csv("./exploration/scripts/Phase1_revision_sen/Dekhordi_2021_genelists.csv")

dekhordi%>%
  group_by(pathway)%>%
  summarize(n=n())

SenPaths <-list(
  "Hu_2021"= c("CDKN2A", "CDKN1A", "CDKN2D", "CASP8", "IL1B", "GLB1", "SERPINE1"), #"Custom senescence signature, expected to incr."
  "CellAge_HAGRv5_OncoInduced" = cellage%>%filter(Type.of.senescence == "Oncogene-induced")%>%pull(Gene.symbol),
  "CellAge_HAGRv5_Replicative" = cellage%>%filter(Type.of.senescence == "Replicative")%>%pull(Gene.symbol),
  "CellAge_HAGRv5_Stress" = cellage%>%filter(Type.of.senescence == "Stress-induced")%>%pull(Gene.symbol),
  "CellAge_HAGRv5_Unclear" = cellage%>%filter(Type.of.senescence == "Unclear")%>%pull(Gene.symbol),
  "Dehkordi_2021_Canonical" = dekhordi%>%filter(pathway=="Canonical")%>%pull(gene),
  "Dehkordi_2021_Response" = dekhordi%>%filter(pathway=="SRP")%>%pull(gene),
  "Dehkordi_2021_Initiating"= dekhordi%>%filter(pathway=="SIP")%>%pull(gene),
  "Dehkordi_2021_CellAge"= dekhordi%>%filter(pathway=="CellAge")%>%pull(gene),
  "Casella_2019_6C"= c("SLCO2B1","CLSTN2","PTCHD4","LINC02154","PURPL"),#instead of up/down --> use 5 genes from 6C which SLCO2B1, CLSTN2 and PTCHD4 mRNAs, as well as LINC02154 and PURPL lncRNAs in the analysis.discriminated scenceses irrespective of CT'
  "SenNet_2024_CNS_all"= c("CDKN2A", "CDKN1A", "IL6", "TNF", "H2AX", "IL1A", "IL1B", "GLB1", "SERPINE1", "TGFB1","TP53", "CCL2", "CCL5", "CXCL1", "CXCL8", "MMP12", "CCL3", "HMGB1", "MMP3", "LMNB1", "BCL2", "CCL4", "CDKN2B", "CSF1", "IGF1", "TIMP2", "PLAUR", "SPP1"),
  "SenNet_2024_CNS_RNA"= c("CDKN2A", "CDKN1A", "IL6", "TNF", "IL1A","IL1B", "SERPINE1","TGFB1", "TP53", "CCL2", "CCL5", "CXCL1", "CXCL8", "MMP12", "CCL3", "HMGB1", "MMP3", "LMNB1", "BCL2", "CCL4", "CSF1","IGF1", "TIMP2", "PLAUR", "SPP1" ),
  "SenNet_2024_CNS_sc"= c("CDKN2A", "CDKN1A", "IL6", "TNF", "IL1A", "IL1B", "SERPINE1","TGFB1", "TP53", "CCL2", "CCL5", "CXCL1", "CXCL8", "MMP12", "CCL3", "HMGB1", "MMP3", "LMNB1", "BCL2", "CCL4", "CSF1","IGF1", "TIMP2", "PLAUR", "SPP1" )
)

#"SenNet_2024_CNS_all": All genes within SenNet consensus that are reported to be in CNS
#"SenNet_2024_CNS_RNA": Genes within SenNet consensus that are reported to be in CNS + detected with ANY RNA-seq modality
#"SenNet_2024_CNS_sc": Genes within SenNet consensus that are reported to be in CNS + detected with sc/snRNA-Seq

saveRDS(SenPaths, file="/data/ADRD/brain_aging/exploration/revision2_senpaths/SenPaths_completed.rds")

###############################
###############################
###############################
names(SenPaths)
poi<-c("CellAge_HAGRv5_OncoInduced","CellAge_HAGRv5_Replicative","CellAge_HAGRv5_Stress","CellAge_HAGRv5_Unclear","Dehkordi_2021_Canonical","Dehkordi_2021_Response","Dehkordi_2021_Initiating","Dehkordi_2021_CellAge","SenNet_2024_CNS_all")
SenPaths_final<-SenPaths[poi]
saveRDS(SenPaths_final, file="/data/ADRD/brain_aging/exploration/revision2_senpaths/SenPaths_completed_final.rds")

###############################
###############################
###############################

library(ComplexUpset)
library(tidyverse)

# 2. Convert the list of gene vectors into a binary presence/absence data frame
# This identifies every unique gene across all lists and checks its presence in each
all_genes <- unique(unlist(SenPaths_final))
upset_data <- data.frame(gene = all_genes)

for (path_name in names(SenPaths_final)) {
  upset_data[[path_name]] <- as.integer(all_genes %in% SenPaths_final[[path_name]])
}

# 3. Define the sets to plot
pathway_names <- names(SenPaths_final)

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

# 5. Save the plot
cairo_pdf("/data/ADRD/brain_aging/exploration/revision2_senpaths/misc/SenPaths_UpSet_Overlap_final.pdf", width = 16, height = 10)
print(p_upset)
dev.off()
