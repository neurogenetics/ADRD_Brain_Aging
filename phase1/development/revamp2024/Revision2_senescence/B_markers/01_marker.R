#libraries
library(Seurat)
library(tidyverse)

setwd("/data/ADRD/brain_aging/")

#Load in senescence pathways as named list of genes - update from 1_variable_feat when Monica finishes; turn it into an RDS and put just hte name in a 0_curation_list
# sen_paths <- list(
#   "Hu_2021"= c(), 
#   "Dehkordi_2021"= c(""), #this may be a separate canonical (CSP); response (SRP);insitiating (SIP)
#   "Casella_2019_canonical" = c(""),
#   "Casella_2019_response" = c(""),
#   "Casella_2019_initiating"= c(""),
#   "Casella_2019_CellAge"= c(""),
#   "Casella_2019_6C"= c("SLCO2B1","CLSTN2","PTCHD4","LINC02154","PURPL"),#instead of up/down --> use 5 genes from 6C which SLCO2B1, CLSTN2 and PTCHD4 mRNAs, as well as LINC02154 and PURPL lncRNAs in the analysis.discriminated scenceses irrespective of CT'
#   "SenNet_2024_CNS_all"= c(""),
#   "SenNet_2024_CNS_RNA"= c(""),
#   "SenNet_2024_CNS_snRNA"= c("")
# )


#load object
obj = readRDS("./exploration/aging.pegasus.leiden_085.subclustered.rds")

###########################################
###########################################
###########################################
###########################################
#get markers per region per ct
cell_counts<-obj@meta.data %>%
  group_by(broad_celltype,Brain_region)%>%
  summarize(n=n()) # I want markers for just broad_celltype (all) then markers for broad_celltype per region
cell_counts

#per ct
# Define your specific gene list: to do GSEA it has to be ALL features. this may require swarm/batch...
my_genes <- rownames(obj@assays$RNA) # all genes in obj; optimized with just varfeat... not suitable for GSEA

Idents(obj) <- "broad_celltype"
cell_types <- levels(Idents(obj))

global_full_list <- list()

for (ct in cell_types) {
  # FindMarkers for ONE cell type vs ALL others
  # Settings thresholds to 0 ensures NO genes are filtered out
  stats <- FindMarkers(
    object = obj,
    ident.1 = ct,
    features = my_genes,
    logfc.threshold = 0,
    min.pct = 0
  )
  
  # Label the metadata
  stats$gene <- rownames(stats)
  stats$cluster <- ct
  
  # Categorize the results
  stats$status <- "nonsig"
  stats$status[stats$p_val_adj < 0.05 & stats$avg_log2FC > 0] <- "up"
  stats$status[stats$p_val_adj < 0.05 & stats$avg_log2FC < 0] <- "down"
  
  global_full_list[[ct]] <- stats
  print(paste0("Done with: ",ct))
}

# Combine into one master table
# all_global_stats <- do.call(rbind, global_full_list)
saveRDS(global_full_list,file="./exploration/revision2_senpaths/B_marker/markers_list_perCT.rds")

#per ct per region
regions <- unique(obj$Brain_region)
region_full_list <- list()

for (r in regions) {
  message("Processing region: ", r)
  obj_sub <- subset(obj, subset = Brain_region == r)
  Idents(obj_sub) <- "broad_celltype"
  
  # Get cell types actually present in this specific region
  cts_in_region <- intersect(Idents(obj_sub), cell_types)
  
  for (ct in cts_in_region) {
    # Skip if there are too few cells to compare
    if(sum(Idents(obj_sub) == ct) < 3) next 
    
    stats <- FindMarkers(
      object = obj_sub,
      ident.1 = ct,
      features = my_genes,
      logfc.threshold = 0,
      min.pct = 0
    )
    
    stats$gene <- rownames(stats)
    stats$cluster <- ct
    stats$region <- r
    
    # Categorize
    stats$status <- "nonsig"
    stats$status[stats$p_val_adj < 0.05 & stats$avg_log2FC > 0] <- "up"
    stats$status[stats$p_val_adj < 0.05 & stats$avg_log2FC < 0] <- "down"
    
    # Create a unique key for the list
    region_full_list[[paste(r, ct, sep = "_")]] <- stats
    print(paste0("Done with: ",paste(r, ct, sep = "_")))
  }
}

# Combine into one massive master table
# all_region_stats <- do.call(rbind, region_full_list)

saveRDS(region_full_list,file="./exploration/revision2_senpaths/B_marker/markers_list_perCTperRegion.rds")

###########################################
###########################################
###########################################
#this is how we do fgsea, I want the pathways to be sen_paths and the DE results to be the markers
# res=readRDS("./data/RNAseq/diffexp_Summer2024_nocovar_resultlist_MJFFtoxOnly.rds")
# 
# names(res)
# 
# #main focus
# hnDAM <- readRDS("../glia_across_NDDs/analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
# #other pathways that can be added - LDAM (humanized), Keren-Shaul (humanized), MGnD (humanized), hnDAM (altered with lower cutoffs)
# # hnDAM4
# # hnDAM3
# # hnDAM2
# #others to add
# micro_signatures <- read.csv("../glia_across_NDDs/analysis/microglia/misc_files/micro_DAM_signatures_to_compare.csv", header = T)
# micro_signatures_list <- lapply(micro_signatures, function(x) x[x != ""])
# names(micro_signatures_list) <- names(micro_signatures)
# #LDAM
# hLDAM = micro_signatures_list[["Haney_2024_iMG_APOE44_LDhi_vs_LDlow"]]
# #KS
# hKS = micro_signatures_list[["Keren.Shaul_2017_DAM_mouse"]]
# #MGnD
# hMGnD = micro_signatures_list[["Krasemann_2017_MGnD_mouse"]]
# #Marshe_GPNMB
# marshe_GPNMB = micro_signatures_list[["Marshe_2025_GPNMB.hi_human"]]
# 
# pathway_id <- "hsa04210"
# 
# # retrieve pathway data
# pathway_data <- keggGet(pathway_id)
# 
# # extract Entrez Gene IDs
# genes_raw <- pathway_data[[1]]$GENE
# entrez_ids <- genes_raw[seq(1, length(genes_raw), 2)]
# 
# # convert Entrez IDs to gene symbols
# gene_symbols <- mapIds(org.Hs.eg.db,
#                        keys = entrez_ids,
#                        column = "SYMBOL",
#                        keytype = "ENTREZID",
#                        multiVals = "first")
# 
# fgsea_res <- data.frame()
# 
# set.seed(12345)
# 
# test=names(res)[21]
# 
# for (test in names(res)) {
#   test_df <- res[[test]]
#   
#   if (any(duplicated(test_df$gene))) {
#     test_df <- test_df[!duplicated(test_df$gene), ]
#   }
#   
#   stats <- -log10(test_df$P.Value) * sign(test_df$logFC)
#   names(stats) <- rownames(test_df)
#   
#   fgseares <- tryCatch(
#     fgsea(
#       pathways = list(
#         "apoptosis" = gene_symbols,
#         "hnDAM" = hnDAM,
#         "KerenShaul_2017_DAM" = hKS,
#         "Haney_2024_LDAM" = hLDAM,
#         "Krasemann_2017_MGnD" = hMGnD,
#         "Marshe_2025_GPNMB" = marshe_GPNMB
#       ),
#       stats = stats
#     ),
#     error = function(e) NULL
#   )
#   
#   if (is.null(fgseares)) {
#     fgseares <- data.frame(
#       NES = NA_real_,
#       pval = NA_real_,
#       pathway = NA_character_,
#       stringsAsFactors = FALSE
#     )
#   }
#   
#   fgsea_res <- rbind(
#     fgsea_res,
#     data.frame(
#       test = test,
#       NES = fgseares$NES,
#       pval = fgseares$pval,
#       pathway = fgseares$pathway,
#       stringsAsFactors = FALSE
#     )
#   )
# }
# 
# 
# fgsea_res$FDR <- p.adjust(fgsea_res$pval, method = "BH")
# 
# # ggplot(fgsea_res,aes(y = test,x = NES, color = pathway))+
# #   geom_point()+
# #   facet_grid(~pathway)
# 
# saveRDS(fgsea_res,file = "./results/RNAseq/gsea_DAM_noCovars.rds")
#fgsea for sen pathways in each ct
saveRDS(fgsea_res,file = "./exploration/revision2_senpaths/B_marker/fgsea_res_perCTperRegion.rds")
#fgsea for sen pathways in each ct per region
saveRDS(fgsea_res,file = "./exploration/revision2_senpaths/B_marker/fgsea_res_perCTperRegion.rds")

###########################################
###########################################
###########################################
# #this is how we visualize fgsea, I want a dotplot colored by direction, sized by -log10(p) with pathways on the y, and de tests on the x (celltypes and region_celltypes) - faceted by overall/per_region
# fgsea_df <- readRDS("./analysis/microglia/cluster_characterization/FGSEA/micro_FGSEA_signatures_output_ANNOTATED.rds")
# 
# ##################################################
# 
# micro_cluster_order <- rev(c("Micro_Homeo",  "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
#                              "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
#                              "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN"))
# 
# fgsea_df$cluster <- factor(fgsea_df$cluster, levels = micro_cluster_order)
# 
# ##################################################
# 
# fgsea_df$signature <- gsub(".hi", "-hi", fgsea_df$signature)
# fgsea_df$signature <- gsub("Keren.Shaul", "Keren-Shaul", fgsea_df$signature)
# fgsea_df$signature <- gsub("APP.PS1", "APP-PS1", fgsea_df$signature)
# fgsea_df$signature <- gsub("Martins.Ferreira", "Martins-Ferreira", fgsea_df$signature)
# 
# signatures_order <- c("Sun_2023_c0_Homeo", "Marshe_2025_CX3CR1-hi", "Marshe_2025_GRID2-hi",
#                       "Martins-Ferreira_2025_c3_Ribo.DAM1", "Martins-Ferreira_2025_c5_Ribo.DAM2", "Mancuso_2024_Xenograft_Ribo",
#                       "Gerrits_2021_c7_AD1", "Gerrits_2021_c9_AD1", "Gerrits_2021_c10_AD1",
#                       "Sun_2023_c4_Lipid", "Martins-Ferreira_2025_c6_Lipo.DAM", "Marshe_2025_GPNMB-hi", "Silvin_2022_Human_DAMs",
#                       "Dolan_2023_iMGL_c8_DAM", "Mancuso_2024_Xenograft_DAM", "Keren-Shaul_2017_5xFAD_DAM", "Krasemann_2017_APP-PS1_MGnD",
#                       "Sun_2023_c6_Stress", "Sun_2023_c7_Glyco", "Marshe_2025_Senescence", "Marshe_2025_Stress", "Martins-Ferreira_2025_c2_DIMs",
#                       "Mancuso_2024_Xenograft_CRM1", "Mancuso_2024_Xenograft_CRM2",
#                       "Chen_2024_PCDH9-hi",
#                       "Sun_2023_c10_Inflamm3",
#                       "Sun_2023_c5_Phago",
#                       "Sun_2023_c12_Cycling", "Dolan_2023_iMGL_c6_Prolif", "Dolan_2023_iMGL_c9_Prolif", "Dolan_2023_iMGL_c10_Prolif")
# 
# fgsea_df$signature <- factor(fgsea_df$signature, levels = signatures_order)
# 
# ##################################################
# 
# fgsea_df$significant[is.na(fgsea_df$significant)] <- FALSE
# fgsea_df <- na.omit(fgsea_df)
# 
# ##################################################
# 
# fgsea_df$log10padj_squished <- pmin(pmax(fgsea_df$log10padj, 0), 50)
# fgsea_df$NES_squished <- pmin(pmax(fgsea_df$NES, -3), 3)
# 
# ##################################################
# 
# p1 <- ggplot(fgsea_df, aes(x = cluster, y = signature)) +
#   geom_point(aes(size = log10padj_squished, fill = NES_squished, color = significant), shape = 21, stroke = 0.5) +
#   scale_size(range = c(1, 10), limits = c(0, 50)) +
#   scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") + 
#   scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
#         text = element_text(family = "Arial"),
#         axis.text.y = element_text(size = 12)) +
#   labs(x = NULL, y = NULL, size = "-log10(padj)", color = "Significant", fill = "Normalized Enrichment") +
#   coord_flip()
# 
# svglite("./analysis/microglia/plots_final/fgsea_against_clusters.svg", width = 14, height = 5.5)
# plot(p1)
# dev.off()

#viz for sen pathways in each ct
pdf("./exploration/revision2_senpaths/B_marker/fgsea_res_perCTperRegion.pdf", width = 14, height = 5.5)
plot(p1)
dev.off()
#viz for sen pathways in each ct per region
pdf("./exploration/revision2_senpaths/B_marker/fgsea_res_perCTperRegion.pdf", width = 14, height = 5.5)
plot(p1)
dev.off()
