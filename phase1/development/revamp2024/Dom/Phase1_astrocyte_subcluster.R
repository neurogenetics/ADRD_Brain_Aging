library(Seurat)
library(harmony)
library(SAHA)
setwd("/data/ADRD/brain_aging/")
obj = readRDS("./exploration/aging.pegasus.leiden_085.subclustered.rds")
astro=subset(obj,broad_celltype=="Astrocyte")

# #function to harmony integrat
astro <- NormalizeData(astro)
astro <- FindVariableFeatures(astro, nfeatures = 2000)
astro <- ScaleData(astro, vars.to.regress = "percent_mito")
astro <- RunPCA(astro)
# 
astro <- RunHarmony(object = astro, reduction = "pca", group.by.vars = "pool_name",
                    reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
print(ElbowPlot(astro, reduction = "harmony", ndims = 3))

n_dims=3

astro <- FindNeighbors(astro, reduction = "harmony", dims = 1:n_dims)
astro <- RunUMAP(astro, reduction = "harmony", dims = 1:n_dims, reduction.name = "umap_harmony")

cluster_res=0.1

astro <- FindClusters(astro, resolution = cluster_res)

astro@meta.data |>
  select(donor_id, RNA_snn_res.0.1) |>
  group_by(donor_id, RNA_snn_res.0.1) |>
  tally(name = "n_cells") |>
  group_by(RNA_snn_res.0.1) |>
  summarize(median_per_cluster = median(n_cells))
# 
# # http://localhost:36207/graphics/plot_zoom?width=409&height=304&scale=1
cairo_pdf("./exploration/plots/Phase1_DJArevisions/astro_Subcluster_HarmonyUMAP.pdf",width = 4.09,height = 3.04)
DimPlot(astro,reduction = "umap_harmony",group.by = "RNA_snn_res.0.1")
dev.off()
cairo_pdf("./exploration/plots/Phase1_DJArevisions/astro_Subcluster_HarmonyUMAP_byRegion.pdf",width = 4.09,height = 3.04)
DimPlot(astro,reduction = "umap_harmony",group.by = "Brain_region")
dev.off()

ast_counts=astro@meta.data%>%
  group_by(RNA_snn_res.0.1,Brain_region)%>%
  summarize(n=n())
p_cts=ggplot(ast_counts,aes(x=RNA_snn_res.0.1,y=n,fill=Brain_region))+
  geom_bar(stat="identity")+
  theme_linedraw()

# http://localhost:40569/graphics/plot_zoom?width=580&height=518&scale=1
cairo_pdf("./exploration/plots/Phase1_DJArevisions/Astro_Subcluster_ctsByRegion.pdf",width = 5.80,height = 5.18)
p_cts
dev.off()


# # A tibble: 14 × 3
# # Groups:   RNA_snn_res.0.1 [4]
# RNA_snn_res.0.1 Brain_region              n
# <fct>           <fct>                 <int>
#   1 0               Entorhinal cortex       771
# 2 0               Middle temporal gyrus    83
# 3 0               Putamen                  15
# 4 0               Subventricular zone    3883
# 5 1               Entorhinal cortex      1531
# 6 1               Middle temporal gyrus   540
# 7 1               Putamen                1645
# 8 1               Subventricular zone     120
# 9 2               Entorhinal cortex      1012
# 10 2               Middle temporal gyrus   244
# 11 2               Putamen                 476
# 12 2               Subventricular zone      80
# 13 3               Entorhinal cortex       884
# 14 3               Subventricular zone       5

# 
saha_input=Seurat2SAHA(obj = astro,output = "AvgExp")
query_avgexp=saha_input$avgexp
# 
# 
# #from gliaNDD astroglia
glia_astro_SAHA <- readRDS("/data/acridj/project/glia_metastudies/gliaNDD_astrocyte_SAHAinputs.rds")
# 
ndd_astro_avgexp <- glia_astro_SAHA$avgexp
# ndd_astro_markers <- glia_astro_SAHA$markers
ndd_astro_varfeat <- glia_astro_SAHA$varfeat
# 
ann<-Create_SAHA_object(query = query_avgexp,db = ndd_astro_avgexp,data_type = "AvgExp")
ann=Initialize_MarkerFree(ann = ann)
ann=Downsample(ann,custom_ds = ndd_astro_varfeat)
ann=NormalizeDS(ann,assay_query = "RNA",norm_method = "across_clust")
ann=CorrelateDS(ann,corr_method = "spearman")
ann=Create_MarkerFree_Viz(ann)
call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp")
# 
# http://localhost:36207/graphics/plot_zoom?width=656&height=441&scale=1
cairo_pdf("./exploration/plots/Phase1_DJArevisions/astro_Subcluster_GliaNDDMarkerFree.pdf",width = 6.56,height = 4.41)
call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp")

#from this https://github.com/neurogenetics/SAHA/blob/main/R/CorrelateDS.R
# CorrelateDS <- function(
    # ann, 
corr_method = "spearman"
# ){
#   # compute Pearson coefficients b/w query and ABC celltype expression profiles
db_columns <- colnames(ann@results$marker_free$norm_merge)[grepl("db", colnames(ann@results$marker_free$norm_merge))]
query_columns <- colnames(ann@results$marker_free$norm_merge)[grepl("query", colnames(ann@results$marker_free$norm_merge))]
#   
correlation.df <- matrix(NA, nrow = length(query_columns), ncol = length(db_columns))
rownames(correlation.df) <- query_columns
colnames(correlation.df) <- db_columns
pval.df <- matrix(NA, nrow = length(query_columns), ncol = length(db_columns))
rownames(pval.df) <- query_columns
colnames(pval.df) <- db_columns
#   
for (i in seq_along(query_columns)) {
  for (j in seq_along(db_columns)) {
    suppressWarnings(cor_test <- cor.test(ann@results$marker_free$norm_merge[[query_columns[i]]], ann@results$marker_free$norm_merge[[db_columns[j]]], method = corr_method))
    correlation.df[i, j] <- cor_test$estimate
    pval.df[i,j]<-cor_test$p.value
  }
}
#   
correlation.df <- as.data.frame(correlation.df)
pval.df <- pval.df %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "query_clust") %>%
  pivot_longer(
    cols = -query_clust,
    names_to = "cell_test",
    values_to = "p_val"
  )
pval.df$fdr <- p.adjust(pval.df$p_val,method = "BH")
write.csv(pval.df,"./exploration/plots/Phase1_DJArevisions/Astro_Subcluster_GliaNDDMarkerFree_pval.csv")
#   ann@results$marker_free$corr=correlation.df
#   ann@params$marker_free$corr_method <- corr_method
#   return(ann)
# }
dev.off()