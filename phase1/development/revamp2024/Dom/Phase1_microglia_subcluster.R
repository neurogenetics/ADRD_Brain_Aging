library(Seurat)
library(harmony)
library(SAHA)
obj = readRDS("./exploration/aging.pegasus.leiden_085.subclustered.rds")
micro=subset(obj,broad_celltype=="Microglia")
#function to harmony integrat
micro <- NormalizeData(micro)
micro <- FindVariableFeatures(micro, nfeatures = 2000)
micro <- ScaleData(micro, vars.to.regress = "percent_mito")
micro <- RunPCA(micro)

micro <- RunHarmony(object = micro, reduction = "pca", group.by.vars = "pool_name",
                         reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
print(ElbowPlot(micro, reduction = "harmony", ndims = 3))

n_dims=3

micro <- FindNeighbors(micro, reduction = "harmony", dims = 1:n_dims)
micro <- RunUMAP(micro, reduction = "harmony", dims = 1:n_dims, reduction.name = "umap_harmony")

cluster_res=0.1

micro <- FindClusters(micro, resolution = cluster_res)

# http://localhost:36207/graphics/plot_zoom?width=409&height=304&scale=1
cairo_pdf("./exploration/plots/Phase1_DJArevisions/Micro_Subcluster_HarmonyUMAP.pdf",width = 4.09,height = 3.04)
DimPlot(micro,reduction = "umap_harmony",group.by = "RNA_snn_res.0.1")
dev.off()

micro_cts=micro@meta.data%>%
  group_by(RNA_snn_res.0.1,Brain_region)%>%
  summarize(n=n())
p_cts=ggplot(micro_cts,aes(x=RNA_snn_res.0.1,y=n,fill=Brain_region))+
  geom_bar(stat="identity")+
  theme_linedraw()

# http://localhost:40569/graphics/plot_zoom?width=580&height=518&scale=1
cairo_pdf("./exploration/plots/Phase1_DJArevisions/Micro_Subcluster_ctsByRegion.pdf",width = 5.80,height = 5.18)
p_cts
dev.off()



saha_input=Seurat2SAHA(obj = micro,output = "AvgExp")
query_avgexp=saha_input$avgexp


#from gliaNDD microglia
glia_micro_SAHA <- readRDS("/data/acridj/project/glia_metastudies/gliaNDD_microglia_SAHAinputs.rds")

ndd_micro_avgexp <- glia_micro_SAHA$avgexp
ndd_micro_markers <- glia_micro_SAHA$markers
ndd_micro_varfeat <- glia_micro_SAHA$varfeat

ann<-Create_SAHA_object(query = query_avgexp,db = ndd_micro_avgexp,data_type = "AvgExp")
ann=Initialize_MarkerFree(ann = ann)
ann=Downsample(ann,custom_ds = ndd_micro_varfeat)
ann=NormalizeDS(ann,assay_query = "RNA",norm_method = "across_clust")
ann=CorrelateDS(ann,corr_method = "spearman")
ann=Create_MarkerFree_Viz(ann)
call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp")
#add stars?? - not in package
ann@results$marker_free$corr
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
  write.csv(pval.df,"./exploration/plots/Phase1_DJArevisions/Micro_Subcluster_GliaNDDMarkerFree_pval.csv")
#   ann@results$marker_free$corr=correlation.df
#   ann@params$marker_free$corr_method <- corr_method
#   return(ann)
# }



# http://localhost:36207/graphics/plot_zoom?width=656&height=441&scale=1
cairo_pdf("./exploration/plots/Phase1_DJArevisions/Micro_Subcluster_GliaNDDMarkerFree.pdf",width = 6.56,height = 4.41)
call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp")
dev.off()