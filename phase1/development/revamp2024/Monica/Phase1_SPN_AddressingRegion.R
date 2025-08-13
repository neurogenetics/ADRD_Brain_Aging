setwd("/data/acridj/project/brain_aging_phase1/")
obj = readRDS("./data/aging.pegasus.leiden_085.subclustered.rds")
levels(obj$broad_celltype)
summary(factor(obj$broad_celltype))
obj = subset(obj, broad_celltype == "SPN")
summary(factor(obj$Brain_region))
#
library(Seurat)
library(harmony)
library(tidyverse)
#### Subclustering of SPN
DefaultAssay(obj)<-"RNA"
# Normalize data
obj <- NormalizeData(obj)

# Find variable features
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

# Scale data
all_genes <- rownames(obj)
obj <- ScaleData(obj, features = all_genes)

obj <- obj %>%
  RunPCA(assay = "RNA") 

obj <- obj %>%
  RunHarmony(assay.use = "RNA", group.by.vars = c("donor_id", "Brain_region"),#add geno back later
             reduction = "pca", reduction.save = 'harmony', plot_convergence = T, lambda = NULL)# Harmony integration grouped vs orig.ident (encoded unique sample info), geno (genetic background), and diff (differentiation protocol)

##############
# Clustering #
##############
# Run UMAP and clustering on integrated data
obj <- obj %>%
  RunUMAP(reduction = "harmony",reduction.name="harmony_10", dims = 1:10)

obj <- obj %>%
  FindNeighbors(reduction = "harmony",dims = 1:10) %>%
  FindClusters(resolution = c(0.1,0.4,0.8))

# Visualize results
colnames(obj@meta.data)
Idents(obj)<-"RNA_snn_res.0.1"
DimPlot(obj, reduction = "harmony_10", label = T)+NoLegend()
# DimPlot(obj, reduction = "harmony_10", split.by = "Brain_region",label = T)
# DimPlot(obj, reduction = "harmony_10", group.by = "donor_id")
# DimPlot(obj, reduction = "harmony_10", group.by = "Age_group", split.by = "Brain_region")
DimPlot(obj, reduction = "harmony_10", group.by = "Brain_region", split.by = "Age_group",cols = c("#832883",#EC 
                                                                                                  "#2d79aa", #MTG
                                                                                                  "#a9254c", #PUT
                                                                                                  "#e4b722"#svz
                                                                                                  ))+NoLegend()

library(scCustomize)


#Features Investigation
goi=c("RBFOX3","DRD1","DRD2","ADORA2A")
DotPlot(obj,features = goi)
FeaturePlot(obj, features = goi,reduction = "harmony_10")



####
# markers=FindAllMarkers(obj,only.pos = T)
# saveRDS(markers,file = "./results/phase1_SPN_subcluster_markers.csv")
# avgexp=AverageExpression(obj,layer = "RNA",slot = "data")
# saveRDS(avgexp,file = "./results/phase1_SPN_subcluster_avgexp.rds")
markers <- readRDS(file = "./results/phase1_SPN_subcluster_markers.csv")
avgexp <- readRDS("./results/phase1_SPN_subcluster_avgexp.rds")
avgexp<-avgexp$RNA



        