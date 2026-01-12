setwd("/data/acridj/project/brain_aging_phase1/")
dat = readRDS("./data/aging.pegasus.leiden_085.subclustered.rds")
levels(dat$broad_celltype)
summary(factor(dat$broad_celltype))
obj = subset(dat, broad_celltype == "SPN")
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
png(filename = "/data/ADRD/brain_aging/exploration/plots/Phase1_DJArevisions/SPN_subclustering.png",width = 5.15,height = 4.00,units = "in",res = 600)
DimPlot(obj, reduction = "harmony_10", label = T)+NoLegend()
dev.off()
# DimPlot(obj, reduction = "harmony_10", split.by = "Brain_region",label = T)
# DimPlot(obj, reduction = "harmony_10", group.by = "donor_id")
# DimPlot(obj, reduction = "harmony_10", group.by = "Age_group", split.by = "Brain_region")
obj$Age_group<-factor(obj$Age_group,levels = c("young","old"))
png(filename = "/data/ADRD/brain_aging/exploration/plots/Phase1_DJArevisions/SPN_subclustering_byagebyRegion.png",width = 10,height = 4.00,units = "in",res = 600)
DimPlot(obj, reduction = "harmony_10", group.by = "Brain_region", split.by = "Age_group",cols = c("#832883",#EC 
                                                                                                  "#2d79aa", #MTG
                                                                                                  "#a9254c", #PUT
                                                                                                  "#e4b722"#svz
))
dev.off()
# library(scCustomize)


#Features Investigation
goi=c("RBFOX3","DRD1","DRD2","ADORA2A")
DotPlot(obj,features = goi)
FeaturePlot(obj, features = goi,reduction = "harmony_10")

###########################
#NEW prettier
###########################


#make barplots of cluster by age_group and another by region, then show dotplot
#question is that it is "too much data" to through out, does it show that it really isn't

obj@meta.data%>%
  group_by(broad_celltype,Brain_region)%>%
  summarize(n=n())





#Asked "what are they" to the ones we remove....
library(SAHA)
#SAHA query: wrong_region (3_ids) vs phase1_avgexp
#IF REMOVING PUT FROM SAHA (sub)
# sub<-subset(obj,broad_celltype=="SPN" & Brain_region != "Putamen")
# sub@meta.data%>%
#   group_by(broad_celltype,Brain_region)%>%
#   summarize(n=n())
# Idents(sub)<-"Brain_region"
#IF not REMOVING PUT FROM SAHA (obj)
sub<-subset(dat,broad_celltype=="SPN")
sub@meta.data%>%
  group_by(broad_celltype,Brain_region)%>%
  summarize(n=n())

Idents(sub)<-"Brain_region"
saha_query=Seurat2SAHA(sub,output = "AvgExp")#sub
#get query from unsubset object
Idents(dat)<-"broad_celltype"
dat<-FindVariableFeatures(dat)
saha_db=Seurat2SAHA(dat,output = "AvgExp")
####SAHA of subcortical regions? need to use mouse?

ann<-Create_SAHA_object(query = saha_query$avgexp,db = saha_db$avgexp%>%dplyr::select(-RNA.Other),data_type = "AvgExp")
ann=Initialize_MarkerFree(ann = ann)
ann=Downsample(ann,custom_ds = saha_db$varfeat)
ann=NormalizeDS(ann,assay_query = "RNA",norm_method = "within_clust")
ann=CorrelateDS(ann,corr_method = "pearson")
ann=Create_MarkerFree_Viz(ann)
# http://localhost:41311/graphics/plot_zoom?width=739&height=352&scale=1
cairo_pdf("/data/ADRD/brain_aging/exploration/plots/Phase1_DJArevisions/SPN_SAHA_wrongRegionVsAll.pdf",width = 7.39,height = 3.52)
call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp")
dev.off()


##########################################
#DENDRO
library(ggplot2)
library(ggdendro)
library(tidyverse)
library(Seurat)
library(cowplot) 
library(extrafont)
# library(scCustomize)
library(patchwork)
library(SAHA)

SPNs_meta <- obj@meta.data

saha_spnsub <- Seurat2SAHA(obj,output = "AvgExp")
SPNs_avgexp <- saha_spnsub$avgexp
SPNs_varfeats <- saha_spnsub$varfeat


colnames(SPNs_avgexp) <- gsub("RNA.g", "", colnames(SPNs_avgexp))
colnames(SPNs_avgexp) <- gsub("\\.", " ", colnames(SPNs_avgexp))

##################################################

# hierarchical clustering

SPNs_avgexp_use <- SPNs_avgexp[rownames(SPNs_avgexp) %in% SPNs_varfeats, ]

SPNs_dendro <- dist(scale(t(SPNs_avgexp_use))) %>% 
  hclust(method = "ward.D2") %>% 
  as.dendrogram()

# get the order of the dendrogram
ggdendrogram(SPNs_dendro, rotate = T)

SPNs_subclass_order <- c("1","0","2","6","5","4","3")


SPNs_dendro <- reorder(SPNs_dendro, 
                       match(colnames(SPNs_avgexp_use), rev(SPNs_subclass_order)), 
                       agglo.FUN=mean)

SPNs_ddata <- dendro_data(SPNs_dendro, type = "rectangle")

SPNs_den_p <- ggplot(segment(SPNs_ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) +
  theme_void() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

########

# number of cells, plotted by dataset

SPNs_ncells_data <- SPNs_meta %>%
  group_by(Brain_region,RNA_snn_res.0.1) %>%
  summarize(ncells = n()) %>%
  ungroup()


SPNs_ncells_p <- ggplot(SPNs_ncells_data, aes(x = ncells, y = RNA_snn_res.0.1, fill = Brain_region)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c(`Entorhinal cortex`="#832883",#EC 
                               `Middle temporal gyrus`="#2d79aa", #MTG
                               `Putamen`="#a9254c", #PUT
                               `Subventricular zone`="#e4b722"#svz
                               )) +
  scale_y_discrete(limits = rev(SPNs_subclass_order)) +
  scale_x_continuous(position = "top", limits = c(0, 16000), breaks = c(0, 4000, 8000, 12000, 16000)) +
  theme_linedraw() +
  theme(axis.text.y = element_text(size = 10, face = "bold", hjust = 1), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90,hjust = 0),
        axis.ticks.x = element_line(linewidth = 0.2), axis.line.x = element_line(linewidth = 0.2),
        text = element_text(family = "Arial"),legend.position = "none")

SPNs_cols = c(`Entorhinal cortex`="#832883",#EC 
                      `Middle temporal gyrus`="#2d79aa", #MTG
                      `Putamen`="#a9254c", #PUT
                      `Subventricular zone`="#e4b722"#svz
)

# number of UMIs/features by cluster

SPNs_nUMIs_p <- ggplot(SPNs_meta, aes(x = nCounts_RNA, y = RNA_snn_res.0.1, fill = RNA_snn_res.0.1)) +
  geom_violin() + 
  scale_fill_manual(values = SPNs_cols) + 
  scale_y_discrete(limits = rev(SPNs_subclass_order)) +
  scale_x_continuous(position = "top", limits = c(0, 85000), breaks = c(0, 25000, 50000, 75000)) +
  theme_void() +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 0), 
        axis.ticks.x = element_line(linewidth = 0.2), axis.line.x = element_line(linewidth = 0.2),
        text = element_text(family = "Arial")) +
  stat_summary(fun = median, geom = "point", size = 2, color="black")


SPNs_nfeats_p <- ggplot(SPNs_meta, aes(x = nFeaturess_RNA, y = RNA_snn_res.0.1, fill = RNA_snn_res.0.1)) +
  geom_violin() + 
  scale_fill_manual(values = SPNs_cols) + 
  scale_y_discrete(limits = rev(SPNs_subclass_order)) +
  scale_x_continuous(position = "top") +
  theme_void() +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 0), 
        axis.ticks.x = element_line(linewidth = 0.2), axis.line.x = element_line(linewidth = 0.2),
        text = element_text(family = "Arial")) +
  stat_summary(fun = median, geom = "point", size = 2, color="black")

########

# dotplot - what we had before + some from Monica's + whatever the one weird one is
# devtools::install_github("immunogenomics/presto")
# markers_cluster3 <- FindMarkers(
#   object = obj,
#   ident.1 = 3
# )

SPNs_genes <- c("RBFOX3", "DRD1","DRD2",
                "ADORA2A", "PPM1E"
)

for_legend <- DotPlot(obj, features = SPNs_genes, cols = c("#d9f4fa", "#64288a"), dot.min = 0.05, dot.scale = 6) +
  scale_x_discrete(limits = SPNs_genes, position = "top") +
  scale_y_discrete(limits = rev(SPNs_subclass_order)) +
  theme(legend.position = "right", axis.title = element_blank(), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.2), axis.ticks = element_line(linewidth = 0.2), 
        axis.text.x = element_text(size = 10, angle = 90,hjust = 0), axis.text.y = element_blank(),
        text = element_text(family = "Arial")) 

SPNs_expr_p <-for_legend+
  NoLegend()

########

# combine plots

# http://localhost:41311/graphics/plot_zoom?width=775&height=294&scale=1
cairo_pdf("/data/ADRD/brain_aging/exploration/plots/Phase1_DJArevisions/SPN_cellsummary.pdf",width = 7.75,height = 2.94)
plot_grid(SPNs_den_p, SPNs_ncells_p, SPNs_nfeats_p,SPNs_nUMIs_p, SPNs_expr_p, 
          nrow = 1, 
          ncol = 5, 
          rel_widths = c(0.5, 
                               1.2, 
                               0.5,
                               0.5,
                               0.7), 
          align = "h", axis = "tb")
dev.off()
# http://localhost:41311/graphics/plot_zoom?width=527&height=381&scale=1
cairo_pdf("/data/ADRD/brain_aging/exploration/plots/Phase1_DJArevisions/SPN_cellsummary_legend.pdf",width = 5.27,height = 3.81)
for_legend
dev.off()

###finally a pie diagram of cells removed
library(tidyverse)

# -----------------------------
# Data + Colors + Plot
# -----------------------------

pie_data <- tibble(
  group = c("Other neurons", "SPN-EC", "SPN-MTG", "SPN-SVZ","SPN-PUT" ),
  count = c(17733 + 28483, 1868, 368, 1309, 29967)
)

cols <- c(
  "Other neurons" = "gray75",
  "SPN-EC" = "#832883",
  "SPN-MTG" = "#2d79aa",
  "SPN-PUT" = "#a9254c",
  "SPN-SVZ" = "#e4b722"
)

# http://localhost:41311/graphics/plot_zoom?width=374&height=176&scale=1
pie_data$group<-factor(pie_data$group,levels = c("Other neurons","SPN-EC","SPN-MTG","SPN-SVZ","SPN-PUT"))
cairo_pdf("/data/ADRD/brain_aging/exploration/plots/Phase1_DJArevisions/SPN_PIE_summary.pdf",width = 3.74,height = 1.76)
ggplot(pie_data, aes(x = "", y = count, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = cols) +
  theme_void() +
  theme(legend.position = "right")
dev.off()
