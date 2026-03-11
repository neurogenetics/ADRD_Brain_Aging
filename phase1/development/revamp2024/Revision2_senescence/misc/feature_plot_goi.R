library(Seurat)

setwd("/data/ADRD/brain_aging/exploration/")
#asked by name by R2
goi<-c("CDKN1A","CDKN2A","SERPINE1")#asked for SASP-related cytokines, holding off on specific ones for now

#load object
obj = readRDS("./aging.pegasus.leiden_085.subclustered.rds")#final sce object converted to seurat
Idents(obj)<-"broad_celltype"

#feature plots
p1<-FeaturePlot(obj,features = goi)
#violin plots
p2<-VlnPlot(obj,features = goi,ncol = 1,pt.size = 0)

# http://localhost:35241/graphics/plot_zoom?width=574&height=507&scale=1
png("./revision2_senpaths/misc/GOI_FeaturePlots.png",width = 5.74,height = 5.07,units = "in",res = 600)
p1
dev.off()
# http://localhost:35241/graphics/plot_zoom?width=485&height=735&scale=1
png("./revision2_senpaths/misc/GOI_ViolinPlots.png",width = 4.85,height = 7.35,units = "in",res = 600)
p2
dev.off()

