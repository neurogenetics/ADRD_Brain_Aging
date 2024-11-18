library(Seurat)
library(SAHA)

phase1=readRDS("{path_to}/aging.pegasus.leiden_085.subclustered.rds")
phase2=readRDS("{path_to}/brain_aging_phase2/phase2_GEXonly.rds")

summary(factor(phase1$broad_celltype))
Idents(phase1)<-"broad_celltype"

summary(factor(phase2$leiden_MultiVI))
Idents(phase2)<-"leiden_MultiVI"


devtools::install_github("neurogenetics/SAHA")

SAHA_phase1_input=Seurat2SAHA(obj = phase1,output = "Both")
SAHA_phase2_input=Seurat2SAHA(obj = phase2,output = "Both")

saveRDS(SAHA_phase1_input,file = "{path_to}/brain_aging_general/SAHA_phase1_list.rds")
saveRDS(SAHA_phase2_input,file = "{path_to}/brain_aging_general/SAHA_phase2_list.rds")
