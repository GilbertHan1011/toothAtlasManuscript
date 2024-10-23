rm(list=ls())
source("script/utils/seurat_utils.R")
jci <- data.table::fread("~/Desktop/disk1/tooth/11.2_jci_tooth/data/GSE142200_E14_M1_Single_cell_matrix.txt.gz")
rownames(jci) <- jci$V1
jci_sub <- jci[,2:ncol(jci)]
rownames(jci_sub) <- jci$V1
seurat <- CreateSeuratObject(jci_sub)
seurat <- qcFun(seurat)
seurat <- runSeurat(seurat)
DimPlot(seurat)
seurat$orig.ident <- "MolarE145_Outie"
saveRDS(seurat,"preprocess_data/MolarE145_Outie.Rds")
seurat@assays$RNA@layers$scale.data <- matrix()

