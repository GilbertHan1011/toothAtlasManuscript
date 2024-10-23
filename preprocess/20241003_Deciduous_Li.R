rm(list=ls())
source("script/utils/seurat_utils.R")
dirName <- "../202409_tooth_raw/Deciduous_Li/10X/"
seurat <- Read10X(dirName)
seurat <- CreateSeuratObject(seurat)
seurat <- qcFun(seurat,Species = "Human")
seurat <- runSeurat(seurat)
seurat@assays$RNA@layers$scale.data <- matrix()
seurat$orig.ident <- "Deciduous_Li"
saveRDS(seurat,"preprocess_data/Deciduous_Li.Rds")
