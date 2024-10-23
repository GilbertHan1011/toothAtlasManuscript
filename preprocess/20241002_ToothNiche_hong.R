rm(list=ls())
source("script/utils/seurat_utils.R")
dirName  = "~/Desktop/disk1/tooth/1.9_tooth_cellreport/"
seuratMerge <- process_seurat_data(dirName)
seuratMerge@assays$RNA@layers$scale.data <- NULL
saveRDS(seuratMerge,"preprocess_data/ToothNiche_Hong.Rds")
