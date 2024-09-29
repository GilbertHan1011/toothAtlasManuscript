rm(list=ls())
source("script/utils/seurat_utils.R")
dirName  = "../202409_tooth_raw/Incisor_Zhang/"
seuratMerge <- process_seurat_data(dirName)
saveRDS(seuratMerge,"preprocess_data/Incisor_Zhang.Rds")
