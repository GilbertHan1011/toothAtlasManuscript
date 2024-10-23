rm(list=ls())
source("script/utils/seurat_utils.R")
dirName <- "../202409_tooth_raw/Zebra_Gruenhagen/"
#files <- list.files("../202409_tooth_raw/CAGE_Chiba/",pattern = "Human*",full.names = T)
seuratMerge <- process_seurat_data(dirName)
seuratMerge@assays$RNA@layers$scale.data <- NULL
saveRDS(seuratMerge,"preprocess_data/Zebra_Gruenhagen.Rds")
#seuratMerge$orig.ident %>% unique
