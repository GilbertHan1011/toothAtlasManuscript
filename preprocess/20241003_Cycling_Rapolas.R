rm(list=ls())
source("script/utils/seurat_utils.R")
dirName <- "../202409_tooth_raw/Cycling_Rapolas/10X/"
#files <- list.files("../202409_tooth_raw/CAGE_Chiba/",pattern = "Human*",full.names = T)
seuratMerge <- process_seurat_data(dirName)
seuratMerge@assays$RNA@layers$scale.data <- NULL
saveRDS(seuratMerge,"preprocess_data/Cycling_Rapolas.Rds")
#seuratMerge$orig.ident %>% unique
