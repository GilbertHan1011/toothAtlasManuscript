rm(list=ls())
source("script/utils/seurat_utils.R")
dirName <- "../202409_tooth_raw/LYPD1_Chiba/LYPD1_Chiba/"
#files <- list.files("../202409_tooth_raw/CAGE_Chiba/",pattern = "Human*",full.names = T)
seurat <- Read10X(dirName)
seurat <- CreateSeuratObject(seurat)
seurat <- qcFun(seurat)
seurat <- runSeurat(seurat)
saveRDS(seurat,"preprocess_data/LYPD1_Chiba.Rds")
#seuratMerge$orig.ident %>% unique
