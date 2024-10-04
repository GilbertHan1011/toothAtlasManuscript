rm(list=ls())
source("script/utils/seurat_utils.R")
dirName  = "~/Desktop/disk1/tooth/12.7_p12Pulp/data/"
seuratMerge <- process_seurat_data(dirName)
saveRDS(seuratMerge,"preprocess_data/MolarP12_Tomoko.Rds")
