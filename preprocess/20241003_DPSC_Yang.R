rm(list=ls())
source("script/utils/seurat_utils.R")
matrix <- read.csv("../202409_tooth_raw/DPSC_Yang/GSE227731/GSE227731_DPSC_PDLSC_scRNAseq_raw_count_matrix.csv",row.names = 1)
seurat <- CreateSeuratObject(matrix)
seurat <- qcFun(seurat,Species = "Human")
seurat <- runSeurat(seurat)

seurat$orig.ident <- "DPSC_Yang"
seurat@assays$RNA@layers$scale.data <- NULL
saveRDS(seurat,"preprocess_data/DPSC_Yang.Rds")

