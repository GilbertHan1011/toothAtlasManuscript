library(Seurat)
library(SeuratData)
library(SeuratDisk)

mes <- readRDS("../../processed_data/integrated_data/20241106_mesenchyme.Rds")
SaveH5Seurat(mes, filename = "../../processed_data/20241118_mes.h5Seurat")
Convert("../../processed_data/20241118_mes.h5Seurat", dest = "h5ad")
