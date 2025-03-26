incisorMes <- readRDS("processed_data/annotation_harmonization/20250119_ncAtlas_incisorMes.Rds")
#incisorMes <- readRDS("processed_data/annotation_harmonization/20250119_ncAtlas_incisorMes.Rds")
nc_sub <- readRDS("../202409_tooth/processed_data/annotation_harmonization/20250119_ncjunjun.Rds")
mesMerge <- merge(nc_sub,c(incisorMes,rat_mes_2))
mesMerge <- JoinLayers(mesMerge)
mesMergeSc <- as.SingleCellExperiment(mesMerge,assay = "RNA")

zellkonverter::writeH5AD(mesMergeSc,"../202409_tooth/processed_data/annotation_harmonization/20250119_mes_anno.h5ad")
