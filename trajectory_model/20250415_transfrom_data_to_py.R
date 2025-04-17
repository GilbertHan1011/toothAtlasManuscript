epi <- readRDS("../202409_tooth/process/lzs_results/processed_data/integrated_data/20241112_epithelium.Rds")


pseudo <- read.csv("process/lzs_results/20241210_pseudotime.csv")
#pseudo <- read.csv("processed_data/trajectory/20241210_enamal_psudotime.csv")
#DimPlot(epi,group.by = "C22_named")
subsetName <- colnames(epi)[epi$C22_named!="mAM"]
subsetName2 <- pseudo$X
namesSub <- intersect(subsetName,subsetName2)
epi <- epi[,namesSub]
rownames(pseudo) <- pseudo$X

pseudo <- pseudo[namesSub,]
epi$pseudo <- pseudo$lightGBM

epiSce <- as.SingleCellExperiment(epi)
zellkonverter::writeH5AD(epiSce,"processed_data/integrated_data/20250414_epi_adata.h5ad")
mes <- readRDS("processed_data/integrated_data/20241106_mesenchyme.Rds")
pseudo <- read.csv("processed_data/trajectory/20241124_pseudotime_predicted.csv",row.names = 1)
mes$pseudotime <- pseudo$lightGBM

mesSce <- as.SingleCellExperiment(mes)
zellkonverter::writeH5AD(mesSce,"processed_data/integrated_data/20250415_mes_adata.h5ad")


