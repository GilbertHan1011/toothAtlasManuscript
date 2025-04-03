mes <- readRDS("processed_data/integrated_data/20250326_mesenchyme.Rds")
mesMeta <- mes@meta.data
write.csv(mesMeta,"processed_data/framework/annotation/20250401_mes_meta.csv")
epi <- readRDS("processed_data/")
mesMeta <- mes@meta.data
write.csv(mesMeta,"processed_data/framework/annotation/20250401_mes_meta.csv")
epi <- readRDS("processed_data/integrated_data/20250328_epithelium.Rds")
epiMeta <- epi@meta.data
write.csv(epiMeta,"processed_data/framework/annotation/20250402_epi_meta.csv")

library(Seurat)
DimPlot(epi)
mesUmap <- mes@reductions$X_umap@cell.embeddings
epiUmap <- epi@reductions$X_umap@cell.embeddings
write.csv(mesUmap,"processed_data/framework/embedding//20250402_mes_umap.csv")
write.csv(epiUmap,"processed_data/framework/embedding//20250402_epi_umap.csv")
