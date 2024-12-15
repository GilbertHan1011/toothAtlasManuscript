dpt <- zellkonverter::readH5AD("../../disk1/limb/important_processed_data/11.16_dpt.h5ad")
dptSeurat <- as.Seurat(dpt,,counts = "counts",data = "X")
head(dptSeurat$pred_dpt)
head(dptSeurat$lineage_chondro)
mesLineage = dptSeurat[,dptSeurat$lineage_mesenchyme]
mesLineage$sample_batch = paste0("Mes_",mesLineage$Sample)
mes_res <- scRNA_2_mat(mes = mesLineage,assay = "originalexp",slot = "count",pseudo_col = "pred_dpt",project_col = "sample_batch")

fibroLineage = dptSeurat[,dptSeurat$lineage_laFibro]
fibroLineage$sample_batch = paste0("Fibro_",fibroLineage$Sample)
fibro_res <- scRNA_2_mat(mes = fibroLineage,assay = "originalexp",slot = "count",pseudo_col = "pred_dpt",project_col = "sample_batch")

leprLineage = dptSeurat[,dptSeurat$lineage_lepr]
leprLineage$sample_batch = paste0("Lepr_",leprLineage$Sample)
lepr_res <- scRNA_2_mat(mes = leprLineage,assay = "originalexp",slot = "count",pseudo_col = "pred_dpt",project_col = "sample_batch")

chondroLineage = dptSeurat[,dptSeurat$lineage_chondro]
chondroLineage$sample_batch = paste0("Chondro_",chondroLineage$Sample)
chondro_res <- scRNA_2_mat(mes = chondroLineage,assay = "originalexp",slot = "count",pseudo_col = "pred_dpt",project_col = "sample_batch")

resList = list(mes_res,fibro_res,lepr_res,chondro_res)
filterList <- lapply(resList, function(x) x[["binned_means_filter"]])
commonGenes <- lapply(filterList,rownames)
common_genes <- Reduce(intersect, commonGenes)

combined_data <- Reduce(cbind, lapply(filterList, function(x) x[common_genes,]))
write.csv(combined_data,"processed_data//trajectory/20241214_input_osteogenicCount.csv")
