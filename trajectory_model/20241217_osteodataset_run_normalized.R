library(abind)

mesLineage = dptSeurat[,dptSeurat$lineage_mesenchyme]
mesLineage$sample_batch = paste0("Mes_",mesLineage$Sample)
mes_res <- scRNA_2_mat(mes = mesLineage,assay = "originalexp",slot = "data",pseudo_col = "pred_dpt",project_col = "sample_batch")

fibroLineage = dptSeurat[,dptSeurat$lineage_laFibro]
fibroLineage$sample_batch = paste0("Fibro_",fibroLineage$Sample)
fibro_res <- scRNA_2_mat(mes = fibroLineage,assay = "originalexp",slot = "data",pseudo_col = "pred_dpt",project_col = "sample_batch")

leprLineage = dptSeurat[,dptSeurat$lineage_lepr]
leprLineage$sample_batch = paste0("Lepr_",leprLineage$Sample)
lepr_res <- scRNA_2_mat(mes = leprLineage,assay = "originalexp",slot = "data",pseudo_col = "pred_dpt",project_col = "sample_batch")

chondroLineage = dptSeurat[,dptSeurat$lineage_chondro]
chondroLineage$sample_batch = paste0("Chondro_",chondroLineage$Sample)
chondro_res <- scRNA_2_mat(mes = chondroLineage,assay = "originalexp",slot = "data",pseudo_col = "pred_dpt",project_col = "sample_batch")

resList = list(mes_res,fibro_res,lepr_res,chondro_res)
filterList <- lapply(resList, function(x) x[["binned_means_filter"]])
commonGenes <- lapply(filterList,rownames)
common_genes <- Reduce(intersect, commonGenes)

#reshape1 <- mes_res$reshaped_data
reshape1 <- mes_res$reshaped_data[,,common_genes]
reshape2 <- fibro_res$reshaped_data[,,common_genes]
reshape3 <- lepr_res$reshaped_data[,,common_genes]
reshape4 <- chondro_res$reshaped_data[,,common_genes]


combined_reshape <- abind(reshape1,reshape2, reshape3,reshape4, along=1)
saveRDS(combined_reshape,"processed_data//trajectory/20241223_input_osteogenicData_reshaped.Rds")


combined_data <- Reduce(cbind, lapply(filterList, function(x) x[common_genes,]))
write.csv(combined_data,"processed_data//trajectory/20241214_input_osteogenicData.csv")

preparedData <- prepare_data_for_gam(resa[,,1616])

# Fit model
fit_test_normal <- bayesian_gam_regression_nb_shape(
  preparedData$x,
  preparedData$y,
  preparedData$array_idx,
  n_knots = 5
)

testDataMat_normal <-  mes_res$reshaped_data

preparedDataInput <- prepare_data_for_gam(inputDataMat[,,1616])
time_start = Sys.time()
fit_test_normal_speed <- bayesian_gam_regression_nb_shape(
  preparedDataInput$x,
  preparedDataInput$y,
  preparedDataInput$array_idx,
  n_knots = 5
)
time_end = Sys.time()
# Fit model
fit_test_normal <- bayesian_gam_regression_nb_shape(
  preparedData$x,
  preparedData$y,
  preparedData$array_idx,
  n_knots = 5
)

pdf("results/trajectory/20241210_conserved_model_validation/20241217_normalized.pdf",width = 10,height = 4)
plot_results_brms(fit_test_normal)
dev.off()
hmData = testDataMat_normal[,,1616]
Heatmap(hmData,cluster_rows = F,cluster_columns = F)
hmData <- t(scale(t(hmData)))
Heatmap(hmData,cluster_rows = F,cluster_columns = F)



preparedData <- prepare_data_for_gam(testDataMat_normal_filter[,,1616])
hmData2 = testDataMat_normal_filter[,,1616]
Heatmap(hmData2,cluster_rows = F,cluster_columns = F)
hmData2 <- t(scale(t(hmData2)))
Heatmap(hmData2,cluster_rows = F,cluster_columns = F)
testDataMat_normal_filter <-  mes_res$reshaped_data
fit_test_normal_filter <- bayesian_gam_regression_nb_shape(
  preparedData$x,
  preparedData$y,
  preparedData$array_idx,
  n_knots = 5
)




pdf("results/trajectory/20241210_conserved_model_validation/20241217_normalized_filter.pdf",width = 10,height = 4)
plot_results_brms(fit_test_normal_filter)
dev.off()

