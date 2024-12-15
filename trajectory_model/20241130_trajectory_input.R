
source("script/utils/trajectory_model_util.R")
mes <- readRDS("processed_data/integrated_data/20241106_mesenchyme.Rds")
pseudo <- read.csv("processed_data/trajectory/20241124_pseudotime_predicted.csv",row.names = 1)
mes <- mes[varGene,rownames(pseudo)]
mes$pseudotime <- pseudo$lightGBM
varGene <- read.csv("processed_data/framework/geneMeta/20241130_mes_vargene2_2000.csv",row.names = 1) %>% unlist()

countAssay <-  as.matrix(mes@assays$originalexp@counts)
#countAssay2 <- as.matrix(mes@assays$originalexp@counts)/mes$size_factors
countAssay2 <- sweep(countAssay, 2, mes$size_factors, FUN = "/")
# countAssay3 <- sweep(mes@assays$originalexp@counts, 2, mes$size_factors, FUN = "/")
countAssay2 <- as(countAssay2,"dgCMatrix")
mes@assays$originalexp@counts <- countAssay2 # divide by size factor

varGeneMat <- scRNA_2_mat(mes,assay = "originalexp",slot = "count",pseudo_col = "pseudotime",project_col = "Project")
genes <- varGeneMat$binned_means_filter %>% rownames()
inputDataMat <- varGeneMat$reshaped_data
#genes2 <- binned_means_filter %>% rownames()

preparedData <- prepare_data_for_gam(inputDataMat[,,1])
fit_count_1_Sgk3 <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)

geneData3_scale <- t(scale(t(inputDataMat[,,1])))
Heatmap(geneData3_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
plot_results_brms(fit_count_1_Sgk3)

preparedData <- prepare_data_for_gam(inputDataMat[,,2])
fit_count_1_Cpa6 <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
plot_results_brms(fit_count_1_Cpa6)

geneData3_scale <- t(scale(t(inputDataMat[,,2])))
Heatmap(geneData3_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)


preparedData <- prepare_data_for_gam(inputDataMat[,,3])
fit_count_1_Prex2 <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
plot_results_brms(fit_count_1_Prex2)
geneData3_scale <- t(scale(t(inputDataMat[,,3])))
Heatmap(geneData3_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

preparedData <- prepare_data_for_gam(inputDataMat[,,434])
fit_count_1_Bglap <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
plot_results_brms(fit_count_1_Bglap)
geneData3_scale <- t(scale(t(inputDataMat[,,434])))
Heatmap(geneData3_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

preparedData <- prepare_data_for_gam(inputDataMat[,,431])
fit_count_1_Bglap2 <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
plot_results_brms(fit_count_1_Bglap2)
geneData3_scale <- t(scale(t(inputDataMat[,,431])))
Heatmap(geneData3_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
Heatmap(inputDataMat[,,431],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)


preparedData <- prepare_data_for_gam(inputDataMat[,,3])
fit_count_2_Prex2 <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
plot_results_brms(fit_count_2_Prex2)
geneData3_scale <- t(scale(t(inputDataMat[,,3])))
Heatmap(inputDataMat[,,3],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)


Heatmap(geneData3_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

varGeneMat2 <- scRNA_2_mat(mes,assay = "originalexp",slot = "data",pseudo_col = "pseudotime",project_col = "Project",batch_thred = 0.5)
inputDataMat <- varGeneMat2$reshaped_data

Heatmap(inputDataMat[,,431],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

preparedData <- prepare_data_for_gam(inputDataMat[,,431])
fit_count_2_Bglap2 <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
plot_results_brms(fit_count_2_Bglap2)
geneData3_scale <- t(scale(t(inputDataMat[,,431])))
Heatmap(geneData3_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
Heatmap(inputDataMat[,,431],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)


plot_results_brms(fit_count_3_Bglap2)


saveRDS(varGeneMat2$reshaped_data,"processed_data/trajectory/20241130_data_varmat.Rds")


pred_data <- data.frame(
  x = x_seq,
  # Use the first array for predictions
  array = factor(rep("Runx2_Shuo", 100))
)
Heatmap(inputDataMat[,,3],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
predictions_3 <- predict(fit_count_1_Prex2$fit, newdata = pred_data)
estimate <- predictions[,1]

plot(preparedData$x,preparedData$y)

preparedData2 <- prepare_data_for_gam(inputDataMat[,,385])
plot(preparedData2$x,preparedData2$y)


varGeneMat2 <- scRNA_2_mat(mes,assay = "originalexp",slot = "data",pseudo_col = "pseudotime",project_col = "Project")
inputDataMat <- varGeneMat2$reshaped_data

preparedData2 <- prepare_data_for_gam(inputDataMat[,,434])
plot(preparedData2$x,preparedData2$y)



preparedData <- prepare_data_for_gam(inputDataMat[,,3])
fit_count_2_Prex2 <- bayesian_gam_regression_fast(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)

