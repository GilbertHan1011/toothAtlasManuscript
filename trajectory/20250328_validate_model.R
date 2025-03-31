setwd("../../../disk2/202409_tooth/")
library(TrajConserve)
savedir <- "process/trajectory/20250328_trajectory_model_test/"
hvg_gene <- read.csv("processed_data/framework/hvg/20250108_hvgEpi.csv",row.names = 1)%>% unlist()
epi_test <- readRDS("process/lzs_results/processed_data/integrated_data/20241112_epithelium.Rds")
epi_pseudo  <- read.csv("process/lzs_results/20241210_pseudotime.csv")
subsetName <- colnames(epi_test)[epi_test$C22_named!="mAM"]
subsetName2 <- epi_pseudo$X
namesSub <- intersect(subsetName,subsetName2)
epi_test <- epi_test[,namesSub]
rownames(epi_pseudo) <- epi_pseudo$X

epi_pseudo <- epi_pseudo[namesSub,]
epi_test@meta.data$pseudo <- epi_pseudo$lightGBM

subset_cell_num <- 8000
subset_cell <- sample(colnames(epi_test),size = subset_cell_num)
small_epi <- epi_test[hvg_gene[1:100],subset_cell]
saveRDS(small_epi,"trajConserve/data/small_example.Rds")
small_epi <- readRDS("trajConserve/data/small_example.Rds")
trajectory_data <- seurat_to_trajectory_array(
  seurat_obj = small_epi,
  assay = "originalexp",
  pseudo_col = "pseudo",
  project_col = "Project"
)
gene_idx <- 10  # Example gene index
model_test <- run_trajectory_model(trajectory_data$reshaped_data, gene_idx)

plot_results_brms(model_test)
library(trajConserve)
TrajConserve::run_multiple_models(trajectory_data$reshaped_data,
                                  gene_indices = 1:20,  # First 20 genes
                                  parallel = TRUE,
                                  n_cores = 20,save_metrics = T,save_metrics_file = paste0(savedir,"test1.hdf5"),
                                  save_plots = TRUE,save_plots_dir = savedir,
                                  save_models = TRUE,save_models_dir =savedir)

h5_test <- rhdf5::H5Fopen( paste0(savedir,"test1.hdf5"))
h5_test$"array_weights"


metric <-  TrajConserve::extract_hdf5_metric(paste0(savedir,"test1.hdf5"))
TrajConserve::plot_hdf5_heatmap(paste0(savedir,"test1.hdf5"))
TrajConserve::calculate_conservation(paste0(savedir,"test1.hdf5"))
