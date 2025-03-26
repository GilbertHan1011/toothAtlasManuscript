source("script/utils/trajectory_model_util.R")


epi <- readRDS("../202409_tooth/process/lzs_results/processed_data/integrated_data/20241112_epithelium.Rds")
pseudo <- read.csv("process/lzs_results/20241210_pseudotime.csv")
pseudo <- read.csv("processed_data/trajectory/20241210_enamal_psudotime.csv")
#DimPlot(epi,group.by = "C22_named")
subsetName <- colnames(epi)[epi$C22_named!="mAM"]
subsetName2 <- pseudo$X
namesSub <- intersect(subsetName,subsetName2)
epi <- epi[,namesSub]
rownames(pseudo) <- pseudo$X

pseudo <- pseudo[namesSub,]
epi$pseudo <- pseudo$lightGBM

#FeaturePlot(epi,"Ambn")

res <- scRNA_2_mat(epi,assay = "originalexp",slot = "data",pseudo_col = "pseudo",project_col = "Sample",n_bin = 100,ensure_tail = T)

# Extract expression matrix
mesDf <- GetAssayData(epi,assay = "originalexp",slot = "data") %>% as.matrix()
pseudotime <- epi@meta.data[["pseudo"]]
batch <- epi@meta.data[["Project"]]

# Bin pseudotime
pseudotime_binned <- .bin_to_100(pseudotime,n_bins = 100)
metaDf <- data.frame(batch, pseudotime_binned)
metaDf$bin <- paste0(metaDf$batch, "_", metaDf$pseudotime_binned)

binned_means <- .calculate_bin_means_fast1(mesDf, metaDf$bin)
colnames(binned_means) <- unique(metaDf$bin)
# Filter genes and batches
geneNum <- round(0.1 * ncol(binned_means))
filteredGene <- rownames(binned_means)[rowSums(binned_means > 0) > geneNum]
head(colnames(binned_means))
#parts <- strsplit(colnames(binned_means), "_")
#prefixes <- sapply(parts, function(x) paste(x[1:2], collapse = "_"))
prefixes <- sapply(strsplit(colnames(binned_means), "_"),
                   function(x) paste(x[-length(x)], collapse = "_"))
numbers <- sapply(strsplit(colnames(binned_means), "_"),
                  function(x) paste(x[length(x)], collapse = "_")) %>% as.numeric()


bath_thred_real <- batch_thred * 100
batchName <- names(table(prefixes) > bath_thred_real)[table(prefixes) > bath_thred_real]
if (TRUE){
  remain_sample <- .examTail(metaDf,n_bin = 100,tail_width = 0.3,tail_num=0.02)
  batchName <- intersect(remain_sample,batchName)
}
binned_means_filter <- binned_means[filteredGene, prefixes %in% batchName]

# Reshape to 3D array


# Process final data
prefixes <- sapply(strsplit(colnames(binned_means_filter), "_"),
                   function(x) paste(x[-length(x)], collapse = "_"))
numbers <- sapply(strsplit(colnames(binned_means_filter), "_"),
                  function(x) paste(x[length(x)], collapse = "_")) %>% as.numeric()
reshaped_data <- .reshape_to_3d(binned_means_filter, prefixes, numbers,n_bins =100)
saveRDS(reshaped_data,"processed_data/trajectory/20241226_epi_reshaped.Rds")
reshaped_data <- readRDS("processed_data/trajectory/20241226_epi_reshaped.Rds")
