rm(list=ls())
library(tidyverse)
source("script/utils/seurat_utils.R")

filesPath <- list.files("preprocess_data/",full.names = T)
fileName <- list.files("preprocess_data/") %>% gsub(".Rds","",.)

#sce <-readH5AD("process/pre-intergration/big_data/20241007_mergeall_filter_gene_step1.h5ad")

runDblFinder <- function(name){
  path <- paste0("processed_data/preprocess_data/",name,".Rds")
  seurat <- readRDS(path)
  try(seurat <- JoinLayers(seurat),silent = TRUE)
  try(seurat@assays$RNA@layers$scale.data <- NULL,silent = TRUE)
  seurat <- dbFinderFun(seurat)
  label <- seurat$scDblFinder_class %>% as.data.frame()
  colnames(label) <- "scDblFinder_class"
  write.csv(label,paste0("process/dblFinder/",name,".csv"))
}
lapply(fileName,runDblFinder)
#runDblFinder("Deciduous_Li")


#lapply(fileName[20:25],runDblFinder)
#
# test <- readRDS("preprocess_data//Deciduous_Li.Rds")
# test
#runDblFinder("Peridontal_Nagata")
