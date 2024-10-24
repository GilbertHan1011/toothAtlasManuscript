# Droplet discovery

## Introduction
For droplet discovery, we recommended you to read these wonderful tutorials:

[Quality Control](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html);

[Droplet detection](https://bioconductor.org/books/3.19/OSCA.advanced/doublet-detection.html);

Here, we adapted scDblFinder for droplet discovery.

## Utils that used in this step
The utils can be found in [github](https://github.com/GilbertHan1011/toothAtlasManuscript/blob/main/script/utils/seurat_utils.R).
```R
#== dbl funciton-------------
dbFinderFun <- function(x){
  require(scDblFinder)
  require(BiocParallel)
  print(unique(x$orig.ident))
  DefaultAssay(x) <- "RNA"
  sce <- as.SingleCellExperiment(x)
  sce <- scDblFinder(sce, samples="orig.ident", BPPARAM=MulticoreParam(10))
  droplet_class = sce$scDblFinder.class
  x$scDblFinder_class <- droplet_class
  return(x)
}
```



## Process Scripts
The utils can be found in [github](https://github.com/GilbertHan1011/toothAtlasManuscript/blob/master/anno_base/20241008_droplet.R).
```R
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
```

## Why this step is necessary
The most bad consequence of remaining droplets is that it will affect the clustering result.
In the follwing figure, we can see that there are several clusters that mainly composed by droplets, which are most likely to be artifacts.
![image](../img/scdbl.png)