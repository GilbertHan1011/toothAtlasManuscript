## Convert files to h5ad
rm(list=ls())
library(Seurat)
library(dplyr)

annoBind <- read.csv("process/annotation/first_round_base/anno/Anno_summary.csv")
meta <- read.csv("data/metadata/metadata.csv",row.names = 1)
projSelect <- meta$Project[meta$Species != "Homo Sapiens"]
proj <- annoBind$project %>% unique()
projSelect[projSelect=="Atlas_Jan"] = "Atlas_Jan_Mouse"
readAllData <- function(projName){
  seurat <- readRDS(paste0("preprocess_data/",projName,".Rds"))
  annotation <- annoBind$Coarse_Label_1[annoBind$project == projName]
  seurat$coarse_anno_1 <- annotation
  return(seurat)
}
seuratList <- lapply(proj,readAllData)
names(seuratList) <- proj
seuratList <- seuratList[unique(projSelect)]

seuratList <- lapply(seuratList,function(x){
  try(x@assays$RNA@layers$scale.data <- NULL,silent = TRUE)
  return(x)
})
seuratMerge <- merge(seuratList[[1]],seuratList[2:17])
seuratMerge <- JoinLayers(seuratMerge)
seuratMergeSce <- as.SingleCellExperiment(seuratMerge,assay = c("RNA") )

zellkonverter::writeH5AD(seuratMergeSce,"process/pre-intergration/big_data/20241007_merge_all_step0.h5ad")
