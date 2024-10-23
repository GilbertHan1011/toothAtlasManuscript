## Convert files to h5ad
rm(list=ls())
library(Seurat)
library(dplyr)
source("script/utils/seurat_utils.R")
annoBind <- read.csv("process/annotation/first_round_base/anno/Anno_summary.csv")
meta <- read.csv("data/metadata/metadata_latest.csv",row.names = 1)
projSelect <- meta$Project[meta$Core_datasets==TRUE]
projSelect <- unique(projSelect)
projSelect <- projSelect[!is.na(projSelect)]
projSelect[projSelect=="Atlas_Jan_Mouse"] = "Atlas_Jan_Mouse_renamed"
readAllData <- function(projName){
  seurat <- readRDS(paste0("processed_data/preprocess_data/",projName,".Rds"))
  annotation <- annoBind$Coarse_Label_1[annoBind$project == projName]
  seurat$coarse_anno_1 <- annotation
  return(seurat)
}
atlasJan <- readRDS("processed_data/preprocess_data/Atlas_Jan_Mouse_renamed.Rds")
seuratList <- lapply(projSelect,readAllData)
names(seuratList) <- projSelect
seuratList[["Atlas_Jan_Mouse"]] <- atlasJan
duplicated(rownames(atlasJan)) %>% sum
names <- make.unique(rownames(atlasJan))
atlasJan <- RenameGenesSeurat(atlasJan,newnames = names)

seuratList[["Atlas_Jan_Mouse"]] <- atlasJan




seuratList <- lapply(seuratList,function(x){
  try(x@assays$RNA@layers$scale.data <- NULL,silent = TRUE)
  return(x)
})

origId <- lapply(seuratList,function(x){
  return(unique(x$orig.ident))
})
origId <- origId %>% unlist()

sampleSelect <- rownames(meta)[meta$Core_datasets==TRUE] %>% unique
sampleSelect <- na.omit(sampleSelect)

setdiff(sampleSelect,origId)
seuratMerge <- merge(seuratList[[1]],seuratList[2:17], add.cell.ids = projSelect)
seuratMerge <- JoinLayers(seuratMerge)
seuratMergeSce <- as.SingleCellExperiment(seuratMerge,assay = c("RNA") )



zellkonverter::writeH5AD(seuratMergeSce,"process/pre-intergration/big_data/20241018_merge_all_step0.h5ad")
