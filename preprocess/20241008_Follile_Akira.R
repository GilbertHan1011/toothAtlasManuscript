library(Seurat)
source("script/utils/seurat_utils.R")
dpCount <- Read10X("../202409_tooth_raw/Follicle_Akira/data/")

dpSeurat <- CreateSeuratObject(dpCount,min.cells = 3, min.features = 500,project = "Follicle_Akira")

hist(dpSeurat$nCount_RNA,breaks = 100)
qcfun <- function(x){
  x <- PercentageFeatureSet(x, "^mt-", col.name = "percent_mito")
  selected_count <- WhichCells(x, expression =( nCount_RNA > 2000 & percent_mito < 15 & nFeature_RNA > 500))
  x <- subset(x, cells = selected_count)
  return(x)
}
dpSeurat <- qcfun(dpSeurat)

runSeurat <- function(x){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  x <- FindNeighbors(x, dims = 1:30)
  x <- FindClusters(x, resolution = 0.5)
  x <- RunUMAP(x, dims = 1:30)
}
dpSeurat <- runSeurat(dpSeurat)
DimPlot(dpSeurat)
FeaturePlot(dpSeurat,c("Sp7","Alpl","Mfap4","Aspn","Pthlh","Dmp1","Acan","Sox9"))
saveRDS(dpSeurat,"../saveData//10.18_dental_Pthrp.Rds")
