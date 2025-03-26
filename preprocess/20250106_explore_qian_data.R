library(Seurat)
qian <- readRDS("../202409_tooth/processed_data/preprocess_data/Molar_Qian.Rds")
source("script/utils/seurat_utils.R")
qian <- runSeurat(qian)


FeaturePlot(qian,c("Krt14","Ambn","Sp7","Ovol2"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))

DimPlot(qian)
FeaturePlot(qian,"Ambn")
FeaturePlot(qian2,"EGFP")
qian2 <- readRDS("../../disk1/tooth/saveData/11.24_liuhuanTooth.Rds")
qian3 <- readRDS("../../disk1/tooth/saveData/11.24_liuhuanTooth.Rds")
DimPlot(qian2)
DimPlot(qian3)


GFPCount <- Read10X("../../disk1/tooth/11.24_liuhuan/data/GFP/")
GFP <- CreateSeuratObject(GFPCount,min.cells = 3, min.features = 200,  project = "GFP")
qcfun <- function(x){
  x <- PercentageFeatureSet(x, "^mt-", col.name = "percent_mito")
  selected_count <- WhichCells(x, expression =( nCount_RNA > 1000 & percent_mito < 10 & nFeature_RNA > 500))
  x <- subset(x, cells = selected_count)
  return(x)
}
GFP <- qcfun(GFP)
dim(GFP)
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
GFP <- runSeurat(GFP)
DimPlot(GFP)

FeaturePlot(GFP,c("Mfap4","Col1a1","Sp7"))
FeaturePlot(GFP,c("Ambn",'Amelx',"Ovol2","Krt14","Bglap","Sp7","Runx2"))
"EGFP" %in% GFP@assays[["RNA"]]@layers$counts@Dimnames[[1]]
