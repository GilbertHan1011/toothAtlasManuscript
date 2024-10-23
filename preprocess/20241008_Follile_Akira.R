library(Seurat)
source("script/utils/seurat_utils.R")
dpCount <- Read10X("../202409_tooth_raw/Follicle_Akira/data/")

dpSeurat <- CreateSeuratObject(dpCount,min.cells = 3, min.features = 500,project = "Follicle_Akira")
dpSeurat <- qcFun(dpSeurat)
dpSeurat <- runSeurat(dpSeurat)

dpSeurat$orig.ident <- "Follicle_Akira"
dpSeurat@assays$RNA@layers$scale.data <- NULL
saveRDS(dpSeurat,"preprocess_data/Follile_Akira.Rds")

markers <- c("Sox9","Vcan",  # Mesenchyme
             "Krt14","Pitx2",# Epithelium
             "C1qa","Napsa", #Immune
             "Cdh5", #Endo
             "Sox10", # neuron
             "Rgs5", # perivasular
             "Prrx1",
             'Hba-a2',# RBC
             'Myog'# muscle
)

FeaturePlot(dpSeurat,markers,label = T)
VlnPlot(dpSeurat,markers,stack = T,flip = T)

new.id <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
            "12")
new.id <- c("Mesenchyme", "Mesenchyme", "Epithelium", "Mesenchyme", "Epithelium", "Mesenchyme", "Immune", "Immune", "Mesenchyme", "Other", "RBC", "RBC",
            "Mesenchyme")
dpSeurat <- renameLabel(dpSeurat,new.id,"coarse_anno_1")

DimPlot(dpSeurat,group.by = "coarse_anno_1")
label <- dpSeurat$coarse_anno_1 %>% as.data.frame()
colnames(label) <- "Coarse_Label_1"

write.csv(label,"process/annotation/first_round_base/anno/Follile_Akira.csv")
#loadLabel <- read.csv(paste0(outputdir,"anno/",baseName,".csv"),row.names = 1)
# saveRDS(dpSeurat)
# hist(dpSeurat$nCount_RNA,breaks = 100)
# qcfun <- function(x){
#   x <- PercentageFeatureSet(x, "^mt-", col.name = "percent_mito")
#   selected_count <- WhichCells(x, expression =( nCount_RNA > 2000 & percent_mito < 15 & nFeature_RNA > 500))
#   x <- subset(x, cells = selected_count)
#   return(x)
# }
# dpSeurat <- qcfun(dpSeurat)
#
# runSeurat <- function(x){
#   x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#   all.genes <- rownames(x)
#   x <- ScaleData(x, features = all.genes)
#   x <- RunPCA(x, features = VariableFeatures(object = x))
#   x <- FindNeighbors(x, dims = 1:30)
#   x <- FindClusters(x, resolution = 0.5)
#   x <- RunUMAP(x, dims = 1:30)
# }
# dpSeurat <- runSeurat(dpSeurat)
# DimPlot(dpSeurat)
# FeaturePlot(dpSeurat,c("Sp7","Alpl","Mfap4","Aspn","Pthlh","Dmp1","Acan","Sox9"))
# saveRDS(dpSeurat,"../saveData//10.18_dental_Pthrp.Rds")
