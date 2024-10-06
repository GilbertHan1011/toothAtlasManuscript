#in R
rm(list =ls())
outputdir <- "process/annotation/first_round_base/"
baseName <- "Epi_Chiba" # You should change this.
print(paste0("preprocess_data/",baseName,".Rds"))
seurat <- readRDS(paste0("preprocess_data/",baseName,".Rds"))
source("script/utils/seurat_utils.R")
seurat[["RNA"]] <- as(seurat[["RNA"]], "Assay")
saveRDS(seurat,paste0("preprocess_data/",baseName,"_v4.Rds"))

#in Rstudio
rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
outputdir <- "process/annotation/first_round_base/"
baseName <- "Epi_Chiba" # You should change this.
print(paste0("preprocess_data/",baseName,"_v4.Rds"))
seurat <- readRDS(paste0("preprocess_data/",baseName,"_v4.Rds"))
source("script/utils/seurat_utils.R")
seurat <- subset(seurat, subset = nFeature_RNA > 1500 & percent_mito < 10)
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst",
                               nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:10)
seurat <- FindClusters(seurat, verbose = T, resolution = 0.5)
seurat <- RunUMAP(seurat,
                  reduction = "pca",
                  dims = 1:10,
                  verbose=TRUE )
markers <- c("Sox9",# Mesenchyme
             "Krt14",# Epithelium
             "C1qa", #Immune
             "Cdh5"#Endothelium
)
FeaturePlot(seurat,markers,label = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))
VlnPlot(seurat,markers,stack = T,flip = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_Vln.pdf"))
dput(levels(seurat))
newID <- c('Mesenchyme','Epithelium','Epithelium','Epithelium', 'Mesenchyme','Mesenchyme','Mesenchyme','Epithelium','Endothelium', 'Epithelium', 'Mesenchyme','Immune')
seurat <- renameLabel(seurat,newID,"coarse_anno_1")
p1 <- DimPlot(seurat,group.by = "coarse_anno_1")
p2 <- DimPlot(seurat,group.by = "orig.ident")
p1|p2
ggsave(paste0(outputdir,"plot/",baseName,"_validation_umap.pdf"))
label <- seurat$coarse_anno_1 %>% as.data.frame()
colnames(label) <- "Coarse_Label_1"

write.csv(label,paste0(outputdir,"anno/",baseName,".csv"))
loadLabel <- read.csv(paste0(outputdir,"anno/",baseName,".csv"),row.names = 1)
loadLabel %>% head()
