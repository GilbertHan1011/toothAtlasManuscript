#in R
rm(list =ls())
outputdir <- "process/annotation/first_round_base/"
baseName <- "HumanCaries_Anunya" # You should change this.
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
baseName <- "HumanCaries_Anunya" # You should change this.
print(paste0("preprocess_data/",baseName,"_v4.Rds"))
seurat <- readRDS(paste0("preprocess_data/",baseName,"_v4.Rds"))
source("script/utils/seurat_utils.R")

markers <- c("OMD","ALPL",  # Mesenchyme
             "KRT14",# Epithelium
             "CD3D",'GNLY', # Immune
             "CDH5", # Endothelium
             "SOX10", # neuron
             "RGS5", # perivasular(fibroblasts)
             'MZB1'# Plasma
) 
FeaturePlot(seurat,markers,label = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))
VlnPlot(seurat,markers,stack = T,flip = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_Vln.pdf"))
dput(levels(seurat))
newID <- c('Mesenchyme','Mesenchyme','Endothelium', 'Mesenchyme','Immune', 'Mesenchyme','Immune', 'perivascular', 'Immune', 'Mesenchyme','neuron', 'neuron', 'Epithelium', 'Plasma')
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
