#in R
rm(list =ls())
outputdir <- "process/annotation/first_round_base/"
baseName <- "MolarP12_Tomoko" # You should change this.
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
baseName <- "MolarP12_Tomoko" # You should change this.
print(paste0("preprocess_data/",baseName,"_v4.Rds"))
seurat <- readRDS(paste0("preprocess_data/",baseName,"_v4.Rds"))
source("script/utils/seurat_utils.R")
#annotation
DimPlot(seurat, reduction = "umap", label=TRUE)
all_markers <- Seurat::FindAllMarkers(seurat,
                                      only.pos = TRUE,
                                      logfc.threshold= 1,
                                      test.use = "wilcox",
)
top_markers <- all_markers %>% group_by(cluster) %>% top_n( n = 30, wt= avg_log2FC )
#fake annotation
markers <- c("Sox9","Vcan",  # Mesenchyme
             "Krt14","Pitx2",# Epithelium
             "C1qa","Napsa", # Immune
             "Cdh5", # Endo
             "Sox10", # neuron
             "Rgs5"# perivasular
) 
FeaturePlot(seurat,markers,label = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))
VlnPlot(seurat,markers,stack = T,flip = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_Vln.pdf"))
dput(levels(seurat))
newID <- c('Mesenchyme', 'Immune', 'neuron', 'perivascular', 'Endothelium', 'neuron', 'Immune','Mesenchyme', 'Immune','Immune','Mesenchyme', 'Endothelium', 'Immune','Immune','Immune','Immune')
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

