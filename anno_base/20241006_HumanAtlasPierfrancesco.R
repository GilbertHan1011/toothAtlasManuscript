#in Rstudio
rm(list =ls())
outputdir <- "process/annotation/first_round_base/"
baseName <- "HumanAtlas_Pierfrancesco" # You should change this.
print(paste0("preprocess_data/",baseName,".Rds"))
seurat <- readRDS(paste0("preprocess_data/",baseName,".Rds"))
source("script/utils/seurat_utils.R")


markers <- c('SOX9', "VCAN",  # Mesenchyme
             'KRT14',"PITX2",# Epithelium
             "C1QA","NAPSA",# Immune
             "CDH5", # Endo
             "SOX10", # neuron
             "RGS5"# perivasular
) 
FeaturePlot(seurat,markers,label = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))
VlnPlot(seurat,markers,stack = T,flip = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_Vln.pdf"))
dput(levels(seurat))
newID <- c('Mesenchyme', 'Endothelium', 'Mesenchyme', 'Endothelium', 'perivascular', 'Epithelium', 'neuron', 'perivascular', 'neuron', 'Epithelium', 'Endothelium','Immune','Mesenchyme')
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
