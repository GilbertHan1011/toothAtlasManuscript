#in Rstudio
rm(list =ls())
outputdir <- "process/annotation/first_round_base/"
baseName <- "Atlas_Jan_Human" # You should change this.
print(paste0("preprocess_data/",baseName,".Rds"))
seurat <- readRDS(paste0("preprocess_data/",baseName,".Rds"))
source("script/utils/seurat_utils.R")

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", 
                               nfeatures = 2000)
seurat <- ScaleData(seurat, features = rownames(seurat))
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
seurat <- RunUMAP(seurat, dims = 1:30, random.seed = 0)
seurat <- FindClusters(seurat, verbose = T, resolution = 0.5,random.seed = 0)

DimPlot(seurat, reduction = "umap", label=TRUE)
markers <- c('SOX9', "VCAN",  # Mesenchyme
             'KRT14',"PITX2",# Epithelium
             "C1QA",# Immune
             "CDH5", # Endo
             "SOX10", # neuron
             "RGS5", # perivasular
             'HBA1'# RBC
) 
FeaturePlot(seurat,markers,label = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))
VlnPlot(seurat,markers,stack = T,flip = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_Vln.pdf"))
dput(levels(seurat))
newID <- c('Mesenchyme', 'Mesenchyme', 'Endothelium', 'perivascular', 'Mesenchyme', 'Epithelium', 'Mesenchyme', 'neuron', 'Endothelium', 'Endothelium', 'perivascular', 'RBC', 'Immune', 'Immune', 'neuron', 'Mesenchyme', 'Mesenchme', 'Immune', 'Mesenchyme')
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
