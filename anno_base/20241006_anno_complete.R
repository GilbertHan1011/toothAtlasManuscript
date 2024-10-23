library(Seurat)
library(dplyr)
outputdir <- "process/annotation/first_round_base/"
fileName = list.files("process/annotation/first_round_base/anno/",pattern = "*.csv")%>% gsub(".csv","",.)
preprocessFile <- list.files("preprocess_data/",pattern = "*.Rds")%>% gsub(".Rds","",.)
setdiff(preprocessFile,fileName)
setdiff(fileName,preprocessFile)

#Deciduous_Li-------------------------------------
baseName <- "Deciduous_Li"
Deciduous_Li <- readRDS("preprocess_data/Deciduous_Li.Rds")
FeaturePlot(Deciduous_Li,features = c("SOX9","KRT14","RGS5"))
FeaturePlot(Deciduous_Li,features = c("VCAN","PRRX1"))
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))
Deciduous_Li$coarse_anno_1 <- "Mesenchyme"

label <- Deciduous_Li$coarse_anno_1 %>% as.data.frame()
colnames(label) <- "Coarse_Label_1"
write.csv(label,paste0(outputdir,"anno/","Deciduous_Li",".csv"))

#DPSC_Yang------------------------
baseName <- "DPSC_Yang"
seurat <- readRDS(paste0("preprocess_data/",baseName,".Rds"))
FeaturePlot(seurat,features = c("VCAN","POSTN"))
FeaturePlot(seurat,features = c("VCAN","IGFBP5","ORC6","CDK8","TUBB","HMGB1"))
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))
FeaturePlot(seurat,features = c("ASPN","S100A4"))
FeaturePlot(seurat,features = c("TUBB","HMGB1"))
DimPlot(seurat,group.by = "orig.ident")

seurat$coarse_anno_1 <- "Unkown"

label <- seurat$coarse_anno_1 %>% as.data.frame()
colnames(label) <- "Coarse_Label_1"
write.csv(label,paste0(outputdir,"anno/",baseName,".csv"))

#==Fibro_Xu-------------

baseName <- "Fibro_Xu"
seurat <- readRDS(paste0("preprocess_data/",baseName,".Rds"))
FeaturePlot(seurat,features = c("VCAN","POSTN","KRT14",
                                "RGS5","DMP1","CDH5",
                                "C1QA","NAPSA"))
FeaturePlot(seurat,features = c("CLU","S100A8"))
FeaturePlot(seurat,features = c("PTPRC"))
FeaturePlot(seurat,features = c("COL1A1"),label = T)
VlnPlot(seurat,"PTPRC")
VlnPlot(seurat,"POSTN")
VlnPlot(seurat,"SP7")
VlnPlot(seurat,"CDH5")
VlnPlot(seurat,"RGS5")
VlnPlot(seurat,"DCN")
VlnPlot(seurat,"MYOG")
VlnPlot(seurat,features = c("VCAN","POSTN","KRT14",
                                "RGS5","DMP1","CDH5",
                                "C1QA","NAPSA"),stack = T,flip = T)
DimPlot(seurat,label = T)
new.id <- c("Fibroblast", "Endothelium", "Immune", "Fibroblast",
            "Perivasular", "Fibroblast", "Fibroblast", "Fibroblast", "Immune",
            "Immune", "Mesenchyme", "Immune", "Endothelium", "Immune", "14", "15", "Immune")
seurat <- renameLabel(seurat,new.id,"coarse_anno_1") # first

#ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))
# FeaturePlot(seurat,features = c("ASPN","S100A4"))
# FeaturePlot(seurat,features = c("TUBB","HMGB1"))
# # DimPlot(seurat,group.by = "orig.ident")
#
# seurat$coarse_anno_1 <- "Unkown"

label <- seurat$coarse_anno_1 %>% as.data.frame()
colnames(label) <- "Coarse_Label_1"
write.csv(label,paste0(outputdir,"anno/",baseName,".csv"))

#==LYPD1_Chiba----------------
baseName <- "LYPD1_Chiba"
seurat <- readRDS(paste0("preprocess_data/",baseName,".Rds"))
markers <- c("Sox9","Vcan",  # Mesenchyme
             "Krt14","Pitx2",# Epithelium
             "C1qa","Napsa", #Immune
             "Cdh5", #Endo
             "Sox10", # neuron
             "Rgs5", # perivasular
             "Myog",
             "Hba-a2"
)
FeaturePlot(seurat,markers,label = T)
VlnPlot(seurat,markers,stack = T,flip = T)
ggsave(paste0(outputdir,"plot/",baseName,"_EDA_Vln.pdf"))
# FeaturePlot(seurat,markers,label = T)
# FeaturePlot(seurat,c("Myog","Hba-a2"),label = T)
# ggsave(paste0(outputdir,"plot/",baseName,"_EDA_FeatUmap.pdf"))

new.id <- c("Mesenchyme", "Mesenchyme", "Mesenchyme", "Mesenchyme", "Mesenchyme", "Mesenchyme", "Epithelium",
            "RBC", "Immune", "muscle", "Endothelium", "Mesenchyme"
)
seurat <- renameLabel(seurat,new.id,"coarse_anno_1") # first
p1 <- DimPlot(seurat,group.by = "coarse_anno_1")
p2 <- DimPlot(seurat,group.by = "orig.ident")
p1|p2

ggsave(paste0(outputdir,"plot/",baseName,"_validation_umap.pdf"))

label <- seurat$coarse_anno_1 %>% as.data.frame()
colnames(label) <- "Coarse_Label_1"
write.csv(label,paste0(outputdir,"anno/",baseName,".csv"))

