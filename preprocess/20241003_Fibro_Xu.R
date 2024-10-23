rm(list=ls())
source("script/utils/seurat_utils.R")
dirName <- "../202409_tooth_raw/Fibro_Xu/10X/"
#files <- list.files("../202409_tooth_raw/CAGE_Chiba/",pattern = "Human*",full.names = T)
mat <- Read10X(dirName)

name1 <- read.table("../202409_tooth_raw/Fibro_Xu/GSM7029372_HPLSC.cellname.txt.gz",header = T)

matM <- data.table::fread("../202409_tooth_raw/Fibro_Xu/GSM7029373_Tooth-M.counts.tsv.gz",header = T) %>% as.data.frame()
rownames(matM) <- matM$gene
matM <- matM[,2:ncol(matM)]

matF <- data.table::fread("../202409_tooth_raw/Fibro_Xu/GSM7029372_HPLSC.counts.tsv.gz",header = T) %>% as.data.frame()
rownames(matF) <- matF$gene
matF <- matF[,2:ncol(matF)]

seurat1 <- CreateSeuratObject(matM,min.cells = 3,min.features = 300,project = "Fibro_Xu_M")
seurat2 <- CreateSeuratObject(matF,min.cells = 3,min.features = 300,project = "Fibro_Xu_F")
seuratList <- list(seurat1,seurat2)
seuratList <- lapply(seuratList, qcFun,Species = "Human")
seuratMerge <- merge(seuratList[[1]],seuratList[[2]])
seuratMerge <- runharmony(seuratMerge)
seuratMerge@assays$RNA@layers$scale.data <- NULL
saveRDS(seuratMerge,"preprocess_data/Fibro_Xu.Rds")
#seuratMerge$orig.ident %>% unique
