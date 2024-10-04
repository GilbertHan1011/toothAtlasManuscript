rm(list=ls())
source("script/utils/seurat_utils.R")
matrix1 <- data.table::fread("../202409_tooth_raw/Premolar_Yin/GSE167251/GSM5100352_SC327_RSEC_MolsPerCell.csv.gz")
matrix1 <- t(matrix1)

colnames(matrix1) <- matrix1[1,]
matrix1 <- matrix1[2:nrow(matrix1),]

#meta <- readxl::read_xlsx("../202409_tooth_raw/Premolar_Yin/GSE167251/GSM5100352_SC327_RSEC_MolsPerCell_P.xlsx")
seurat <- CreateSeuratObject(matrix1)
matrix2 <- data.table::fread("../202409_tooth_raw/Premolar_Yin/GSE167251/GSM5100353_SC328_RSEC_MolsPerCell.csv.gz")
matrix2 <- t(matrix2)

colnames(matrix2) <- matrix2[1,]
matrix2 <- matrix2[2:nrow(matrix2),]

seurat <- CreateSeuratObject(matrix1,project = "Premolar_Yin_1")
seurat2 <- CreateSeuratObject(matrix2,project = "Premolar_Yin_2")
seuratList <- list(seurat,seurat2)
seuratList <- lapply(seuratList,qcFun,Species = "Human")
seuratMerge <- merge(seuratList[[1]],seuratList[[2]])
seuratMerge <- runharmony(seuratMerge)


saveRDS(seuratMerge,"preprocess_data/Premolar_Yin.Rds")
