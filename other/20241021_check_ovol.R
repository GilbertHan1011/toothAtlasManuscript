library(Seurat)
junjun <- readRDS("processed_data/preprocess_data/ToothNc_Junjun.Rds")
FeaturePlot(junjun,c("Ovol2","Sp7"))

hong <- readRDS("processed_data/preprocess_data/ToothNiche_Hong.Rds")
FeaturePlot(hong,c("Ovol2","Sp7"))
