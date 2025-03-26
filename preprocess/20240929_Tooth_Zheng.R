library(Seurat)
source("script/utils/seurat_utils.R")
M1 <- Read10X("../202409_tooth_raw/Tooth_Zheng//M1/")
M2 <- Read10X("../202409_tooth_raw/Tooth_Zheng//M2/")
M1_seurat <- CreateSeuratObject(M1)
M2_seurat <- CreateSeuratObject(M2)
seurat_list <- list(M1_seurat,M2_seurat)
seurat_list <- lapply(seurat_list,qcFun)
seurat_list[[1]]$orig.ident <- "Tooth_Zheng_M1"
seurat_list[[2]]$orig.ident <- "Tooth_Zheng_M2"
seuratMerge <- merge(seurat_list[[1]],seurat_list[[2]])
seuratMerge <- runharmony(seuratMerge)

DimPlot(seuratMerge, group.by = "orig.ident")
markers <-
c("Sox9",
  "Msx2",
  "Vcan",
  "Kit",
  "Pclo",
  "Igfbp5",
  "Postn",
  "Fst",
  "Igfbp2",
  "Smoc2",
  "Sfrp2",
  "Igfbp3",
  "Smpd3",
  "Gsc",
  "Wnt10a",
  "Cdk1",
  "Pitx2",
  "Krt14",
  "Dsp",
  "C1qa",
  "Aif1",
  "Napsa",
  "Cdh5",
  "Rgs5",
  "Acp5"
)

FeaturePlot(seuratMerge,markers[1:9])
FeaturePlot(seuratMerge,markers[10:18])
FeaturePlot(seuratMerge,"Wnt10a")
saveRDS(seuratMerge,"preprocess_data/Tooth_Zheng.Rds")


rat_mes_2 <- seuratMerge[,colnames(rat_mes)]
DimPlot(rat_mes_2, group.by = "curate_anno",label = T)
ggsave("results/annotation/anno_harmonize/20250119_mes/20250119_rat_combine.pdf")

rat_mes_2$curate_anno <- rat_mes$label
