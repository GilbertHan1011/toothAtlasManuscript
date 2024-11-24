mes <- readRDS("../../disk1/zjl_sc/data/5.5_mesenchymal.Rds")
FeaturePlot(mes,"Notch2")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))


chai <- readRDS("../../disk1/limb/data/cranioAtlas/Mandible2020_Chai_early.Rds")
FeaturePlot(chai,"Notch2")
wt <- readRDS("../../disk1/limb/important_processed_data/5.4_wtintegrate_full_seurat.Rds")

FeaturePlot(wt,"Notch2")
chai <- wt[,wt$Project=="Mandible2020_Chai"]
FeaturePlot(chai,"Notch2")



chai2 <- readRDS("../../disk1/limb/data/cranioAtlas/Mandible2020_Chai_early.Rds")

source("../../disk1/limb/function/seurat_utils.R")

chai2 <- runSeurat(chai2)
