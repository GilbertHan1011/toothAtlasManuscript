source("script/utils/trajectory_model_util.R")
library(abind)
library(RColorBrewer)

epi <- readRDS("../202409_tooth/process/lzs_results/processed_data/integrated_data/20241112_epithelium.Rds")
pseudo <- read.csv("process/lzs_results/20250106_pseudotime.csv",row.names = 1)
pseudo <- pseudo[colnames(epi),]
epi$pseudo <- pseudo$lightGBM
DimPlot(epi,group.by = "C22_named")
epi <- epi[,!is.na(epi$pseudo)]
FeaturePlot(epi,"pseudo")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))

epi <- epi[,epi$C22_named != "OEE Anxa1"]
#DimPlot(epi_sub,group.by = "C22_named")
#pseudo <- read.csv("processed_data/trajectory/20241210_enamal_psudotime.csv")
# #DimPlot(epi,group.by = "C22_named")
# subsetName <- colnames(epi)[epi$C22_named!="mAM"]
# subsetName2 <- pseudo$X
# namesSub <- intersect(subsetName,subsetName2)
# epi <- epi[,namesSub]
# rownames(pseudo) <- pseudo$X
#
# pseudo <- pseudo[namesSub,]
# epi$pseudo <- pseudo$lightGBM

#FeaturePlot(epi,"Ambn")

res <- scRNA_2_mat(epi,assay = "originalexp",slot = "data",pseudo_col = "pseudo",project_col = "Project",n_bin = 100,ensure_tail = T)
saveRDS(res$reshaped_data,"processed_data/trajectory/20250106_epi_reshaped_v2.Rds")
epiMat <- res$reshaped_data

mes <- readRDS("processed_data/integrated_data/20241106_mesenchyme.Rds")
pseudo <- read.csv("processed_data/trajectory/20241124_pseudotime_predicted.csv",row.names = 1)
#mes <- mes[varGene,rownames(pseudo)]
mes <- mes[,rownames(pseudo)]
pseudo <- pseudo[colnames(mes),]
mes$pseudotime <- pseudo$lightGBM
#varGene <- read.csv("processed_data/framework/geneMeta/20241130_mes_vargene2_2000.csv",row.names = 1) %>% unlist()
resMes <- scRNA_2_mat(mes,assay = "originalexp",slot = "data",pseudo_col = "pseudotime",project_col = "Project",n_bin = 100,ensure_tail = T)
mesBin <- resMes$reshaped_data

osteoMat <- readRDS("processed_data/trajectory/20241223_input_osteogenicData_reshaped.sync-conflict-20241224-163620-POOXVVO.Rds")
osteoSample <- sample(dimnames(osteoMat)[[1]], 20)
osteoMatSample <- osteoMat[osteoSample,,]

geneName1 <- dimnames(osteoMat)[[3]]
geneName2 <- dimnames(epiMat)[[3]]
geneName3 <- dimnames(mesBin)[[3]]


common_genes <- Reduce(intersect, list(geneName1,geneName2,geneName3))

reshape1 <- osteoMatSample[,,common_genes]
reshape2 <- epiMat[,,common_genes]
reshape3 <- mesBin[,,common_genes]
combined_reshape <- abind(reshape1,reshape2, reshape3, along=1)
dimnames(combined_reshape)[[1]][21:30] <- paste0("Epi_",dimnames(combined_reshape)[[1]][21:30])
dimnames(combined_reshape)[[1]][31:39] <- paste0("Mes_",dimnames(combined_reshape)[[1]][31:39])

saveRDS(combined_reshape,"processed_data/trajectory/20250106_all_combined.Rds")
