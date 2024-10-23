df <- readRDS("preprocess_data/Follile_Akira.Rds")
df$orig.ident <- "Follile_Akira"
saveRDS(df, "preprocess_data/Follile_Akira.Rds")

peri <- readRDS("preprocess_data/Periontal_Nagta.Rds")
DimPlot(peri)
peri$orig.ident <- "Periontal_Nagta"
saveRDS(peri, "preprocess_data/Follile_Akira.Rds")


#== Atlas_Jan-----------------------------
toothIncisorMolar <- readRDS("~/Desktop/disk1/tooth/10.19_nchuman/data/incisor_molar_data.rds")
incisorCombine <- readRDS("~/Desktop/disk1/tooth/saveData/incisorMolarConbine.Rds")

incisorCombine$orig.ident <- Idents(incisorCombine)

incisorCombine$orig.ident <- paste0("Atlas_Jan_",
                                    incisorCombine$orig.ident)



saveRDS(incisorCombine,"preprocess_data/Atlas_Jan_Mouse.Rds")

human <- readRDS("~/Desktop/disk1/tooth/10.19_nchuman/data/human_data.rds")
humanTooth <- readRDS("~/Desktop/disk1/tooth/saveData/humanTooth.Rds")

humanTooth$orig.ident <- paste0("Atlas_Jan_",humanTooth$X3)
saveRDS(humanTooth,"preprocess_data/Atlas_Jan_Human.Rds")


#== Epi_Chiba
Epi <- readRDS("~/Desktop/disk1/tooth/saveData/10.20_toothEpi.Rds")

Epi$orig.ident <- "Epi_Chiba"
saveRDS(Epi,"preprocess_data/Epi_Chiba.Rds")

#== Runx2
runx2 <- readRDS("~/Desktop/disk1/tooth/saveData/toothRunx2.Rds")
runx2$orig.ident <- "Runx2_Shuo"
saveRDS(runx2,"preprocess_data/Runx2_Shuo.Rds")


#== Liu
lh <- readRDS("~/Desktop/disk1/tooth/saveData/11.24_liuhuanTooth.Rds")
lh$orig.ident <- "Molar_Qian"
saveRDS(lh,"preprocess_data/Molar_Qian.Rds")
