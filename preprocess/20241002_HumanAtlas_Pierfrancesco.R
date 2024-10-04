tooth_iscience <- readRDS("~/Desktop/disk1/tooth/saveData/10.30_tooth_iscience.Rds")
DimPlot(tooth_iscience)
ident <- tooth_iscience$orig.ident %>% gsub("data\\/\\/","",.)
ident <- paste0("HumanAtlas_Pierfrancesco_",ident)
tooth_iscience$orig.ident <- ident
saveRDS(tooth_iscience,'preprocess_data/HumanAtlas_Pierfrancesco.Rds')


