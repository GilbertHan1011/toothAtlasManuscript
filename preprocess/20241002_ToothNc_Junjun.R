ncJunjun <- readRDS("~/Desktop/disk1/tooth/saveData/Nctooth_Junjun_Merge.Rds")
ncJunjun@assays$RNA@scale.data <- matrix()
ncJunjun@assays$integrated <- NULL
idents <- ncJunjun$orig.ident
idents <- gsub("\\.\\/","",idents)
idents <- paste0("ToothNc_Junjun_",idents)
ncJunjun$orig.ident <- idents
saveRDS(ncJunjun,"preprocess_data/ToothNc_Junjun.Rds")
