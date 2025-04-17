epi <- readRDS("../202409_tooth/process/lzs_results/processed_data/integrated_data/20241112_epithelium.Rds")
epi$C12_named

epi$trajMarker <- "Other"
epi$trajMarker[epi$C12_named%in%c("sAM", "PA", "mAM","PA and eAM", "eAM and sAM")] = "Ameloblast"
epi$trajMarker[epi$C12_named%in%c("DEP", "OEE lineage", "SR")] = "Stem"
Idents(epi) <- epi$trajMarker
epiMarker <- FindMarkers(epi,ident.1 = "Ameloblast",ident.2 = "Stem")
write.csv(epiMarker,"process/trajectory/20250415_GO_downstream/20250415_epimarker_seurat.csv")
