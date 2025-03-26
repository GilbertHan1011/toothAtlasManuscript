mesAtlas <- readRDS("../202409_tooth_raw/nc_atlas/mesenchyme_data.rds")
label <- as.character(mesAtlas$clusters)
mesCount <- mesAtlas$pagoda.object$varinfo$mat
mesCount <- data.table::fread("../202409_tooth_raw/nc_atlas/mesenchymal (1).txt")
mesCount <- as.matrix(mesCount)  # Convert to matrix
rownames(mesCount) <- mesCount[,1]  # Set rownames
mesCount <- mesCount[,-1]
mesSeurat <- CreateSeuratObject(mesCount)
mesKeep <- mesAtlas$emb %>%  rownames
mesSeurat <- mesSeurat[,mesKeep]
mesSeurat$label <- label


incisor <-  data.table::fread("../202409_tooth_raw/nc_atlas/counts_SS2_mouse_incisor.txt")
incisor <- as.matrix(incisor)  # Convert to matrix
rownames(incisor) <- incisor[,1]  # Set rownames
incisor <- incisor[,-1]
incisorSeurat <- CreateSeuratObject(incisor)
annoLabel <- data.table::fread("../202409_tooth_raw/nc_atlas/annotation_mouse_incisor.txt",header = F)
annoLabel <- as.data.frame(annoLabel)
rownames(annoLabel) <- annoLabel$V1
annoLabel <- annoLabel[colnames(incisor),]
incisorSeurat$label <- annoLabel$V2
incisorMes <- incisorSeurat[,incisorSeurat$label%in%c("Distal pulp", "Apical pulp" ,
                                                     "Maturing pulp","Pre-odontoblasts","Dental follicle 1",
                                                     "Dental follicle 2","Alveolar osteo.")]
incisorMes <- runSeurat(incisorMes)
DimPlot(incisorMes,group.by = "label",label = T)
ggsave("results/annotation/anno_harmonize/20250119_mes/20250119_incisor_combine.pdf")
saveRDS(incisorMes,"processed_data/annotation_harmonization/20250119_ncAtlas_incisorMes.Rds")


