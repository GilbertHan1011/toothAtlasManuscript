rm(list=ls())
nc2 <- readRDS("../../disk1/tooth/saveData/Nctooth_Junjun_Merge.Rds")
nc_seurat <- readRDS("processed_data/preprocess_data/ToothNc_Junjun.Rds")
FeaturePlot(nc_seurat,c("Prrx1","Krt14","Epcam"))
FeaturePlot(nc_seurat,c("Twist1","Acan"))
DimPlot(nc_seurat,label = T)
selected_cluster <- edit(nc_seurat$seurat_clusters %>% unique)
nc_sub <- nc_seurat[,nc_seurat$seurat_clusters%in%selected_cluster]
DimPlot(nc_sub)
## E13.5----------------
nc_E13.5 <- nc_sub[,nc_sub$orig.ident=="ToothNc_Junjun_E13.5"]
nc_E13.5 <- runSeurat(nc_E13.5)
FeaturePlot(nc_E13.5,c("Prrx1","Krt14","Epcam"))
FeaturePlot(nc_E13.5,c("Tfap2b","Msx1","Acta2","Ibsp","Col2a1","Cdh5","Top2a","Fibin","Akap12"))
new_id <- c("Cycling cells", "Dermal fibroblast", "Dental mesenchyme", "Dermal fibroblast",
            "Dermal fibroblast", "Dermal fibroblast",
            "Dermal fibroblast", "Dental mesenchyme", "Dental mesenchyme", "Dental mesenchyme",
            "Others", "Others", "Chondrogenic cells", "Others", "Osteogenic cells", "Dermal fibroblast",
            "Cycling cells")
names(new_id) <- levels(nc_E13.5)
nc_E13.5 <- RenameIdents(nc_E13.5,new_id)

## E14.5----------------
nc_E14.5 <- nc_sub[,nc_sub$orig.ident=="ToothNc_Junjun_E14.5"]
nc_E14.5 <- runSeurat(nc_E14.5)
FeaturePlot(nc_E14.5,c("Prrx1","Krt14","Epcam"))
FeaturePlot(nc_E14.5, c("Tfap2b","Top2a","Fibin",
                        "Igf1","Zfhx3","Scx",
                        "Alpl","Ibsp"),label = T)
FeaturePlot(nc_E14.5, c("Krt14"),label = T)
Idents(nc_E14.5) <- nc_E14.5$seurat_clusters
FeaturePlot(nc_E14.5, c("Postn","Hmgb2","Igfbp3","Mfap4","Maf","Fos"),label = T)
#FeaturePlot(nc_E14.5, c("Postn","Hmgb2","Igfbp3","Mfap4","Maf","Fos"),label = T)
dput(levels(nc_E14.5))
new.id <- c("Dermal fibroblast", "Cycling cell", "Dental mesenchyme", "Dental mesenchyme",
            "Dermal fibroblast", "Dermal fibroblast", "Dermal fibroblast",
            "Dermal fibroblast", "Dermal fibroblast", "Dental mesenchyme", "Tenogenic cells", "Others",
            "Dermal fibroblast", "Dermal fibroblast", "Others", "Others", "Osteogenic cells", "Others")

names(new.id) <- levels(nc_E14.5)
nc_E14.5 <- RenameIdents(nc_E14.5,new.id)
DimPlot(nc_E14.5)


## E16.5----------------
nc_E16.5 <- nc_sub[,nc_sub$orig.ident=="ToothNc_Junjun_E16.5"]
nc_E16.5 <- runSeurat(nc_E16.5)
FeaturePlot(nc_E16.5, c("Tfap2b","Top2a","Fibin",
                        "Igf1","Zfhx3","Scx",
                        "Alpl","Ibsp"),label = T)
FeaturePlot(nc_E16.5, c("Aspn","Lmo1","Msx1",
                        "Mfap5","Eln","Top2a"),label = T)
FeaturePlot(nc_E16.5, c("Lmo1","Fst","Lepr","Aldh1a2"),label = T)
FeaturePlot(nc_E16.5, c("Phex"),label = T)

newid_E16.5 <- c("Coronal papilla", "Dermal fibroblast", "Apical follicle",
                 "Dermal fibroblast", "Apical follicle", "Lateral follicle", "Lateral follicle",
                 "Dermal fibroblast", "Dermal fibroblast", "Dermal fibroblast",
                 "Odontoblast cells", "Dermal fibroblast",
                 "Others", "Dermal fibroblast", "Dermal fibroblast",
                 "Tenogenic cells", "Dermal fibroblast", "Dermal fibroblast", "Dermal fibroblast")

names(newid_E16.5) <- levels(nc_E16.5)
nc_E16.5 <- RenameIdents(nc_E16.5,newid_E16.5)

nc_E16.5 <- FindSubCluster(nc_E16.5,cluster = "Coronal papilla",graph.name = "RNA_snn")
DimPlot(nc_E16.5,group.by  = "sub.cluster")
Idents(nc_E16.5) <- nc_E16.5$sub.cluster
newid_E16.5_2 <- c("Dermal fibroblast", "Lateral follicle", "Coronal papilla",
                   "Coronal papilla", "Apical follicle", "Odontoblast cells",
                   "Apical papilla", "Others", "Tenogenic cells", "Apical papilla",
                   "Coronal papilla", "Coronal papilla", "Coronal papilla"
)

names(newid_E16.5_2) <- levels(nc_E16.5)
nc_E16.5 <- RenameIdents(nc_E16.5,newid_E16.5_2)
DimPlot(nc_E16.5,label = T)

#== P35---------------
nc_P3.5 <- nc_sub[,nc_sub$orig.ident=="ToothNc_Junjun_P3.5"]
nc_P3.5 <- runSeurat(nc_P3.5)
DimPlot(nc_P3.5,label = T)
FeaturePlot(nc_P3.5,c("Nnat","Enpp6","Pcp4",
                      "Tnmd","Dspp","Smoc2"),label = T)
FeaturePlot(nc_P3.5,c("Snai2"),label = T)
FeaturePlot(nc_P3.5,c("Slc1a3"),label = T)

new_id_P3.5 <- c("Coronal papilla", "Middle papilla",
                 "Apical papilla",
                 "Middle papilla", "Apical follicle",
                 "Apical papilla", "Lateral follicle",
                 "Apical papilla", "Others", "Others",
                 "Others", "Coronal papilla",
                 "Lateral follicle")
names(new_id_P3.5) <- levels(nc_P3.5)
nc_P3.5 <- RenameIdents(nc_P3.5,new_id_P3.5)
DimPlot(nc_P3.5)

#== P75----------------

nc_P7.5 <- nc_sub[,nc_sub$orig.ident=="ToothNc_Junjun_P7.5"]
nc_P7.5 <- runSeurat(nc_P7.5)
FeaturePlot(nc_P7.5,c("Nnat","Enpp6","Pcp4",
                      "Tnmd","Dspp","Smoc2"),label = T)
DimPlot(nc_P7.5,label = T)
FeaturePlot(nc_P7.5,c("Nnat","Enpp6","Pcp4",
                      "Tnmd","Dspp","Smoc2"),label = T)
newid_P7.5 <- c("Middle papilla", "Coronal papilla", "Apical follicle",
                "Lateral follicle", "Odontoblast", "Apical papilla",
                "Apical papilla", "Others", "Others", "Lateral follicle", "Lateral follicle")
names(newid_P7.5) <- levels(nc_P7.5)
nc_P7.5 <- RenameIdents(nc_P7.5,newid_P7.5)
DimPlot(nc_P7.5)
cluster1 <- Idents(nc_E13.5) %>% as.data.frame()
cluster2 <- Idents(nc_E14.5) %>% as.data.frame()
cluster3 <- Idents(nc_E16.5) %>% as.data.frame()
cluster4 <- Idents(nc_P3.5) %>% as.data.frame()
cluster5 <- Idents(nc_P7.5) %>% as.data.frame()
write.csv(cluster1,"process/annotation/20250114_annotation_harmonized/20250114_NC_junjun_E135_harmonized.csv")
write.csv(cluster2,"process/annotation/20250114_annotation_harmonized/20250114_NC_junjun_E145_harmonized.csv")
write.csv(cluster3,"process/annotation/20250114_annotation_harmonized/20250114_NC_junjun_E165_harmonized.csv")
write.csv(cluster4,"process/annotation/20250114_annotation_harmonized/20250114_NC_junjun_P35_harmonized.csv")
write.csv(cluster5,"process/annotation/20250114_annotation_harmonized/20250114_NC_junjun_P75_harmonized.csv")

clusterDf <- list(cluster1,cluster2,cluster3,cluster4,cluster5)
clusterBind <- do.call(rbind, clusterDf)
clusterBind$. <- as.character(clusterBind$.)
clusterBind[clusterBind == "Odontoblast cells"] = "Osteoblast"
clusterBind[clusterBind == "Cycling cell"] = "Cycling cells"

nc_sub$curate_anno <- clusterBind[colnames(clusterBind)]


DimPlot(nc_sub, group.by = "curate_anno",label = T,raster = T)
ggsave("results/annotation/anno_harmonize/20250119_mes/20250119_ncJunjun_combine.pdf",width = 6,height = 6)
saveRDS(nc_sub,"processed_data/annotation_harmonization/20250119_ncjunjun.Rds")
#
# all_cluster <- c()
#
#
# #mes <- readRDS("processed_data/integrated_data/20241106_mesenchyme.Rds")
# FeaturePlot(mes,c("Ibsp","Col2a1","Cdh5","Top2a","Fibin","Akap12"))
