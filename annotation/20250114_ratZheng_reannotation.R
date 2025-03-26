rat <- readRDS("processed_data/preprocess_data/Tooth_Zheng.Rds")
FeaturePlot(rat,"Msx1")
FeaturePlot(rat,c("Igfbp5","Kit","Gsc",
                  "Postn","Smoc2","Smpd3"),label = T,)
rat <- FindSubCluster(rat,cluster = "3",
                      graph.name = "RNA_snn")
Idents(rat) <- rat$sub.cluster
DimPlot(rat,group.by = "sub.cluster")

newid <- c("Peridontal ligment", "Dental follicles",
           "Odontoblast", "Others",
           "Dental pulp", "Dental pulp", "Odontoblast",
           "Apical pulp", "Dental follicles", "Others", "Others",
  "Dental follicles", "Pre-odontoblasts", "Pre-odontoblasts", "Dental pulp",
  "Pre-odontoblasts", "Others")
names(newid) <- levels(rat)
rat <- RenameIdents(rat, newid)
DimPlot(rat, label = T)
rat$label <- Idents(rat)
#rat_mes <- rat[,rat$label!= "Others"]
rat_mes <- subset(x = rat, subset = label != "Others")
#rat_mes <- JoinLayers(rat_mes)

rat_cluster <- rat_mes$label %>% as.data.frame()
write.csv(rat_cluster,"process/annotation/20250114_annotation_harmonized/20250118_rat_annotation.csv")
#rat_cluster <- read.csv("process/annotation/20250114_annotation_harmonized/20250118_rat_annotation.csv",row.names = 1)
rat_mes = rat[,rownames(rat_cluster)]
rat_mes$curate_anno <- rat_cluster$.
#rat_mes <- rat[,rat$label!= "Others"]
DimPlot(rat_mes,group.by = "label")
DimPlot(rat,group.by = "label")
write.csv(as.matrix(rat@reductions$umap), )
