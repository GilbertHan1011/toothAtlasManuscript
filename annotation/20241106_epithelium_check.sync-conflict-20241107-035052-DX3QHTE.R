## check epithelium
outdir <- "results/annotation/epi_anno_explore/"
dir.create(outdir)
epi <- readRDS("processed_data/integrated_data/20241106_epithelum.Rds")
DimPlot(epi)


clusterObj <- readRDS("process/annotation/epi_annotation/annotation_mid_file/20241030_epi_cluster_object.Rds")
mrtree <- readRDS("process/annotation/epi_annotation/annotation_mid_file/epi_mrtree.Rds")
cluster_res <- readRDS("process/annotation/epi_annotation/annotation_mid_file/tooth_epi_harmonized__pruned_mrtree_clustering_results.rds")

epi$C5 <- cluster_res$labelmat[,1]
epi$C12 <- cluster_res$labelmat[,2]
epi$C22 <- cluster_res$labelmat[,3]
epi$C31 <- cluster_res$labelmat[,4]
epi$C43 <- cluster_res$labelmat[,5]
epi$C73 <- cluster_res$labelmat[,6]
epi <- NormalizeData(epi,normalization.method = "LogNormalize")
DimPlot(epi,group.by = "C5")

FeaturePlot(epi,"Odam")
FeaturePlot(epi,"Cd55")
FeaturePlot(epi,"Dspp")
FeaturePlot(epi,"Jph4")


## C5---------
p1 <- DimPlot(epi,group.by = "C5")
p2 <- FeaturePlot(epi,"Ambn")
p1|p2
ggsave(paste0(outdir,"C5_ambn.pdf"),width = 8,height = 4)

## C12--------------
mycolor <- colorRampPalette(brewer.pal(8,'Spectral'))(16)
p1 <- DimPlot(epi,group.by = "Development.stage")
p2 <- DimPlot(epi,group.by = "Age")+scale_color_manual(values=mycolor)
p3 <- DimPlot(epi,group.by = "C12",cells.highlight = colnames(epi)[epi$C12=="C12-1"])
p2|p3
ggsave(paste0(outdir,"C12-1.pdf"),width = 8,height = 4)


epi$C12 <- cluster_res$labelmat[,2]

FeaturePlot(epi, c("Sfrp5",
                   "Vwde",
                   "Cyp26a1",
                   "Shh",
                   "Igfbpl1"
))
FeaturePlot(epi, c("Cdh6","Lrp11"
))

FeaturePlot(epi, c("Fos","Egr1"
))

FeaturePlot(epi, c("Vat1l","Fam19a4"
))
FeaturePlot(epi, c("Thbd","Gnrh1"
))
FeaturePlot(epi, c("Krt17"
))
