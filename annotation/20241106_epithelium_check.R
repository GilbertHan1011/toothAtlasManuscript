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

p1 <- FeaturePlot(epi,c("Slco4a1","Th","Thbd","Jph4","Sfrp5"),raster = T)
p2 <- DimPlot(epi,group.by = "C12",cells.highlight = colnames(epi)[epi$C12=="C12-2"],raster = T)
p1+p2
ggsave(paste0(outdir,"C12-2.pdf"),width = 6,height = 8)


p2 <- FeaturePlot(epi,c("Vat1l","Fam19a4"))
p3 <- DimPlot(epi,group.by = "C12",cells.highlight = colnames(epi)[epi$C12=="C12-3"])
p2+p3
ggsave(paste0(outdir,"C12-3.pdf"),width = 6,height = 6)

p2 <- FeaturePlot(epi,c("Col22a1","Vwde","Kif5c"))
p3 <- DimPlot(epi,group.by = "C12",cells.highlight = colnames(epi)[epi$C12=="C12-5"])
p2+p3
ggsave(paste0(outdir,"C12-5.pdf"),width = 6,height = 6)


p2 <- FeaturePlot(epi,c("Odam","Klk4","Amtn"))
p3 <- DimPlot(epi,group.by = "C12",cells.highlight = colnames(epi)[epi$C12=="C12-10"])
p2+p3
ggsave(paste0(outdir,"C12-10.pdf"),width = 6,height = 6)


p2 <- FeaturePlot(epi,c("Dspp"))
p3 <- DimPlot(epi,group.by = "C22",cells.highlight = colnames(epi)[epi$C22=="C22-12"])
p2+p3
ggsave(paste0(outdir,"C22-12.pdf"),width = 6,height = 6)


FeaturePlot(epi,c("Slc5a8","Gm17660","Ptpn22"))

p2 <- FeaturePlot(epi,c("Enam","Cd55","Galnt12"))
p3 <- DimPlot(epi,group.by = "C12",cells.highlight = colnames(epi)[epi$C12=="C12-9"])
p2+p3
ggsave(paste0(outdir,"C12-9.pdf"),width = 6,height = 6)

p2 <- FeaturePlot(epi,c("Dspp"))
p3 <- DimPlot(epi,group.by = "C22",cells.highlight = colnames(epi)[epi$C22=="C22-12"])
p2+p3
ggsave(paste0(outdir,"C22-12.pdf"),width = 6,height = 6)


p1 <- FeaturePlot(epi,c("Slco4a1","Th"),raster = T)
p2 <- DimPlot(epi,group.by = "C22",cells.highlight = colnames(epi)[epi$C22=="C22-7"],raster = T)
p1+p2
ggsave(paste0(outdir,"C22-2.pdf"),width = 6,height = 6)


p1 <- FeaturePlot(epi,c("Thbd","Jph4"),raster = T)
p2 <- DimPlot(epi,group.by = "C22",cells.highlight = colnames(epi)[epi$C22=="C22-5"],raster = T)
p1+p2
ggsave(paste0(outdir,"C22-5.pdf"),width = 6,height = 6)

p1 <- FeaturePlot(epi,c("Ibsp","Enpp2","Cldn10"),raster = T)
p2 <- DimPlot(epi,group.by = "C22",cells.highlight = colnames(epi)[epi$C22=="C22-4"],raster = T)
p1+p2
ggsave(paste0(outdir,"C22-4.pdf"),width = 6,height = 6)


FeaturePlot(epi,c("Rab3il1",
                  "Pmch",
                  "Cyp2s1",
                  "Psmb10"
))

FeaturePlot(epi,c("Enpp2",
                  "C1qb",
                  "Ibsp",
                  "Rhov",
                  "Nrarp",
                  "Tactd2"
))
