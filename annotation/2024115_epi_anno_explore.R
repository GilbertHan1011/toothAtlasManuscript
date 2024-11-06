rm(list=ls())
library(Seurat)
library(tidyverse)

source("script/utils/anno_utils.R")
source("script/utils/stratified_wilcoxon_functions.R")
source("script/utils/anno_integrate.R")
process_dir <- "process/annotation/epi_annotation/annotation_mid_file/"
load_json("script/annotation/config/epithelium_parameter.json")
full_seurat <- readRDS("processed_data/integrated_data/20241030_epithelium.Rds")
updated_cluster_object <- readRDS(paste0(parameter_list$harmonization_folder_path,"annotation_mid_file/",
                                         parameter_list$new_name_suffix,"_pruned_mrtree_clustering_results",".rds"))
edgelist_updated_new_labels <- updated_cluster_object$edgelist
labelmat_updated_new_labels <- updated_cluster_object$labelmat
genes_to_include <- read.csv("process/annotation/epi_annotation/annotation_mid_file/20241030_gene_include.csv",row.names = 1) %>% unlist

full_seurat$C5 <- labelmat_updated_new_labels[,1]
full_seurat$C12 <- labelmat_updated_new_labels[,2]
full_seurat$C22 <- labelmat_updated_new_labels[,3]
full_seurat$C31 <- labelmat_updated_new_labels[,4]
full_seurat$C43 <- labelmat_updated_new_labels[,5]
full_seurat$C73 <- labelmat_updated_new_labels[,6]
DimPlot(full_seurat,group.by = "C5") + DimPlot(full_seurat,group.by = "C12") + DimPlot(full_seurat,group.by = "C22") + DimPlot(full_seurat,group.by = "C31") + DimPlot(full_seurat,group.by = "C43") + DimPlot(full_seurat,group.by = "C73")

DimPlot(full_seurat,group.by = "C5") + FeaturePlot(full_seurat,c("Slco4a1","Th","Amer1")) + FeaturePlot(full_seurat,c("Vat1l","Fam19a4","Hey2"))
DimPlot(full_seurat,group.by = "C5") + FeaturePlot(full_seurat,c("Rab3il1","Pmch","Cyp2s1",'Psmb10','C1qb','Ibsp'))
DimPlot(full_seurat,group.by = "C5") + FeaturePlot(full_seurat,c("Amelx","Ambn"))
DimPlot(full_seurat,group.by = "C5") +FeaturePlot(full_seurat,'Vwde')
DimPlot(full_seurat,group.by = "C12") + DimPlot(full_seurat,group.by = "Stage")
DimPlot(full_seurat,group.by = "C12") + FeaturePlot(full_seurat,c("Slco4a1","Th","Amer1")) + FeaturePlot(full_seurat,c("Vat1l","Fam19a4","Hey2"))
DimPlot(full_seurat,group.by = "C12") + DimPlot(full_seurat,group.by = "Development.stage")

DimPlot(full_seurat,group.by = "C12", label = TRUE) + FeaturePlot(full_seurat,c('Rab3il1', 'Pmch', 'Cyp2s1', 'Psmb10', 'C1qb', 'Ibsp'))
DimPlot(full_seurat,group.by = "C5", label = TRUE)
DimPlot(full_seurat,group.by = "C12", label = TRUE) 
DimPlot(full_seurat,group.by = "C22", label = TRUE)
DimPlot(full_seurat,group.by = "C31", label = TRUE)
DimPlot(full_seurat,group.by = "C43", label = TRUE)
DimPlot(full_seurat,group.by = "C73", label = TRUE)

DimPlot(full_seurat,group.by = "C12", label = TRUE) + FeaturePlot(full_seurat,c("Amelx","Ambn"))
DimPlot(full_seurat,group.by = "C12", label = TRUE) + FeaturePlot(full_seurat,c("Ryr2", "Mylk", "Sox5"))
DimPlot(full_seurat,group.by = "Development.stage", label = TRUE)

DimPlot(full_seurat,group.by = "C22", label = TRUE) + FeaturePlot(full_seurat,Outer_enamel_epithelium <- c('Foxp1', 'Ntrk2', 'Tbx1'))

FeaturePlot(full_seurat,c('Isl1', 'Cldn10', 'Sfrp5', 'Pthlh', 'Grp', 'Shisa2', 'Crabp1', 'Fam19a4', 'Shh', 'Aplf', 'Krt18'))
FeaturePlot(full_seurat,c('Isl1', 'Cldn10', 'Pthlh', 'Grp', 'Shisa2'), ncol = 3)
FeaturePlot(full_seurat,c("Rab3il1", "Pmch", "Cyp2s1", "Psmb10", "C1qb", "Ibsp"), ncol = 3, cols = c("lightgrey", "blue", "red"), max.cutoff = "q99", min.cutoff = "q90")
FeaturePlot(full_seurat, c('Mmp20', 'Enam', 'Lama2', 'Nrn1l', 'Cd55', 'Plod2', 'Galnt12', 'Cd24a'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
DimPlot(full_seurat,group.by = "Sample", label = TRUE)
FeaturePlot(full_seurat, c('Klk4', 'Amtn', 'Odam', 'Csn3'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
FeaturePlot(full_seurat, c('Dspp'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
DimPlot(full_seurat,group.by = "C22", label = TRUE) + FeaturePlot(full_seurat, c("Col22a1", "Vwde", "Kif5c"), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
DimPlot(full_seurat,group.by = "C22", label = TRUE) + FeaturePlot(full_seurat, c('Notch1', 'Notch2'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
DimPlot(full_seurat,group.by = "C22", label = TRUE) + FeaturePlot(full_seurat, c('Pthlh'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")

FeaturePlot(full_seurat, c('Krt17', 'Tacstd2', 'Sfn', 'Tagln', 'Igfbp3', 'Pthlh'), cols = c("lightgrey", "blue", "red"))
FeaturePlot(full_seurat, "Atf3", cols = c("lightgrey", "blue", "red"))
FeaturePlot(full_seurat, c('Tacstd2'))
FeaturePlot(full_seurat, c('Tacstd2'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
FeaturePlot(full_seurat,  c("Pcdh7"))
FeaturePlot(full_seurat,  c('Hba-x'))

umap_data <- as.data.frame(Embeddings(full_seurat, "X_umap"))
umap_data$C22G3 <- ifelse(full_seurat@meta.data[["C22"]] == "C22-3", "C22-3", "Unselected")

ggplot(umap_data, aes(x = Xumap_1, y = Xumap_2)) +
  geom_point(aes(color = C22G3), size = 0.5) +
  scale_color_manual(values = c("Unselected" = "grey", "C22-3" = "red")) +  
  labs(title = "highlight", color = "") +
  theme_minimal() +
  theme(legend.position = "right") + FeaturePlot(full_seurat, c('Scube3'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")

umap_data$Stage <- ifelse(full_seurat@meta.data[["Stage"]] == "Young Adult", "Young Adult", "Unselected")

ggplot(umap_data, aes(x = Xumap_1, y = Xumap_2)) +
  geom_point(aes(color = Stage), size = 0.5) +
  scale_color_manual(values = c("Unselected" = "grey", "Young Adult" = "red")) +  
  labs(title = "highlight", color = "") +
  theme_minimal() +
  theme(legend.position = "right") 

library(SeuratWrappers)
cds <- as.cell_data_set(full_seurat)
head(colData(cds))

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

full_seurat@active.ident
list.cluster <- setNames(as.vector(full_seurat@meta.data[["leiden_clusters_level_3"]]), rownames(full_seurat@meta.data))
# list.cluster <- full_seurat@meta.data[["leiden_clusters_level_3"]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- full_seurat@reductions$X_umap@cell.embeddings
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, reduction_method = 'UMAP',
                                 group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
# cds <- order_cells(cds)

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 15]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)


#Level1
DimPlot(full_seurat,group.by = "Stage", label = TRUE)
#OEE
FeaturePlot(full_seurat, c("Slco4a1", "Th"), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
#SI
FeaturePlot(full_seurat, c("Thbd", "Jph4", 'Notch1'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
#IEE
FeaturePlot(full_seurat, c("Sfrp5", 'Cldn10'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
#SR
FeaturePlot(full_seurat,  c("Scube3", "Fam19a4"), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
#pAM
FeaturePlot(full_seurat, c("Col22a1", "Vwde", "Kif5c"), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
#eAM
FeaturePlot(full_seurat, c('Dspp'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
#sAM
FeaturePlot(full_seurat, c('Enam', 'Cd55', 'Galnt12'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
#mAM
FeaturePlot(full_seurat, c('Klk4', 'Amtn', 'Odam'), cols = c("lightgrey", "blue", "red"), max.cutoff = "q99")
#Evidence for Low Quality cells
DimPlot(full_seurat,group.by = "Sample", label = TRUE) + DimPlot(full_seurat,group.by = "C5", label = TRUE)
DotPlot(full_seurat, group.by = "C22", features = 'Hba-x') + RotatedAxis()
#VlnPlot(full_seurat, group.by = "C22", features = 'Hba-x') + NoLegend()
#Questions
DimPlot(full_seurat,group.by = "Stage", label = TRUE) + DimPlot(full_seurat,group.by = "Sample", label = TRUE)
