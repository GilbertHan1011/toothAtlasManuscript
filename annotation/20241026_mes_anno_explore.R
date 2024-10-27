fig_dir <- "results/annotation/anno_explore/"
dir.create(fig_dir)
# C2-----------------------
p1 <- DimPlot(full_seurat,group.by = "Stage")
p2 <- DimPlot(full_seurat,group.by = "C2")
p1|p2
ggsave(paste0(fig_dir,"C2_age.pdf"),height = 4,width = 8)

full_seurat@assays$originalexp@data
FeaturePlot(full_seurat,"Prrx1",,slot = "data")
#saveRDS(full_seurat,"processed_data/integrated_data/20241026_mesenchyme.Rds")

# C9---------------------------
anno_marker <- read.csv("data/annotation/annotation_marker.csv")
anno_marker <- anno_marker[,c("Clusters","Marker.genes")]
split_list <- split(anno_marker, anno_marker$Clusters)

# Convert each data frame in the list to a vector of marker genes
result_list <- lapply(split_list, function(x) x$Marker.genes)
result_list <- result_list[2:58]
for (i in names(result_list)) {
  genes <- result_list[[i]] %>% unlist
  res_dir <- paste0(fig_dir, i)
  if (!dir.exists(res_dir)) {
    dir.create(res_dir, recursive = TRUE)
  }
  for (gene in genes) {
    tryCatch({
      p <- FeaturePlot(full_seurat, features = gene)
      ggsave(paste0(res_dir, "/", gene, ".pdf"), plot = p, width = 4, height = 4)
    }, error = function(e) {
      message(paste("Error processing gene:", gene, "- Error:", e$message))
    })
  }
}

DimPlot(full_seurat,group.by = "C9")
FeaturePlot(full_seurat,c("Meg3","Osr1","Egfl6"))

### C9-8 odontoblast

