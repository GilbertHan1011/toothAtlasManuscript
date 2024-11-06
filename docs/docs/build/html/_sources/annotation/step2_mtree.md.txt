# Step 2: MRTree

## Introduction
In this step, we will use the MRTree algorithm {cite:p}`pengCellTypeHierarchy2021` to build a hierarchical clustering tree from the input labels. [MRTree](https://github.com/pengminshi/MRtree) is a powerful tool using multiple resolution clustering results to build a hierarchical tree.

## Input
- `20241025_mes_anno_mrtree_input_label.csv`: This file contains the input labels for the MRTree algorithm, which is derived from the [Step1](./step1_choose_optimal_cluster_level.md).
- `mesenchyme_parameter.json`: This file contains the parameters for the MRTree algorithm.

## Output
- `annotation_mid_file/mes_mrtree.Rds`: This file contains the raw result of the MRTree algorithm.
- `20241025_mes_cluster_object.Rds`: This file contains the final cluster object, including the label matrix, edge list, node list, and data tree.

## Script
```R
#== Set Environment ------------------
rm(list = ls())
library(mrtree)
library(dplyr)
library(data.tree)
library(data.table)
library(jsonlite)

source("script/utils/mrtree_fun.R")

process_dir <- "process/annotation/mes_annotation/annotation_mid_file/"
mrtree_input_labels <- fread("processed_data/metadata/20241025_mes_anno_mrtree_input_label.csv")

# Convert input labels to numeric matrix
cluster_matrix_for_mrtree <- as.matrix(mrtree_input_labels)
cluster_matrix_for_mrtree <- apply(cluster_matrix_for_mrtree, 2, as.numeric)

message(Sys.time(), ": Build mrtree using matrix with: ", dim(cluster_matrix_for_mrtree)[1], " cells and ", dim(cluster_matrix_for_mrtree)[2], " cluster levels.")

# Feed into mrtree
mrtree_res <- mrtree(
  cluster_matrix_for_mrtree,
  prefix = "leiden_clusters_level_",
  suffix = NULL,
  max.k = Inf,
  consensus = FALSE,
  sample.weighted = FALSE,
  augment.path = FALSE,
  verbose = TRUE,
  n.cores = 5
)

saveRDS(mrtree_res, paste0(process_dir, "mes_mrtree.Rds"))

# Load and process parameters
parameter_list <- read_json("script/annotation/config/mesenchyme_parameter.json")
parameter_list <- lapply(parameter_list, function(x) {
  if (is.list(x)) {
    return(unlist(x))
  } else {
    return(x)
  }
})

message(Sys.time(), ": Saved raw result of length ", length(mrtree_res), " to: ", paste0(parameter_list$harmonization_folder_path, parameter_list$new_name_suffix, "_mrtree_res_raw.rds"))

message(Sys.time(), ": Create mrtree output list")

# Make label matrix with column names included
if (parameter_list$use_recon_labelmat) {
  labelmat <- mrtree_res$labelmat.recon
  Ks.recon <- apply(labelmat, 2, function(y) length(unique(y)))
  unique.idx <- 1:length(Ks.recon)
  colnames(labelmat) <- paste0("K", Ks.recon[unique.idx])
} else {
  labelmat <- mrtree_res$labelmat.mrtree
}

# Build dataframe with labels
n <- nrow(labelmat)
backup <- colnames(labelmat)
labelmat <- matrix(paste(matrix(rep(colnames(labelmat), each = n), nrow = n), labelmat, sep = '-'), nrow = n)
colnames(labelmat) <- backup
df <- as.data.frame(unique(labelmat), stringsAsFactors = FALSE)

# Save in data.tree format
df$pathString <- apply(df, 1, function(x) paste(c('all', x), collapse = '/'))
tree.datatree <- as.Node(df)

# Export edgelist and nodelist from data.tree
edges <- ToDataFrameNetwork(tree.datatree, "isLeaf", "level", "count", "totalCount", "height")
nodes <- data.frame(id = c("all", as.character(unique(edges$to))), label = c("all", as.character(unique(edges$to))))
nodes <- rbind(c("all", FALSE, 1, 2, 228, max(edges$height) + 1), edges[, 2:ncol(edges)]) %>% 
  rename(id = to) %>% 
  mutate(label = id)

# Make cluster object (list)
cluster_object <- list(
  labelmat = labelmat,
  edgelist = edges,
  nodelist = nodes,
  data_tree = tree.datatree,
  mrtree_output = mrtree_res
)

# Save cluster object
saveRDS(cluster_object, paste0(process_dir, "20241025_mes_cluster_object.Rds"))

message("Finalized")
```
