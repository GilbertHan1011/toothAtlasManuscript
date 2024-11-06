# Step 5: Annotation

## Introduction
In this step, we will annotate the mesenchymal cell clusters identified in the previous steps. The annotation process involves updating the Seurat object with new cluster labels, finding marker genes for the updated clusters, and calculating specificity for comparisons. We will also manually annotate the clusters and update marker scores.

## Input
- `processed_data/integrated_data/20241024_mesenchyme.Rds`: The Seurat object containing integrated mesenchymal cell data.
- `process/annotation/mes_annotation/annotation_mid_file/20241025_gene_include.csv`: CSV file containing genes to include in the analysis.
- `process/annotation/mes_annotation/annotation_mid_file/[new_name_suffix]_pruned_mrtree_clustering_results.rds`: The pruned clustering results from Step 4.
- `script/annotation/config/mesenchyme_parameter.json`: Configuration file with parameters for annotation.

## Output
- `process/annotation/mes_annotation/annotation_mid_file/comparisons_all_update.csv`: CSV file containing all marker comparisons.
- `process/annotation/mes_annotation/annotation_mid_file/comparisons_sibling_update.csv`: CSV file containing sibling marker comparisons.
- `result/4.26_annotation/4.26_annotation_orig.csv`: CSV file containing the original annotation data.
- `processed_data/integrated_data/20241106_mesenchyme.Rds`: The updated Seurat object with annotations.
- `process/annotation/mes_annotation/annotation_mid_file/[new_name_suffix]_[marker_suffix]_annotation_result.txt`: Text file containing the annotation results.
- `process/annotation/mes_annotation/annotation_mid_file/[new_name_suffix]_[marker_suffix]_annotation_labelmat.txt`: Text file containing the annotation label matrix.
- `process/annotation/mes_annotation/annotation_mid_file/[new_name_suffix]_[marker_suffix]_annotation_markers_filtered.txt`: Text file containing the filtered annotation markers.
- `process/annotation/mes_annotation/annotation_mid_file/[new_name_suffix]_[marker_suffix]_anno_bind_df.csv`: CSV file containing the bound annotations.

## Script
```R
#== Change age label---------------------------------
# Clear the workspace
rm(list = ls())

# Load required libraries
library(Seurat)
library(tidyverse)

# Load utility functions
source("script/utils/anno_utils.R")
source("script/utils/stratified_wilcoxon_functions.R")
source("script/utils/anno_integrate.R")

# Define the directory for processed data
process_dir <- "process/annotation/mes_annotation/annotation_mid_file/"

# Load JSON parameters
load_json("script/annotation/config/mesenchyme_parameter.json")

# Load the Seurat object containing integrated mesenchymal cell data
full_seurat <- readRDS("processed_data/integrated_data/20241024_mesenchyme.Rds")

# Load the updated cluster object
updated_cluster_object <- readRDS(paste0(parameter_list$harmonization_folder_path, "annotation_mid_file/",
                                         parameter_list$new_name_suffix, "_pruned_mrtree_clustering_results.rds"))

# Extract edgelist and labelmat from the updated cluster object
edgelist_updated_new_labels <- updated_cluster_object$edgelist
labelmat_updated_new_labels <- updated_cluster_object$labelmat

# Load genes to include
genes_to_include <- read.csv("process/annotation/mes_annotation/annotation_mid_file/20241025_gene_include.csv", row.names = 1) %>% unlist

# Update Seurat object with new labels
full_seurat$C2 <- labelmat_updated_new_labels[, 1]
DimPlot(full_seurat, group.by = "C2")

#== Find markers--------------------------
# Find markers for the updated clusters
all_markers_mrtree_list_update <- findMarkers_tree2(
  full_seurat,
  edgelist = edgelist_updated_new_labels[, 1:2],
  labelmat = labelmat_updated_new_labels,
  n_cores = 15,
  assay = assay_markers,
  slot = assay_slot,
  test.use = test.use,
  batch_var = batch_var,
  latent_vars_seurat = latent_vars_seurat,
  genes_to_include = genes_to_include,
  logfc.threshold = logfc.threshold,
  min.pct = min.pct,
  min.diff.pct = min.diff.pct,
  max.cells.per.ident = max.cells.per.ident,
  min.cells.feature = min.cells.feature,
  min.cells.group = min.cells.group,
  base = base,
  only.pos = only.pos
)

# Extract comparisons from the marker results
comparisons_all_update <- all_markers_mrtree_list_update$comparisons_All
comparisons_siblings_update <- all_markers_mrtree_list_update$comparisons_Sibling

# Calculate specificity for comparisons
comparisons_all_update$specificity <- ((comparisons_all_update$pct.1 + specificity_base) /
                                        (comparisons_all_update$pct.2 + specificity_base)) * comparisons_all_update$avg_log2FC
comparisons_siblings_update$specificity <- ((comparisons_siblings_update$pct.1 + specificity_base) /
                                             (comparisons_siblings_update$pct.2 + specificity_base)) * comparisons_siblings_update$avg_log2FC

# View and save comparisons
View(comparisons_siblings_update)
write.csv(comparisons_all_update, paste0(process_dir, "comparisons_all_update.csv"))
write.csv(comparisons_siblings_update, paste0(process_dir, "comparisons_sibling_update.csv"))

# Load saved comparisons
comparisons_siblings <- read.csv(paste0(process_dir, "comparisons_sibling_update.csv"), row.names = 1)
comparisons_all <- read.csv(paste0(process_dir, "comparisons_all_update.csv"), row.names = 1)

#== Label--------------------------------
# Update Seurat object with new labels for different cluster levels
full_seurat$C2 <- labelmat_updated_new_labels[, 1]
full_seurat$C9 <- labelmat_updated_new_labels[, 2]
full_seurat$C19 <- labelmat_updated_new_labels[, 3]
full_seurat$C29 <- labelmat_updated_new_labels[, 4]
full_seurat$C57 <- labelmat_updated_new_labels[, 5]
full_seurat$C92 <- labelmat_updated_new_labels[, 6]

# Plot the clusters
DimPlot(full_seurat, group.by = "C2", label = TRUE)
DimPlot(full_seurat, group.by = "C9", label = TRUE)
DimPlot(full_seurat, group.by = "C19", label = TRUE)
DimPlot(full_seurat, group.by = "C29", label = TRUE)
DimPlot(full_seurat, group.by = "C57", label = TRUE)
DimPlot(full_seurat, group.by = "C92", label = TRUE)

# Define groups for annotation binding
group_select <- c("Mandibular_Maxillary", "Molar_Incisor", "Histology", "Treatment",
                  "Age", "Stage", "Development.stage", "coarse_anno_1")

# Bind annotations
annoBindDf <- anno_bind(full_seurat, c("C2", "C9", "C19", "C29", "C57", "C92"), group_select)
percentPlotFun(full_seurat, "C19", "Age")

# Create a table of cluster counts by age
df <- table(full_seurat$C6, full_seurat$Age.In.Detail.) %>% t %>% as.matrix()
df <- t(apply(df, 1, function(x) x / sum(x)))
df <- as.data.frame(df)

# Convert the table to long format
df_long <- df %>% rownames_to_column("Age") %>%
  pivot_longer(cols = -Age, names_to = "Cluster")

# Identify genes to remove based on patterns
other_genes_remove <- rownames(full_seurat@assays$originalexp@counts)[grepl("RP|Gm|Hb[ab]|Rik|-ps", rownames(full_seurat@assays$originalexp@counts))]
all_exclusion_genes <- unique(c(other_genes_remove))

##########
### Load parameters and packages
##########

# Load the MRTree clustering results
mrtree_result <- updated_cluster_object

# Extract edgelist and labelmat
edgelist <- mrtree_result$edgelist
labelmat <- mrtree_result$labelmat

# Process edgelist and cluster levels
edgelist <- edgelist[, c("from", "to", "level")]
cluster_levels <- as.data.frame(labelmat) %>%
  tidyr::pivot_longer(everything(), names_to = "clusterlevel", values_to = "cluster") %>%
  dplyr::group_by(cluster) %>%
  dplyr::add_count(name = "ncells") %>%
  dplyr::distinct(clusterlevel, cluster, ncells)
edgelist <- dplyr::left_join(edgelist, cluster_levels, by = c("to" = "cluster")) %>%
  dplyr::arrange(level)
all_nodes <- unique(edgelist$to)

# Load marker comparisons
markers_comparisons_siblings <- comparisons_siblings
markers_comparisons_all <- comparisons_all

##########
### Run annotation
##########

# Print starting message and parameters
message("Starting annotation: ")
print(parameter_list$manual_names_annotation_abbr)

# Annotate the tree
annotation_results <- annotate_tree(
  edgelist = edgelist,
  labelmat = labelmat_updated_new_labels,
  markers_comparisons_all = markers_comparisons_all,
  markers_comparisons_siblings = markers_comparisons_siblings,
  preferred_genes = character(0),
  manual_names = parameter_list$manual_names_annotation_abbr,
  overwrite_with_manual = TRUE,
  manual_exclude_genes = all_exclusion_genes,
  max_pval_adj = parameter_list$max_pval_adj,
  min_specificity = parameter_list$min_specificity,
  min_specificity_sibling_children = parameter_list$min_specificity_sibling_children,
  scale_preferred = 1,
  limit_factor = parameter_list$limit_factor,
  max_score_siblings_children = parameter_list$max_score_siblings_children,
  reverse_order = parameter_list$reverse_order
)

# Print message for formatting results
message("Formating annotation results")

# Get annotation results
annotation_df <- annotation_results$annotation_df
annotation_df$clean_names_withID <- paste0(annotation_df$cluster_id, ": ", annotation_df$clean_names)

# Update labelmat with cell IDs
labelmat_updated <- cbind(Cell_ID = full_seurat@meta.data$Cell_ID, mrtree_result$labelmat)

# Create a per cell cluster annotation
cell_cluster_map <- labelmat_updated %>%
  as.data.frame() %>%
  tidyr::gather(-Cell_ID, key = "clusterlevel", value = "cluster_id")
cell_cluster_map$clusterlevel <- gsub("_pruned", "", cell_cluster_map$clusterlevel)

# Create a wide version of the labelmat for Seurat metadata
annotation_df_wide <- annotation_df %>%
  dplyr::left_join(cell_cluster_map, by = c("clusterlevel" = "clusterlevel", "cluster_id" = "cluster_id")) %>%
  dplyr::select(Cell_ID, clusterlevel, clean_names_withID) %>%
  tidyr::spread(key = clusterlevel, value = clean_names_withID)

# Update column names and sort them
colnames(annotation_df_wide)[2:ncol(annotation_df_wide)] <- paste0(colnames(annotation_df_wide)[2:ncol(annotation_df_wide)], "_named")
vec_with_numbers <- as.numeric(stringr::str_extract(colnames(annotation_df_wide), "[0-9]+"))
names(vec_with_numbers) <- colnames(annotation_df_wide)
sorted_colnames <- names(sort(vec_with_numbers, na.last = FALSE))
annotation_df_wide <- annotation_df_wide[, sorted_colnames]
annotation_df_wide <- annotation_df_wide[match(full_seurat@meta.data$Cell_ID, annotation_df_wide$Cell_ID),]

# Clean descriptive markers
descriptive_markers_df <- annotation_results$descriptive_markers_df
descriptive_markers_df <- dplyr::left_join(descriptive_markers_df, annotation_df[, c("cluster_id", "clean_names", "clean_names_withID")], by = c("cluster_id" = "cluster_id"))

##########
### Save annotation
##########

# Save annotation results
data.table::fwrite(annotation_df, file = paste0(parameter_list$harmonization_folder_path,
                                               parameter_list$new_name_suffix, "_", parameter_list$marker_suffix,
                                               "_annotation_result.txt"), sep = "\t")

data.table::fwrite(annotation_df_wide, file = paste0(parameter_list$harmonization_folder_path,
                                                    parameter_list$new_name_suffix, "_", parameter_list$marker_suffix,
                                                    "_annotation_labelmat.txt"), sep = "\t")

data.table::fwrite(descriptive_markers_df, file = paste0(parameter_list$harmonization_folder_path,
                                                        parameter_list$new_name_suffix, "_",
                                                        parameter_list$marker_suffix,
                                                        "_annotation_markers_filtered.txt"), sep = "\t")

#== Annotate manually----------------------------------
# Prepare annotation data for manual annotation
annotation_orig <- annotation_df_wide
rownames(annotation_orig) <- annotation_orig$Cell_ID
annotation_orig <- annotation_orig[-1]
annotation_orig <- annotation_orig %>%
  mutate_each(function(x) stringr::str_replace(x, "C[:digit:]+-[:digit:]+: ", ""))

# Add metadata to Seurat object
full_seurat <- AddMetaData(full_seurat, annotation_orig)
DimPlot(full_seurat, group.by = "C2_named")

# Save annotation data
write.csv(annotation_orig, "result/4.26_annotation/4.26_annotation_orig.csv", row.names = TRUE, quote = FALSE)
full_seurat@misc$annotation_result <- annotation_df
full_seurat@misc$clustering_edgelist <- edgelist

#== Update score-------------------------
# Update marker scores
marker_update <- descriptive_markers_df %>%
  group_by(cluster_id) %>%
  top_n(1, wt = avg_score) %>%
  select(cluster_id, gene)
marker_update <- marker_update %>% column_to_rownames("cluster_id")

# Bind annotations and save
annoBindDf <- anno_bind(full_seurat, c("C2", "C6", "C12", "C21", "C41"), group_select)
annoBindDf <- merge(annoBindDf, marker_update, by = "row.names")
write.csv(annoBindDf, paste0(parameter_list$harmonization_folder_path,
                            parameter_list$new_name_suffix, "_", parameter_list$marker_suffix,
                            "_anno_bind_df.csv"))

# Save the updated Seurat object
saveRDS(full_seurat, "processed_data/integrated_data/20241106_mesenchyme.Rds")
```
