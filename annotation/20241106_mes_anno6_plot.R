###############################################################################
# Tooth Mesenchyme Annotation and Visualization
# Script: Cluster Tree and Heatmap Visualization for Annotated Cell Types
#
# This script generates cluster tree visualizations with heatmaps showing:
# - Age distribution
# - Histology distribution
# - Cell count information
#
# Author: Modified by Gilbert
# Date: 2024-03-26
###############################################################################

# Load required libraries
library(aplot)
library(Seurat) # For Seurat object manipulation
library(dplyr)
library(ggplot2)
library(forcats) # For fct_relevel
source("script/utils/fun_plot_figure.R")
source("script/utils/fun_plotclustertree.R")

# Create results directory
results_path_figure2 <- "results/annotation/20241106_anno_step6_plot"
dir.create(results_path_figure2, showWarnings = FALSE, recursive = TRUE)

# Load color settings
load_colors()

# Define visualization parameters
leaf_level_column <- "C57"
leaf_level <- 6

###############################################################################
# Data Preparation: Filter Clusters and Process Edge List
###############################################################################

# Filter out lowquality clusters from the edge list
exclusion_clusters <- c("C9-5", "C9-9", "C9-3", "C9-4", "C19-9",
                         "C19-10", "C19-11", "C19-14",
                         "C19-18", "C19-19", "C29-14", "C29-15", "C29-22",
                         "C29-23", "C29-28", "C29-29", "C57-29", "C57-30",
                         "C57-31", "C57-39", "C57-40", "C57-57", "C92-92",
                         "C92-64", "C92-65", "C92-66", "C92-67", "C92-68", "C57-56",
                         "C92-49", "C92-90", "C92-91", "C92-48")

# Set cell identities
Idents(full_seurat) <- full_seurat$C19

# Backup the original edge list
edgelist_bk <- full_seurat@misc$clustering_edgelist

# Filter the edge list to remove excluded clusters
edgelist <- full_seurat@misc$clustering_edgelist
edgelist <- edgelist %>%
  dplyr::filter(!(from %in% exclusion_clusters | to %in% exclusion_clusters))

# Filter the annotation dataframe similarly
annotation_df <- annotation_df %>%
  dplyr::filter(!(cluster_id %in% exclusion_clusters))

###############################################################################
# Histology Annotation
###############################################################################

# Read metadata and prepare histology lookup
metadata <- read.csv("data/metadata/20250326_metadata.csv", row.names = 1)

# Filter metadata based on full_seurat samples
metadata <- metadata[rownames(metadata) %in% full_seurat$Sample, ]

# Create a unique dataframe for Project and Histology
metaAnno <- unique(metadata[c("Project", "Histology")])

# Create a lookup vector for Histology based on Project
histology_lookup <- setNames(metaAnno$Histology, metaAnno$Project)

# Use match to ensure the order is preserved
indices <- match(full_seurat$Project, names(histology_lookup))
annotated_histology <- histology_lookup[indices]

# Assign the annotated histology to full_seurat
full_seurat$Histology <- as.character(annotated_histology)

###############################################################################
# Heatmap Matrix Generation Functions
###############################################################################

# Function to create a percentage heatmap matrix for a given grouping variable
create_percent_heatmap <- function(seurat_obj, grouping_var, cluster_column) {
  # Create contingency data
  heatmap_data <- seurat_obj@meta.data %>%
    dplyr::select(Cell_ID, !!sym(grouping_var), !!sym(cluster_column)) %>%
    dplyr::group_by(!!sym(cluster_column), !!sym(grouping_var)) %>%
    dplyr::count(name = "presence") %>%
    dplyr::group_by(!!sym(cluster_column)) %>%
    dplyr::mutate(presence = presence / sum(presence) * 100) %>%
    dplyr::ungroup() %>%
    tidyr::spread(key = !!sym(grouping_var), value = presence) %>%
    as.data.frame()

  # Convert to matrix
  matrix_data <- as.matrix(heatmap_data[, 2:ncol(heatmap_data)])
  rownames(matrix_data) <- heatmap_data[[cluster_column]]
  matrix_data[is.na(matrix_data)] <- 0

  return(matrix_data)
}

# Function to create a cell count percentage matrix
create_count_heatmap <- function(seurat_obj, cluster_column) {
  heatmap_data <- seurat_obj@meta.data %>%
    dplyr::select(Cell_ID, !!sym(cluster_column)) %>%
    dplyr::group_by(!!sym(cluster_column)) %>%
    dplyr::count(name = "ncells") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ncells_pct = ncells / sum(ncells) * 100) %>%
    as.data.frame()

  matrix_data <- as.matrix(heatmap_data[, "ncells_pct", drop = FALSE])
  colnames(matrix_data) <- "%"
  rownames(matrix_data) <- heatmap_data[[cluster_column]]

  return(matrix_data)
}

###############################################################################
# Generate Heatmap Matrices
###############################################################################

# Factor Age variable for proper ordering
full_seurat$Age <- factor(full_seurat$Age)

# Create Age-based heatmap matrix
heatmap_matrix <- create_percent_heatmap(full_seurat, "Age", leaf_level_column)

# Get proper order for displaying projects by age
newOrder <- full_seurat@meta.data %>%
  dplyr::select(Project, Age) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq != 0) %>%
  arrange(Age)

# Create column name mapping for reference
colnames_overview <- data.frame(
  number = 1:ncol(heatmap_matrix),
  Dataset = colnames(heatmap_matrix)
)
data.table::fwrite(colnames_overview,
                  paste0(results_path_figure2, "/circular_tree_colnames.txt"),
                  sep = "\t")

# Store original column names before replacing (we'll use these for display)
original_colnames <- colnames(heatmap_matrix)

# Create cell count heatmap matrix
heatmap_matrix2 <- create_count_heatmap(full_seurat, leaf_level_column)

# Create histology-based heatmap matrix
full_seurat$Histology <- factor(full_seurat$Histology)
heatmap_matrix3 <- create_percent_heatmap(full_seurat, "Histology", leaf_level_column)

# Get proper order for projects by histology
newOrder2 <- full_seurat@meta.data %>%
  dplyr::select(Project, Histology) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq != 0) %>%
  arrange(Histology)

# Create column name mapping for reference
colnames_overview3 <- data.frame(
  number = 1:ncol(heatmap_matrix3),
  Dataset = colnames(heatmap_matrix3)
)
data.table::fwrite(colnames_overview3,
                  paste0(results_path_figure2, "/circular_tree_colnames_histology.txt"),
                  sep = "\t")

# Store original histology column names
original_colnames3 <- colnames(heatmap_matrix3)

###############################################################################
# Prepare Annotation Data
###############################################################################

# Create annotation dataframe for cluster tree
anno_df <- annotation_df %>%
  dplyr::select(cluster_id, clusterlevel, cluster_name = clean_names)

# Add first part of cluster name for labeling
anno_df$first_cluster_name <- sapply(anno_df$cluster_name, function(x) {
  strsplit(x, "\\.")[[1]][1]
})

###############################################################################
# Plot Cluster Tree
###############################################################################

# Set tree appearance parameters
tree_color <- "grey70"

# Create base cluster tree
circular_tree <- plot_cluster_tree(
  edgelist = edgelist,
  leaf_level = leaf_level,
  anno_df = anno_df,
  metadata = full_seurat@meta.data,
  label_size = 4,
  show_genes = TRUE,
  vjust_label = -0.5,
  edge_color = tree_color,
  node_color = tree_color
)

# Rotate tree for better visualization
circular_tree <- rotate_tree(circular_tree, -90)

# Save the basic tree
ggsave(paste0(results_path_figure2, "/20250220_circletree_full_layer.pdf"),
       plot = circular_tree, width = 8, height = 8)

###############################################################################
# Add Heatmaps to Cluster Tree
###############################################################################

# Add Age heatmap
circular_tree_heat <- add_heatmap(
  circular_tree = circular_tree,
  heatmap_matrix = heatmap_matrix,
  heatmap_colors = c(bg_col, "darkred"),
  heatmap_colnames = TRUE,
  legend_title = "Pct Age",
  matrix_offset = 2.2,
  matrix_width = 0.8,
  colnames_angle = 0,
  legend_text_size =2,
  hjust_colnames = 0.5,   # Adjusted for better alignment
  na_color = "white",
  heatmap_text_size = 1   # Increased from 0 to make text visible
)

# Add Histology heatmap
circular_tree_heat <- add_heatmap(
  circular_tree = circular_tree_heat,
  heatmap_matrix = heatmap_matrix3,
  heatmap_colors = c(bg_col, "darkgreen"),
  heatmap_colnames = TRUE,
  legend_title = "Pct Tissue",
  matrix_offset =12,
  matrix_width = 0.4,
  colnames_angle = 0,
  legend_text_size = 2,
  hjust_colnames = 0.5,   # Adjusted for better alignment
  na_color = "white",
  heatmap_text_size = 1   # Increased from 0 to make text visible
)


# Save the circular tree with heatmaps
ggsave(paste0(results_path_figure2, "/20250220_circletree_full_heatmaps.pdf"),
       plot = circular_tree_heat, width = 12, height = 12)

###############################################################################
# Create Linear Tree with Heatmaps
###############################################################################

# Create linear tree
normaltree <- ggtree(circular_tree$data) +
  geom_nodelab(aes(x = branch, label = first_cluster_name),
               color = "darkred", vjust = -0.25, size = 2.5) +
  geom_tiplab(ggplot2::aes(x = branch, label = first_cluster_name),
              color = "darkred", vjust = -0.25, size = 2.5)

# Add Age heatmap to linear tree
normal_tree_heat <- add_heatmap(
  circular_tree = normaltree,
  heatmap_matrix = heatmap_matrix,
  heatmap_colors = c(bg_col, "darkred"),
  heatmap_colnames = TRUE,
  legend_title = "Pct Age",
  matrix_offset = 0.5,
  matrix_width = 0.57,
  colnames_angle = 0,
  legend_text_size = 8,
  hjust_colnames = 0,
  na_color = "white",
  scale_limits = c(0, 40),
  heatmap_text_size = 2  # Increased from 0 to make text visible
)

# Add Histology heatmap to linear tree
normal_tree_heat <- add_heatmap(
  circular_tree = normal_tree_heat,
  heatmap_matrix = heatmap_matrix3,
  heatmap_colors = c(bg_col, "darkorange"),
  heatmap_colnames = TRUE,
  legend_title = "Pct Tissue",
  matrix_offset = 3.35,
  matrix_width = 0.2,
  colnames_angle = 0,
  legend_text_size = 8,
  hjust_colnames = 0,
  na_color = "white",
  heatmap_text_size = 2  # Increased from 0 to make text visible
)

# Display linear tree with heatmaps
normal_tree_heat

# Save the linear tree with heatmaps
ggsave(paste0(results_path_figure2, "/20250220_lineartree_heatmaps.pdf"),
       plot = normal_tree_heat, width = 8, height = 6)

###############################################################################
# Generate DotPlot of Marker Genes
###############################################################################

# Set identities to the desired clusters
Idents(full_seurat) <- full_seurat@meta.data[[leaf_level_column]]
selectCluster <- levels(full_seurat)

# Sort clusters numerically
index <- order(as.numeric(sub("C57-", "", selectCluster)))
selectCluster <- selectCluster[index]

# Get the order of clusters based on the tree
clusterOrder <- normaltree$data %>%
  dplyr::filter(clusterlevel == "C57") %>%
  arrange(y) %>%
  select(label) %>%
  unlist()

# Get top marker genes for each cluster
comparisons_all_update <- comparisons_all
comparisonGene <- comparisons_all_update %>%
  dplyr::filter(cluster_id %in% selectCluster & !(gene %in% all_exclusion_genes)) %>%
  group_by(cluster_id) %>%
  arrange(desc(specificity)) %>%
  top_n(1) %>%
  arrange(fct_relevel(cluster_id, clusterOrder)) %>%
  ungroup() %>%
  select(gene) %>%
  unlist()
clusterOrder <- c(clusterOrder,setdiff(levels(full_seurat),as.character(clusterOrder)))
# Order the levels based on tree structure
levels(full_seurat) <- as.character(clusterOrder)

# Create dotplot of marker gene expression
p2 <- DotPlot(
  full_seurat,
  features = unique(comparisonGene),
  assay = "originalexp",
  cluster.idents = TRUE,
  cols = c("lightgrey", "tan1")
) +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 70, size = 10, face = "bold", hjust = 1)
  )

# Combine tree heatmap with dotplot
p3 <- p2 %>%
  insert_left(normal_tree_heat, width = 0.8)

# Display combined plot
p3

# Save the final combined visualization
ggsave(paste0(results_path_figure2, "/20250220_tree_dotplot_combined.pdf"),
       plot = p3, height = 8, width = 12)

# Save key objects for future use
#saveRDS(p3, paste0(results_path_figure2, "/20250220_tree_dotplot_combined.Rds"))
#saveRDS(circular_tree_heat, paste0(results_path_figure2, "/20250220_circular_heatmap.Rds"))
saveRDS(full_seurat, "processed_data/integrated_data/20250326_mesenchyme.Rds")

# Print completion message
cat("Visualization script completed successfully.\n")
