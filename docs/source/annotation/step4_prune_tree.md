# Step 4: Prune Tree

## Introduction
This step involves pruning the hierarchical clustering tree generated by MRTree in previous steps. The pruning process aims to refine the cluster structure by merging small or poorly differentiated clusters based on marker gene expression and cell counts. This helps to create a more robust and biologically meaningful cluster hierarchy.

## Input
- `20241025_mes_cluster_object.Rds`: The MRTree clustering results from [Step 2](./step2_mtree.md).
- `20241024_mesenchyme.Rds`: The Seurat object containing integrated mesenchymal cell data.
- `mesenchyme_parameter.json`: Configuration file with parameters for tree pruning.
- `annotation_mid_file/[new_name_suffix]_[start_node]_markers_siblings_[marker_suffix].tsv`: Marker genes identified by comparing sibling clusters.
- `annotation_mid_file/[new_name_suffix]_[start_node]_markers_all_[marker_suffix].tsv`: Marker genes identified by comparing each cluster to all other clusters.

## Output
- `[new_name_suffix]_pruned_mrtree_clustering_results.rds`: A list containing the pruned clustering results, including updated labelmat, edgelist, nodelist, and data tree structure.

```R
## Script

##########
### Load parameters and packages
##########

# Clear the workspace
rm(list = ls())

# Load utility functions
source("script/utils/anno_integrate.R")
message(Sys.time(),": Starting mrtree pruning .." )

# Load required packages
require(tidyverse)
require(Seurat)
require(Matrix)

# Define the directory for processed data
processed_dir <- "process/annotation/mes_annotation/"

# Load the MRTree clustering results
mrtree_result <- readRDS("process/annotation/mes_annotation/annotation_mid_file/20241025_mes_cluster_object.Rds")

# Load JSON parameters
source("script/utils/anno_utils.R")
load_json("script/annotation/config/mesenchyme_parameter.json")

# Load the Seurat object containing integrated mesenchymal cell data
full_seurat <- readRDS("processed_data/integrated_data/20241024_mesenchyme.Rds")

# Load marker genes identified by comparing sibling clusters
markers_comparisons_siblings <- read.table(paste0(parameter_list$harmonization_folder_path,"annotation_mid_file/",
                                                  parameter_list$new_name_suffix,"_",parameter_list$start_node,"_markers_siblings_",
                                                  parameter_list$marker_suffix,".tsv"),sep="\t",header = T)

# Load marker genes identified by comparing each cluster to all other clusters
markers_comparisons_all <-  read.table(paste0(parameter_list$harmonization_folder_path,"annotation_mid_file/",
                                              parameter_list$new_name_suffix,"_",parameter_list$start_node,
                                              "_markers_all_",parameter_list$marker_suffix,".tsv"),sep="\t",header = T)

##########
### Further clean up marker genes to remove batch expressed genes
##########

if(parameter_list$run_additional_gene_removal){
  # Extract gene expression data
  gene_expr_dataset = full_seurat@assays[[assay_markers]]@data
  max_cells = 30000
  
  # Sample cells if the dataset is too large
  if(ncol(gene_expr_dataset) > max_cells){
    set.seed(1234)
    idx = sample(1:ncol(gene_expr_dataset),max_cells)
    gene_expr_dataset = gene_expr_dataset[,idx]
    dataset_factor = full_seurat@meta.data$Project[idx]
  } else {
    dataset_factor = full_seurat@meta.data$Project
  }
  
  # Binarize gene expression data
  gene_expr_dataset[gene_expr_dataset != 0] = 1
  gene_expr_dataset = as.matrix(gene_expr_dataset)

  # Calculate per dataset occurrence
  per_dataset_occ = data.frame(t(rowsum(t(gene_expr_dataset),dataset_factor)))
  per_dataset_occ2 = data.frame(t(apply(per_dataset_occ,1,function(x){return(round((x/sum(x)),5))})))
  per_dataset_occ2 = apply(per_dataset_occ2,2,function(x){x[is.nan(x)] = 0; return(x)})

  # Filter genes that are not occurring in at least 5 datasets
  thresh = 0.000001
  per_dataset_occ_binary = as.matrix(per_dataset_occ2)
  per_dataset_occ_binary[per_dataset_occ_binary > thresh] = 1
  per_dataset_occ_binary[per_dataset_occ_binary < 1] = 0

  # Calculate entropy
  per_dataset_occ_entropy = apply(per_dataset_occ2,1,entropy_fun)
  per_dataset_occ_min = data.frame(occ=rowSums(per_dataset_occ_binary),entropy = per_dataset_occ_entropy,gene = rownames(per_dataset_occ_binary))

  # Exclude genes based on occurrence and entropy
  balance_exclude_genes = per_dataset_occ_min$gene[per_dataset_occ_min$occ <= 4 | per_dataset_occ_min$entropy <= 1.25 ]
  message("Excluding additional: ",length(balance_exclude_genes)," genes")

  # Update marker comparisons by removing excluded genes
  markers_comparisons_all = markers_comparisons_all[!markers_comparisons_all$gene %in% balance_exclude_genes,]
  markers_comparisons_siblings = markers_comparisons_siblings[!markers_comparisons_siblings$gene %in% balance_exclude_genes,]

  # Save the list of excluded genes
  writeList_to_JSON(list(exlude_genes_pruning = balance_exclude_genes),filename = paste0(parameter_list$harmonization_folder_path,"annotation_mid_file/",parameter_list$new_name_suffix,parameter_list$marker_suffix,"_additionally_removed_markers.json"))
}

##########
### Prepare raw edgelists
##########

# Extract key objects from the MRTree result
edgelist = mrtree_result$edgelist
labelmat = mrtree_result$labelmat

# Add label matrix to Seurat metadata
full_seurat@meta.data = cbind(full_seurat@meta.data,labelmat)

##########
### Run mrtree pruning based on markers
##########

message(Sys.time(),": Prune clusters .." )

# Initialize edgelist and all_nodes
edgelist = edgelist[,c("from","to","level")]
cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
  dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
all_nodes = unique(edgelist$to)

# Remove columns from Seurat metadata that start with "C"
full_seurat@meta.data <- full_seurat@meta.data[, !grepl("^C", names(full_seurat@meta.data))]

# Iterate over all nodes to perform pruning
merge_list = list()
for(n in 1:length(all_nodes)){
  current_node = all_nodes[n]
  parent_node = edgelist$from[edgelist$to == current_node]
  sibling_nodes = edgelist$to[edgelist$from == parent_node & edgelist$to != current_node]
  current_level = edgelist$clusterlevel[edgelist$to == current_node]
  
  # Skip nodes below the minimum prune level
  if(edgelist$level[edgelist$to == current_node] < parameter_list$min_prune_level){next}

  # Filter siblings by cell count
  sibling_nodes = sibling_nodes[sibling_nodes %in% edgelist$to[edgelist$ncells >= parameter_list$min_cells]]
  if(length(sibling_nodes) == 0){next}
  
  # Get cell count of the current node
  ncells_current_node = edgelist$ncells[edgelist$to == current_node]

  # Get number of sibling markers
  sibling_markers = markers_comparisons_siblings %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>%
    dplyr::filter(p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)

  # Get global markers
  global_markers = markers_comparisons_all %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>%
    dplyr::filter(p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)

  # Get parent global markers
  parent_global_markers = markers_comparisons_all %>% dplyr::filter(cluster_id == parent_node) %>% dplyr::arrange(desc(specificity)) %>%
    dplyr::filter(p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)

  # Check if the number of sibling markers or cell count is below the threshold
  if(nrow(sibling_markers) < parameter_list$min_sibling_markers | ncells_current_node < parameter_list$min_cells){
    if(nrow(parent_global_markers) > 3){
      # Get the expression of the parent markers in the current cluster
      current_gene_expression = FetchData(full_seurat,vars = parent_global_markers$gene,cells = full_seurat@meta.data$Cell_ID[full_seurat@meta.data[,current_level] == current_node])
      current_gene_expression_mean = colMeans(current_gene_expression)

      # Iterate over all sibling markers to find the best merge candidate
      intersection_lengths = c()
      coexpr = c()
      for(sib in sibling_nodes){
        sibling_gene_expression = FetchData(full_seurat,vars = parent_global_markers$gene,cells = full_seurat@meta.data$Cell_ID[full_seurat@meta.data[,current_level] == sib])
        sibling_gene_expression_mean = colMeans(sibling_gene_expression)
        coexpr[sib] = cor(current_gene_expression_mean,sibling_gene_expression_mean)
      }

      # Select the node with the highest correlation for merging
      merge_nodes = names(coexpr)[coexpr == max(coexpr)]
    } else {
      merge_nodes = sibling_nodes[1] # Fallback to the first sibling node
    }
    message("Merging node ",current_node," (",n,"/",length(all_nodes),") into node(s) ",paste0(merge_nodes,sep=", "))
    all_nodes_to_merge = c(merge_nodes,current_node)
    merge_list[[paste0(all_nodes_to_merge,collapse = "_")]] = all_nodes_to_merge
  }
}

##########
### Create pruned edgelist
##########

message(Sys.time(),": Update labels and save .." )

# Initialize updated edgelist and labelmat
edgelist_updated = edgelist
labelmat_updated = labelmat
merge_list2 = merge_list

# Update edgelist and labelmat based on merge_list
for(i in 1:length(merge_list2)){
  nodes_to_merge = merge_list2[[i]]
  print(nodes_to_merge)
  
  # Only update if there are truly 2 unique labels remaining
  if(length(unique(nodes_to_merge)) > 1){
    merge_node = as.character(sort((nodes_to_merge))[1])
    print(paste0(" >> merge to ",merge_node))
    
    # Update edgelist
    edgelist_updated$from[edgelist_updated$from %in% nodes_to_merge] = merge_node
    edgelist_updated$to[edgelist_updated$to %in% nodes_to_merge] = merge_node
    edgelist_updated = edgelist_updated %>% distinct(from,to,clusterlevel) # Merge duplicate rows to one entry and update the cell total
    
    # Update labelmat
    current_level = edgelist$clusterlevel[edgelist$to == merge_node]
    labelmat_updated[labelmat_updated[,current_level] %in% nodes_to_merge,current_level] = merge_node
    
    # Change all occurrences in merge_list
    merge_list2 = lapply(merge_list2,function(x,remove_nodes,new_node){x[x %in% remove_nodes] = new_node;return(x)},remove_nodes=nodes_to_merge,new_node=merge_node)
  }
}

##########
### Update labels in edgelist and labelmat
##########

# Define old and new prefixes for cluster levels
old_prefix = parameter_list$old_prefix # "K"
new_prefix = parameter_list$new_prefix # "C"

# Update cluster levels in the edgelist
all_cluster_levels_updated = edgelist_updated %>% dplyr::group_by(clusterlevel) %>% dplyr::count()
all_cluster_levels_updated$clusterlevel_new = paste0(all_cluster_levels_updated$clusterlevel %>% stringr::str_extract(old_prefix)  %>% stringr::str_replace(old_prefix, new_prefix),all_cluster_levels_updated$n)

# Create new labels using the pruned numbers of labels per level
all_rows = list()
for(c in 1:nrow(all_cluster_levels_updated)){
  edgelist_updated_current_level = edgelist_updated[edgelist_updated$clusterlevel == all_cluster_levels_updated$clusterlevel[c],]
  for(i in 1:nrow(edgelist_updated_current_level)){
    new_node_label = data.frame(old_node = edgelist_updated_current_level$to[i],
                                new_node = paste0(all_cluster_levels_updated$clusterlevel_new[c],"-",i),
                                new_cluster_level = all_cluster_levels_updated$clusterlevel_new[c])
    all_rows[[new_node_label$new_node]] = new_node_label
  }
}
new_labels = as.data.frame(do.call(rbind,all_rows))

# Make new edgelist with updated labels
edgelist_updated_new_labels = dplyr::left_join(edgelist_updated,new_labels[,c("old_node","new_node")],by=c("from"="old_node"))
edgelist_updated_new_labels = dplyr::left_join(edgelist_updated_new_labels,new_labels,by=c("to"="old_node"))
edgelist_updated_new_labels = edgelist_updated_new_labels %>% dplyr::select(from = new_node.x, to = new_node.y,clusterlevel = new_cluster_level)
edgelist_updated_new_labels$from[is.na(edgelist_updated_new_labels$from)] = "all"

# Make new labelmat with updated labels
labelmat_updated_new_labels = apply(labelmat_updated,2,function(x,new_labels){
  new_labels$new_node[match(x,new_labels$old_node)]
},new_labels=new_labels)

# Update column names of the new labelmat
colnames(labelmat_updated_new_labels) = stringr::str_extract(labelmat_updated_new_labels[1,],pattern = paste0(parameter_list$new_prefix,"[0-9]+"))

##########
### Create pruned result version
##########

# Save the pruned result in data.tree format
require(data.tree)
df = as.data.frame(unique(labelmat_updated_new_labels), stringsAsFactors = F)
df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
tree.datatree = data.tree::as.Node(df)

# Export the edgelist
edges = data.tree::ToDataFrameNetwork(tree.datatree,"isLeaf","level","count","totalCount","height")
nodes = data.frame(id = c("all",as.character(unique(edges$to))),label=c("all",as.character(unique(edges$to))))
nodes = rbind(c("all",FALSE,1,2,208,max(edges$height)+1),edges[,2:ncol(edges)]) %>% dplyr::rename(id = to) %>% dplyr::mutate(label = id)

# Create the updated cluster object and add to misc
updated_cluster_object = list(labelmat = labelmat_updated_new_labels,
                              edgelist = edges,
                              nodelist = nodes,
                              data_tree = tree.datatree)

##########
### Save
##########

# Save the updated cluster object
saveRDS(updated_cluster_object,paste0(parameter_list$harmonization_folder_path,"annotation_mid_file/",parameter_list$new_name_suffix,"_pruned_mrtree_clustering_results",".rds"))

message("Finalized pruning")
```