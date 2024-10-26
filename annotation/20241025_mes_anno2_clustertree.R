#== set environment------------------
rm(list = ls())
library(mrtree)
library(dplyr)
source("script/utils//mrtree_fun.R")
process_dir = "process/annotation/mes_annotation/annotation_mid_file/"
mrtree_input_labels <- data.table::fread("processed_data/metadata//20241025_mes_anno_mrtree_input_label.csv")

cluster_matrix_for_mrtree = as.matrix(mrtree_input_labels)
#cluster_matrix_for_mrtree <- apply(cluster_matrix_for_mrtree, 2, as.character)
cluster_matrix_for_mrtree <- apply(cluster_matrix_for_mrtree, 2, as.numeric)
#print(head(cluster_matrix_for_mrtree))

message(Sys.time(),": Build mrtree using matrix with: ",dim(cluster_matrix_for_mrtree)[1]," cells and ",dim(cluster_matrix_for_mrtree)[2]," cluster levels." )
# feed into mrtree
# function from 'mrtree_functions.R' which is copied from original repo
mrtree_res <- mrtree(cluster_matrix_for_mrtree,
                     prefix = "leiden_clusters_level_",
                     suffix = NULL,
                     max.k = Inf,
                     consensus = FALSE,
                     sample.weighted = FALSE, # parameter_list$mr_tree_weighted,
                     augment.path = FALSE,
                     verbose = TRUE,
                     n.cores = 5)

saveRDS(mrtree_res,paste0(process_dir,"mes_mrtree.Rds"))

parameter_list = jsonlite::read_json("script/annotation/config/mesenchyme_parameter.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

message(Sys.time(),": Saved raw result of length ",length(mrtree_res)," to: ",paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_mrtree_res_raw",".rds"))

message(Sys.time(),": Create mrtree ouput list" )
# make labelmat with colnames included
if(parameter_list$use_recon_labelmat){
  labelmat=mrtree_res$labelmat.recon
  Ks.recon = apply(labelmat, 2, function(y) length(unique(y)))
  unique.idx = 1:length(Ks.recon) # don#t remove last col
  colnames(labelmat) = paste0("K", Ks.recon[unique.idx])
}else{
  labelmat=mrtree_res$labelmat.mrtree
}

# build dataframe with labels
n=nrow(labelmat)
backup = colnames(labelmat)
labelmat = matrix(paste(matrix(rep(colnames(labelmat),each=n), nrow = n), labelmat, sep='-'), nrow = n)
colnames(labelmat)=backup
df = as.data.frame(unique(labelmat), stringsAsFactors = F)

# save in data.tree format
require(data.tree)
df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
tree.datatree = data.tree::as.Node(df)

# export edgelist  and nodelist from data.tree
edges= data.tree::ToDataFrameNetwork(tree.datatree,"isLeaf","level","count","totalCount","height")
nodes = data.frame(id = c("all",as.character(unique(edges$to))),label=c("all",as.character(unique(edges$to))))
nodes = rbind(c("all",FALSE,1,2,228,max(edges$height)+1),edges[,2:ncol(edges)]) %>% dplyr::rename(id = to) %>% dplyr::mutate(label = id)

# make cluster object (list)
cluster_object = list(labelmat = labelmat,
                      edgelist = edges ,
                      nodelist = nodes,
                      data_tree = tree.datatree,
                      mrtree_output = mrtree_res)
#seurat_object_harmonized@misc$mrtree_clustering_results = cluster_object


labelmat=mrtree_res$labelmat.recon
Ks.recon = apply(labelmat, 2, function(y) length(unique(y)))
unique.idx = 1:length(Ks.recon) # don#t remove last col
colnames(labelmat) = paste0("K", Ks.recon[unique.idx])


n=nrow(labelmat)
backup = colnames(labelmat)
labelmat = matrix(paste(matrix(rep(colnames(labelmat),each=n), nrow = n), labelmat, sep='-'), nrow = n)
colnames(labelmat)=backup
df = as.data.frame(unique(labelmat), stringsAsFactors = F)



require(data.tree)
df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
tree.datatree = data.tree::as.Node(df)

edges= data.tree::ToDataFrameNetwork(tree.datatree,"isLeaf","level","count","totalCount","height")
nodes = data.frame(id = c("all",as.character(unique(edges$to))),label=c("all",as.character(unique(edges$to))))
nodes = rbind(c("all",FALSE,1,5,max(edges$height)+1),edges[,2:ncol(edges)]) %>% dplyr::rename(id = to) %>% dplyr::mutate(label = id)

cluster_object = list(labelmat = labelmat,
                      edgelist = edges ,
                      nodelist = nodes,
                      data_tree = tree.datatree,
                      mrtree_output = mrtree_res)


##########
###Save
##########
saveRDS(cluster_object,paste0(process_dir,"20241025_mes_cluster_object.Rds"))
# # save object as rds
# saveRDS(cluster_object,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))

message("Finalized")
