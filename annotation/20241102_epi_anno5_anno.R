#== Change age label---------------------------------
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

#== find marker--------------------------

all_markers_mrtree_list_update = findMarkers_tree2(full_seurat,
                                                   edgelist=edgelist_updated_new_labels[,1:2],
                                                   labelmat = labelmat_updated_new_labels,
                                                   n_cores = 15,
                                                   assay=assay_markers,
                                                   slot=assay_slot,
                                                   test.use = test.use,
                                                   batch_var=batch_var,
                                                   latent_vars_seurat = NULL,
                                                   genes_to_include = genes_to_include,
                                                   logfc.threshold = logfc.threshold,
                                                   min.pct = min.pct,
                                                   min.diff.pct = min.diff.pct,
                                                   max.cells.per.ident=max.cells.per.ident,
                                                   min.cells.feature = min.cells.feature,
                                                   min.cells.group = min.cells.group,
                                                   base = base,
                                                   only.pos = only.pos)

comparisons_all_update <- all_markers_mrtree_list_update$comparisons_All
comparisons_siblings_update <- all_markers_mrtree_list_update$comparisons_Sibling
comparisons_all_update$specificity = ((comparisons_all_update$pct.1 + specificity_base) /
                                        (comparisons_all_update$pct.2 + specificity_base)) * comparisons_all_update$avg_log2FC
comparisons_siblings_update$specificity = ((comparisons_siblings_update$pct.1 + specificity_base) /
                                             (comparisons_siblings_update$pct.2 + specificity_base)) * comparisons_siblings_update$avg_log2FC

View(comparisons_siblings_update)
write.csv(comparisons_all_update,paste0(process_dir,"comparisons_all_update.csv"))
write.csv(comparisons_siblings_update,paste0(process_dir,"comparisons_siblings_update.csv"))


#== label--------------------------------
# full_seurat$C2 <- labelmat_updated_new_labels[,1]
# full_seurat$C9 <- labelmat_updated_new_labels[,2]
# full_seurat$C19 <- labelmat_updated_new_labels[,3]
# full_seurat$C29 <- labelmat_updated_new_labels[,4]
# full_seurat$C57 <- labelmat_updated_new_labels[,5]
# full_seurat$C92 <- labelmat_updated_new_labels[,6]
# #full_seurat$C137 <- labelmat_updated_new_labels[,7]
# DimPlot(full_seurat,group.by = "C2",label = T)
# DimPlot(full_seurat,group.by = "C9",label = T)
# DimPlot(full_seurat,group.by = "C19",label = T)
# DimPlot(full_seurat,group.by = "C29",label = T)
# DimPlot(full_seurat,group.by = "C57",label = T)
# DimPlot(full_seurat,group.by = "C92",label = T)
# #p1+p2
# #DimPlot(full_seurat)


# identify indices of elements that meet the condition
idx <- which(C36_pred > 0.4, arr.ind = TRUE)

# extract row and column names
rows <- rownames(C36_pred)[idx[, 1]]
cols <- colnames(C36_pred)[idx[, 2]]

# print row and column names
for (i in 1:length(rows)) {
  cat("Value", C36_pred[rows[i], cols[i]], "in", cols[i], "column and", rows[i], "row exceeds the threshold of", 0.5, "\n")
}

# identify indices of elements that meet the condition
idx <- which(C36_pred > 0.4, arr.ind = TRUE)

# extract row and column names
rows <- rownames(C36_pred)[idx[, 1]]
cols <- colnames(C36_pred)[idx[, 2]]

# print row and column names
for (i in 1:length(rows)) {
  cat("Value", C36_pred[rows[i], cols[i]], "in", cols[i], "column and", rows[i], "row exceeds the threshold of", 0.5, "\n")
}


prediction_large <- multilevel_predict(full_seurat,"X_scANVI",level = train_level,pred_level = "C2",subsample_size = 1500)
prediction_large <- prediction_large[, -which(names(prediction_large) == "ident")]
full_seurat$merge_id_level3
group_select <- c("Organ", "Tissue", "Tissue.Specific.", "Stage","Age",
                  "Age.In.Detail.","Species","Origin", "coarse_label",train_level)

annoBindDf <- anno_bind(full_seurat,c("C2","C7","C19","C36","C49","C90","C137"),group_select)
Idents(full_seurat) <- full_seurat$C19
cell_19_7 <- WhichCells(full_seurat,idents = "C19-5")
DimPlot(full_seurat,cells.highlight  = cell_19_7)
full_seurat$merge_id_level3
cell_late_mes <- colnames(full_seurat)[full_seurat$merge_id_level3=="C_Sox18+ Dermis"]
DimPlot(full_seurat,cells.highlight  = cell_late_mes)

cell_test <- WhichCells(full_seurat,idents = "LA_Tspan11.Dpt.Fibroblast")
DimPlot(full_seurat,cells.highlight  = cell_test)

percentPlotFun(full_seurat, "C36", "Organ")

Idents(full_seurat) <- full_seurat$C90
cell_90_33 <- WhichCells(full_seurat,idents = "C90-33")
DimPlot(full_seurat,cells.highlight  = cell_90_33)

Idents(full_seurat) <- full_seurat$anno_level_5
cell_test <- WhichCells(full_seurat,idents = "LE_Aldh1a3.Hypertrophic.Chondrocyte")
DimPlot(full_seurat,cells.highlight  = cell_test)
percentPlotFun(full_seurat,cluster_id = "C6",group_id = "Tissue")
DimPlot(full_seurat,group.by = "K65",label = T)
FeaturePlot(full_seurat,"Col10a1")
FeaturePlot(full_seurat,c("Prg4","Igfbp6"))


df <- table(full_seurat$C6,full_seurat$Age.In.Detail.)%>%t%>%as.matrix()

df <- t(apply(df, 1, function(x) x/sum(x)))
df <- as.data.frame(df)

df_long <- df%>%rownames_to_column("Age")%>%
  pivot_longer(cols = -Age,names_to = "Cluster")


#
# # test on full
# # parameter_list = jsonlite::read_json("data/parameters_annotation_v2_1.json")
# # parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
#
# # read features to excludes
features_exclude_list= balance_exclude_genes
# # features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))
# #features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
#
# # load seurat
# full_seurat = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))
#
# # load mrtree clustering
# mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))
#
# # mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_mrtree_clustering_results",".rds"))
#
# message("All markers for: ",length(unique(markers_comparisons_all$cluster_id))," clusters available")
# message("Sibling markers for: ",length(unique(markers_comparisons_siblings$cluster_id))," clusters available")
#
#
# # parameter_list$marker_suffix,
# additional_remove_genes = unlist(jsonlite::read_json(unlist(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"raw","_additionally_removed_markers.json"))))
# # gather some other genes that are not informative during annoation:
other_genes_remove = rownames(full_seurat@assays$originalexp@counts)[grepl("RP|Gm|Hb[ab]|Rik|-ps",rownames(full_seurat@assays$originalexp@counts))]
# make list of all genes that should be removed:
all_exclusion_genes = unique(c(features_exclude_list,other_genes_remove))



##########
### Load parameters and packages
##########

mrtree_result <- updated_cluster_object
# get key objects
edgelist = mrtree_result$edgelist
labelmat = mrtree_result$labelmat

edgelist = edgelist[,c("from","to","level")]
cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
  dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
all_nodes = unique(edgelist$to)

markers_comparisons_siblings <- comparisons_siblings
markers_comparisons_all <- comparisons_all
##########
### Run annotation
##########
command_args<-commandArgs(TRUE)
param_file = "script/annotation/config/parameters_wtintegrate.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
message("Starting annotation: ")
print(parameter_list$manual_names_annotation_abbr)

annotation_results = annotate_tree(edgelist = edgelist,#[1:32,],#edgelist = edgelist[1:291,],
                                   labelmat = labelmat_updated_new_labels,
                                   markers_comparisons_all = comparisons_all_update,
                                   markers_comparisons_siblings = comparisons_siblings_update,
                                   preferred_genes=character(0),
                                   manual_names= parameter_list$manual_names_annotation_abbr,
                                   overwrite_with_manual = TRUE,
                                   manual_exclude_genes=all_exclusion_genes,
                                   max_pval_adj= parameter_list$max_pval_adj,
                                   min_specificity = parameter_list$min_specificity,
                                   min_specificity_sibling_children= parameter_list$min_specificity_sibling_children,
                                   scale_preferred=1,
                                   limit_factor=parameter_list$limit_factor,
                                   max_score_siblings_children = parameter_list$max_score_siblings_children,
                                   reverse_order = parameter_list$reverse_order)

message("Formating annotation results")

## get annotation
annotation_df = annotation_results$annotation_df
# update cluster names with ID
annotation_df$clean_names_withID = paste0(annotation_df$cluster_id,": ",annotation_df$clean_names)
# get labelmat and add cell ids
labelmat_updated = cbind(Cell_ID=full_seurat@meta.data$cellid,mrtree_result$labelmat)
# make a per cell cluster annotation
cell_cluster_map = labelmat_updated %>% as.data.frame() %>% tidyr::gather(-Cell_ID,key="clusterlevel",value="cluster_id")
cell_cluster_map$clusterlevel = gsub("_pruned","",cell_cluster_map$clusterlevel)
# make wide version of labelmat ( that can be added to seurat metadata)
annotation_df_wide = annotation_df %>% dplyr::left_join(cell_cluster_map,by=c("clusterlevel"="clusterlevel","cluster_id"="cluster_id")) %>%
  dplyr::select(Cell_ID,clusterlevel,clean_names_withID)  %>% tidyr::spread(key = clusterlevel,value = clean_names_withID)
colnames(annotation_df_wide)[2:ncol(annotation_df_wide)] = paste0(colnames(annotation_df_wide)[2:ncol(annotation_df_wide)],"_named")
# sort columns in mrtree annotation wide
vec_with_numbers = as.numeric(stringr::str_extract(colnames(annotation_df_wide),"[0-9]+"))
names(vec_with_numbers) = colnames(annotation_df_wide)
sorted_colnames = names(sort(vec_with_numbers,na.last = FALSE))
annotation_df_wide = annotation_df_wide[,sorted_colnames]
# ensure order makes sense
annotation_df_wide = annotation_df_wide[match(full_seurat@meta.data$cellid,annotation_df_wide$Cell_ID),]

# remove existing if running again on same object:
# full_seurat@meta.data = full_seurat@meta.data[,!grepl("C*\\_named",colnames(full_seurat@meta.data))]
# # add to seurat
# full_seurat@meta.data= dplyr::left_join(full_seurat@meta.data,annotation_df_wide,by="Cell_ID")
# rownames(full_seurat@meta.data) =  full_seurat@meta.data$Cell_ID

## clean markers
descriptive_markers_df = annotation_results$descriptive_markers_df
descriptive_markers_df = dplyr::left_join(descriptive_markers_df,annotation_df[,c("cluster_id","clean_names","clean_names_withID")],by=c("cluster_id"="cluster_id"))

##########
### Save annotation
##########

data.table::fwrite(annotation_df,file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_result.txt") ,sep="\t")

data.table::fwrite(annotation_df_wide,file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_labelmat.txt") ,sep="\t")

data.table::fwrite(descriptive_markers_df,file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,"_annotation_markers_filtered.txt") ,sep="\t")


# #== annoated_manual----------------------------------
annotation_orig <- annotation_df_wide
rownames(annotation_orig) <- annotation_orig$Cell_ID
annotation_orig <- annotation_orig[-1]
annotation_orig <- annotation_orig%>%
  mutate_each(function(x) stringr::str_replace(x, "C[:digit:]+-[:digit:]+: ",""))

full_seurat <- AddMetaData(full_seurat,annotation_orig)
DimPlot(full_seurat,group.by = "C2_named")
write.csv(annotation_orig,"result/4.26_annotation/4.26_annotation_orig.csv",row.names = T,quote = F)
full_seurat@misc$annotation_result <- annotation_df
full_seurat@misc$clustering_edgelist <- edgelist

#== update score-------------------------
marker_update <- descriptive_markers_df%>%
  group_by(cluster_id)%>%
  top_n(1,wt = avg_score)%>%
  select(cluster_id,gene)
marker_update <- marker_update%>%column_to_rownames("cluster_id")
annoBindDf <- anno_bind(full_seurat,c("C2","C6","C12","C21","C41"),group_select)
annoBindDf <- merge(annoBindDf,marker_update, by = "row.names")
write.csv(annoBindDf,"result/4.8_annotation/4.9_anno_bind_df.csv")



mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(38)
DimPlot(full_seurat,group.by = "Age.In.Detail.",cols = mycolor)


