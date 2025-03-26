library(Seurat)

mes <- readRDS("processed_data/integrated_data/20241106_mesenchyme.Rds")
DimPlot(mes,group.by = "C9_named")
FeaturePlot(mes,"Fibin")
library(tidyverse)


source("script/utils/anno_utils.R")
source("script/utils/stratified_wilcoxon_functions.R")
source("script/utils/anno_integrate.R")
process_dir <- "process/annotation/mes_annotation/annotation_mid_file/"
load_json("script/annotation/config/mesenchyme_parameter.json")
full_seurat <- mes
updated_cluster_object <- readRDS(paste0(parameter_list$harmonization_folder_path,"annotation_mid_file/",
                                         parameter_list$new_name_suffix,"_pruned_mrtree_clustering_results",".rds"))
edgelist_updated_new_labels <- updated_cluster_object$edgelist
labelmat_updated_new_labels <- updated_cluster_object$labelmat
genes_to_include <- read.csv("process/annotation/mes_annotation/annotation_mid_file/20241025_gene_include.csv",row.names = 1) %>% unlist
full_seurat$C2 <- labelmat_updated_new_labels[,1]
DimPlot(full_seurat,group.by = "C2")
#== find marker--------------------------
all_markers_mrtree_list_update = findMarkers_tree2(full_seurat,
                                                   edgelist=edgelist_updated_new_labels[,1:2],
                                                   labelmat = labelmat_updated_new_labels,
                                                   n_cores = 15,
                                                   assay=assay_markers,
                                                   slot=assay_slot,
                                                   test.use = test.use,
                                                   batch_var=batch_var,
                                                   latent_vars_seurat = latent_vars_seurat,
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
write.csv(comparisons_all_update,paste0(process_dir,"comparisons_sibling_update.csv"))

comparisons_siblings <- read.csv(paste0(process_dir,"comparisons_sibling_update.csv"),row.names = 1)
comparisons_all <- read.csv(paste0(process_dir,"comparisons_all_update.csv"),row.names = 1)



#== label--------------------------------
full_seurat$C2 <- labelmat_updated_new_labels[,1]
full_seurat$C9 <- labelmat_updated_new_labels[,2]
full_seurat$C19 <- labelmat_updated_new_labels[,3]
full_seurat$C29 <- labelmat_updated_new_labels[,4]
full_seurat$C57 <- labelmat_updated_new_labels[,5]
full_seurat$C92 <- labelmat_updated_new_labels[,6]
#full_seurat$C137 <- labelmat_updated_new_labels[,7]
DimPlot(full_seurat,group.by = "C2",label = T)
DimPlot(full_seurat,group.by = "C9",label = T)
DimPlot(full_seurat,group.by = "C19",label = T)
DimPlot(full_seurat,group.by = "C29",label = T)
DimPlot(full_seurat,group.by = "C57",label = T)
DimPlot(full_seurat,group.by = "C92",label = T)
#p1+p2
#DimPlot(full_seurat)


group_select <- c("Mandibular_Maxillary", "Molar_Incisor", "Histology",  "Treatment",
                  "Age", "Stage", "Development.stage", "coarse_anno_1")

annoBindDf <- anno_bind(full_seurat,c("C2","C9","C19","C29","C57","C92"),group_select)

percentPlotFun(full_seurat, "C19", "Age")



df <- table(full_seurat$C6,full_seurat$Age.In.Detail.)%>%t%>%as.matrix()

df <- t(apply(df, 1, function(x) x/sum(x)))
df <- as.data.frame(df)

df_long <- df%>%rownames_to_column("Age")%>%
  pivot_longer(cols = -Age,names_to = "Cluster")


other_genes_remove = rownames(full_seurat@assays$originalexp@counts)[grepl("LOC|RGD|AA|AY|RP|Gm|Hb[ab]|Rik|-ps",rownames(full_seurat@assays$originalexp@counts))]
# make list of all genes that should be removed:
all_exclusion_genes = unique(c(other_genes_remove))



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
# command_args<-commandArgs(TRUE)
# param_file = "script/annotation/config/parameters_wtintegrate.json"
# # read all parameters and filepaths
# parameter_list = jsonlite::read_json(param_file)
# # if some fields are lists --> unlist
# parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
message("Starting annotation: ")
print(parameter_list$manual_names_annotation_abbr)

allgenes <- c(markers_comparisons_all$gene, markers_comparisons_siblings$gene) %>% unique()
excluded_genes1 <- grep("mt-" ,allgenes,value = T)
HbGene = grep("^Hb[a-b]", allgenes, value = T)
RbGene = grep("^Rp[ls]", allgenes, value = T)
RbMGene = grep("^Mrp[ls]", allgenes, value = T)
excluded_genes2 = grep("LOC", allgenes, value = T)
excluded_genes3 = grep("\\.", allgenes, value = T)
excluded_genes4 = allgenes[str_length(allgenes) > 13]
excluded_genes5 <- grep("Mt-" ,allgenes,value = T)
all_excluded_genes <- c(excluded_genes1,excluded_genes2,excluded_genes3,
                        excluded_genes4,excluded_genes5,HbGene,RbGene,RbMGene)

markers_comparisons_siblings <- markers_comparisons_siblings[!markers_comparisons_siblings$gene%in%all_excluded_genes,]
markers_comparisons_all <- markers_comparisons_all[!markers_comparisons_all$gene%in%all_excluded_genes,]



annotation_results = annotate_tree(edgelist = edgelist,#[1:32,],#edgelist = edgelist[1:291,],
                                   labelmat = labelmat_updated_new_labels,
                                   markers_comparisons_all = markers_comparisons_all,
                                   markers_comparisons_siblings = markers_comparisons_siblings,
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
labelmat_updated = cbind(Cell_ID=full_seurat@meta.data$Cell_ID,mrtree_result$labelmat)
# make a per cell cluster annotation
cell_cluster_map = labelmat_updated %>% as.data.frame() %>% tidyr::gather(-Cell_ID, key="clusterlevel",value="cluster_id")
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
annotation_df_wide = annotation_df_wide[match(full_seurat@meta.data$Cell_ID,annotation_df_wide$Cell_ID),]

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

data.table::fwrite(annotation_df,file = paste0(parameter_list$harmonization_folder_path,
                                               parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,
                                               "_annotation_result.txt") ,sep="\t")

data.table::fwrite(annotation_df_wide,file = paste0(parameter_list$harmonization_folder_path,
                                                    parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,
                                                    "_annotation_labelmat.txt") ,sep="\t")

data.table::fwrite(descriptive_markers_df,file = paste0(parameter_list$harmonization_folder_path,
                                                        parameter_list$new_name_suffix,"_",
                                                        parameter_list$marker_suffix,
                                                        "_annotation_markers_filtered.txt") ,sep="\t")


# #== annoated_manual----------------------------------
annotation_orig <- annotation_df_wide
rownames(annotation_orig) <- annotation_orig$Cell_ID
annotation_orig <- annotation_orig[-1]
annotation_orig <- annotation_orig%>%
  mutate_each(function(x) stringr::str_replace(x, "C[:digit:]+-[:digit:]+: ",""))

full_seurat <- AddMetaData(full_seurat,annotation_orig)
DimPlot(full_seurat,group.by = "C2_named")
write.csv(annotation_orig,paste0(process_dir,"20250217_annotation_meta.csv"),row.names = T,quote = F)
full_seurat@misc$annotation_result <- annotation_df
full_seurat@misc$clustering_edgelist <- edgelist

#== update score-------------------------
marker_update <- descriptive_markers_df%>%
  group_by(cluster_id)%>%
  top_n(1,wt = avg_score)%>%
  select(cluster_id,gene)
marker_update <- marker_update%>%column_to_rownames("cluster_id")
anno_bind <- function(object, cluster_select, group_select, threshold_add=0.2) {
  # Add error checking
  annolist <- lapply(cluster_select, function(x) {
    result <- tryCatch({
      anno_dict(object, cluster_id = x, group_select=group_select, threshold_add = threshold_add)
    }, error = function(e) {
      message(sprintf("Error processing cluster %s: %s", x, e$message))
      return(NULL)
    })

    # Check if result is valid before transposing
    if (!is.null(result) && length(result) > 0) {
      return(t(result))
    } else {
      return(NULL)
    }
  })

  # Remove NULL results
  annolist <- Filter(Negate(is.null), annolist)

  # Only proceed if we have valid results
  if (length(annolist) > 0) {
    annoBind <- do.call(rbind, annolist)
    return(as.data.frame(annoBind))
  } else {
    stop("No valid results to combine")
  }
}
annoBindDf <- anno_bind(full_seurat,c("C2","C9","C19","C29","C57","C92"),group_select)
annoBindDf <- merge(annoBindDf,marker_update, by = "row.names")
write.csv(annoBindDf,paste0(parameter_list$harmonization_folder_path,
                            parameter_list$new_name_suffix,"_",parameter_list$marker_suffix,
                            "_anno_bind_df.csv"))

saveRDS(full_seurat,"processed_data/integrated_data/20241106_mesenchyme.Rds")

