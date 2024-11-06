library(aplot)
source("script/utils/fun_plot_figure.R")
source("script/utils/fun_plotclustertree.R")
results_path_figure2="results/annotation/20241106_anno_step6_plot"
dir.create(results_path_figure2)
load_colors()
leaf_level_column = "C92"
leaf_level = 6
full_seurat <-

full_seurat_bk <- full_seurat
Idents(full_seurat) <- full_seurat$C19
#full_seurat <- full_seurat[,full_seurat$C19 != "C19-7"]
edgelist_bk <- full_seurat@misc$clustering_edgelist

edgelist <- edgelist%>%dplyr::filter(!(from %in% c("C19-7","C36-13","C49-16","C90-31","C90-32","C90-33") | to %in%  c("C19-7","C36-13","C49-16","C90-31","C90-32","C90-33")))
edgelist[26,1] <- "C7-5"
annotation_df <- annotation_df%>%dplyr::filter(!(cluster_id %in%  c("C19-7","C36-13","C49-16","C90-31","C90-32","C90-33") ))


#== change order of column-------------------

full_seurat$Age.In.Detail. <-factor(full_seurat$Age.In.Detail.)

newOrder <- full_seurat@meta.data %>% dplyr::select(Project,Age.In.Detail.)%>%
  table()%>%as.data.frame()%>%
  filter(Freq!=0)%>%
  arrange(Age.In.Detail.)
dput(as.character(newOrder$Project))

heatmap_data <- heatmap_data[,c("C137",newLevel)]


heatmap_matrix = as.matrix(heatmap_data[,2:ncol(heatmap_data)])
rownames(heatmap_matrix) = heatmap_data[,leaf_level_column]
heatmap_matrix[is.na(heatmap_matrix)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix)),Dataset = colnames(heatmap_matrix))
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames.txt"),sep="\t")
colnames(heatmap_matrix) = colnames_overview$number

# make data for second heatmap with n cells
heatmap_data2 = full_seurat@meta.data %>% dplyr::select(cellid,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column)) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "ncells")  %>% dplyr::ungroup()  %>% dplyr::mutate(ncells_pct = ncells / sum(ncells)*100)  %>% as.data.frame()
heatmap_matrix2 = as.matrix(heatmap_data2[,"ncells_pct",drop=FALSE])
colnames(heatmap_matrix2) = "%"
rownames(heatmap_matrix2) = heatmap_data2[,leaf_level_column]
heatmap_matrix2[is.na(heatmap_matrix2)] = 0

# make data for third heatmap with regions
# heatmap_data3 = full_seurat@meta.data %>% dplyr::select(cellid,Region_summarized,!!sym(leaf_level_column)) %>%
#   dplyr::group_by(!!sym(leaf_level_column),Region_summarized) %>% dplyr::count() %>% dplyr::group_by(!!sym(leaf_level_column)) %>%
#   dplyr::top_n(n = 1,wt = n) %>% ungroup() %>%
#   dplyr::distinct(!!sym(leaf_level_column),Region_summarized,.keep_all=TRUE) %>% as.data.frame()
# heatmap_matrix3 = as.matrix(heatmap_data3[,"Region_summarized",drop=F])
# rownames(heatmap_matrix3) = heatmap_data3[,leaf_level_column]
# colnames(heatmap_matrix3) = "R"


#== age matrix--------------
heatmap_data3 = full_seurat@meta.data %>% dplyr::select(cellid,Age.In.Detail.,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Age.In.Detail.) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  tidyr::spread(key = Age.In.Detail.,value=presence) %>% as.data.frame()



heatmap_matrix_age = as.matrix(heatmap_data3[,2:ncol(heatmap_data3)])
rownames(heatmap_matrix_age) = heatmap_data3[,leaf_level_column]
heatmap_matrix_age[is.na(heatmap_matrix_age)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix_age)),Dataset = colnames(heatmap_matrix_age))
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames_age.txt"),sep="\t")
colnames(heatmap_matrix_age) = colnames_overview$number



#== tissue matrix--------------
heatmap_data4 = full_seurat@meta.data %>% dplyr::select(cellid,Tissue.Specific.,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Tissue.Specific.) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  tidyr::spread(key = Tissue.Specific.,value=presence) %>% as.data.frame()


heatmap_matrix_tissue = as.matrix(heatmap_data4[,2:ncol(heatmap_data4)])
rownames(heatmap_matrix_tissue) = heatmap_data4[,leaf_level_column]
heatmap_matrix_tissue[is.na(heatmap_matrix_tissue)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix_tissue)),Dataset = colnames(heatmap_matrix_tissue))
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames_tissue.txt"),sep="\t")
colnames(heatmap_matrix_tissue) = colnames_overview$number

#== Organ matrix--------------
heatmap_data5 = full_seurat@meta.data %>% dplyr::select(cellid,Organ,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Organ) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  tidyr::spread(key = Organ,value=presence) %>% as.data.frame()


heatmap_matrix_organ = as.matrix(heatmap_data5[,2:ncol(heatmap_data5)])
rownames(heatmap_matrix_organ) = heatmap_data5[,leaf_level_column]
heatmap_matrix_organ[is.na(heatmap_matrix_organ)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix_organ)),Dataset = colnames(heatmap_matrix_organ))
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames_organ.txt"),sep="\t")
colnames(heatmap_matrix_organ) = colnames_overview$number


rgbList <- apply(heatmap_matrix_organ,1,function(x) rgb(x[[1]]/100,x[[2]]/100,x[[3]]/100))
heatmap_matrix3 = as.matrix(rgbList)
rownames(heatmap_matrix3) = rownames(heatmap_matrix_organ)
colnames(heatmap_matrix3) = "O"
rownames(heatmap_matrix3) = heatmap_data5[,leaf_level_column]

color_value_vector <- rgbList
names(color_value_vector) <- color_value_vector
circular_tree_heat = circular_tree_heat + ggplot2::scale_fill_manual(values = color_value_vector,na.value = "grey80")+
  guides(fill = FALSE)


## make annotation data frame
anno_df = annotation_df %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})

# plot cluster tree:

tree_color = "grey70"
circular_tree = plot_cluster_tree(edgelist = edgelist,
                                  leaf_level=8,
                                  anno_df = anno_df ,
                                  metadata=full_seurat@meta.data,
                                  label_size = 4,
                                  show_genes = TRUE,
                                  vjust_label = -0.5,
                                  edge_color = tree_color,
                                  node_color = tree_color)

circular_tree = rotate_tree(circular_tree, -90)
circular_tree
ggsave("result/4.26_annotation/5.4_annotree_circle.pdf",width = 8,height = 8)
# circular_tree = plot_cluster_tree(edgelist = edgelist,
#                                   leaf_level=6,
#                                   anno_df = anno_df ,
#                                   metadata=full_seurat@meta.data,
#                                   label_size = 2,
#                                   show_genes = TRUE,
#                                   vjust_label = -0.25,
#                                   edge_color = tree_color,
#                                   node_color = tree_color)

circular_tree = rotate_tree(circular_tree, -90)
# circular_tree
# circular_tree$data$tip.label <- circular_tree$data$cluster_name
circular_tree_heat = add_heatmap(circular_tree=circular_tree,
                                 heatmap_matrix = heatmap_matrix_age,
                                 heatmap_colors=c(bg_col,"darkred"),
                                 #scale_limits = c(0,100),
                                 heatmap_colnames =T,
                                 legend_title = "Pct Age",
                                 matrix_offset = 2.2,
                                 matrix_width =0.8,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=8,
                                 na_color = "white",
                                 heatmap_text_size=0)

circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix_tissue,
                                 heatmap_colors=c(bg_col,"darkorange"),
                                 #scale_limits = c(0,100),
                                 heatmap_colnames =T,
                                 legend_title = "Pct Tissue",
                                 matrix_offset = 7.8,
                                 matrix_width =0.35,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=20,
                                 na_color = "white",
                                 heatmap_text_size=0)
circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix_organ,
                                 heatmap_colors=c(bg_col,"dark green"),
                                 #scale_limits = c(0,100),
                                 heatmap_colnames =T,
                                 legend_title = "Pct Tissue",
                                 matrix_offset = 10.3,
                                 matrix_width =0.08,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=20,
                                 na_color = "white",
                                 heatmap_text_size=0)



heatmap_matrix3 = as.matrix(rgbList)
rownames(heatmap_matrix3) = rownames(heatmap_matrix_organ)
colnames(heatmap_matrix3) = "O"
circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix3,
                                 heatmap_colors=c(bg_col,"darkred"),
                                 heatmap_colnames =TRUE,
                                 legend_title = "Pct Dataset",
                                 matrix_offset = 10.8,
                                 matrix_width =0.05,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=12,
                                 na_color = "grey80",
                                 heatmap_text_size=3
)
circular_tree_heat + ggplot2::scale_fill_manual(values = color_value_vector,na.value = "grey80")+
  guides(fill = FALSE)


# plot tree with heatmaps
circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix2,
                                 heatmap_colors=c(bg_col,"gold"),
                                 scale_limits = c(0,2),
                                 heatmap_colnames =T,
                                 legend_title = "Cells",
                                 matrix_offset = 6.5,
                                 matrix_width =0.06,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=10,
                                 na_color = "white",
                                 heatmap_text_size=0)
circular_tree_heat


#== dotplot----------------------


#
#
# circular_tree = plot_cluster_tree(edgelist = edgelist,
#                                   leaf_level=6,
#                                   anno_df = anno_df ,
#                                   metadata=full_seurat@meta.data,
#                                   label_size = 4,
#                                   show_genes = TRUE,
#                                   vjust_label = -0.5,
#                                   edge_color = tree_color,
#                                   node_color = tree_color)

# p1 <- normaltree
#
# library(aplot)
# p3 <- p2%>%
#   insert_left(p1,width = 1)
# p3
# library(cowplot)
# plot_grid(p1, p2, nrow = 1, rel_widths = c(0.5,2), align = 'h')
#

normaltree <- ggtree(circular_tree$data)+
  geom_nodelab(aes(x=branch, label=first_cluster_name), color="darkred",vjust=-0.25,size=2.5)+
  geom_tiplab(ggplot2::aes(x=branch, label=first_cluster_name),color="darkred",vjust=-0.25,size=2.5)
normal_tree_heat = add_heatmap(circular_tree=normaltree,
                               heatmap_matrix = heatmap_matrix_age,
                               heatmap_colors=c(bg_col,"darkred"),
                               #scale_limits = c(0,100),
                               heatmap_colnames =T,
                               legend_title = "Pct Age",
                               matrix_offset = 0.5,
                               matrix_width =0.57,
                               colnames_angle=0,
                               legend_text_size = 8,
                               hjust_colnames=0,
                               na_color = "white",scale_limits = c(0,40),
                               heatmap_text_size=0)

normal_tree_heat = add_heatmap(circular_tree=normal_tree_heat,
                               heatmap_matrix = heatmap_matrix_tissue,
                               heatmap_colors=c(bg_col,"darkorange"),
                               #scale_limits = c(0,100),
                               heatmap_colnames =T,
                               legend_title = "Pct Tissue",
                               matrix_offset = 3.2,
                               matrix_width =0.25,
                               colnames_angle=0,
                               legend_text_size = 8,
                               hjust_colnames=-0.5,
                               na_color = "white",
                               heatmap_text_size=0)
normal_tree_heat = add_heatmap(circular_tree=normal_tree_heat,
                               heatmap_matrix = heatmap_matrix_organ,
                               heatmap_colors=c(bg_col,"dark green"),
                               #scale_limits = orangec(0,100),
                               heatmap_colnames =T,
                               legend_title = "Pct Tissue",
                               matrix_offset = 5.0,
                               matrix_width =0.03,
                               colnames_angle=0,
                               legend_text_size = 8,
                               hjust_colnames=-0.5,
                               na_color = "white",
                               heatmap_text_size=0)
normal_tree_heat = add_heatmap(circular_tree=normal_tree_heat,
                               heatmap_matrix = heatmap_matrix3,
                               heatmap_colors=c(bg_col,"darkred"),
                               heatmap_colnames =TRUE,
                               legend_title = "Pct Dataset",
                               matrix_offset =5.2,
                               matrix_width =0.02,
                               colnames_angle=0,
                               legend_text_size = 8,
                               hjust_colnames=12,
                               na_color = "grey80",
                               heatmap_text_size=3
)
normal_tree_heat=normal_tree_heat + ggplot2::scale_fill_manual(values = color_value_vector,na.value = "grey80")+
  guides(fill = FALSE)

Idents(full_seurat) <- full_seurat$C137
selectCluster <-levels(full_seurat)
index <- order(as.numeric(sub("C137-", "", selectCluster)))

# Use the index to get the sorted vector
selectCluster <- selectCluster[index]

normal_tree_heat
clusterOrder <- normaltree$data%>%
  dplyr::filter(clusterlevel=="C137")%>%
  arrange(y)%>%
  select(label)%>%
  unlist()



comparisonGene <- comparisons_all_update%>%
  dplyr::filter(cluster_id%in%selectCluster&!(gene%in%all_exclusion_genes))%>%
  group_by(cluster_id)%>%
  arrange(desc(specificity))%>%
  top_n(1)%>%
  arrange(fct_relevel(cluster_id, clusterOrder))%>%
  ungroup()%>%
  select(gene)%>%
  unlist()
levels(full_seurat) <- clusterOrder

p2 <- DotPlot(full_seurat,
              features = unique(comparisonGene),assay = "originalexp",
              cluster.idents = T,cols = c("lightgrey", "tan1"))+
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 70,size = 10,face = "bold",hjust  = 1)
  )
p3 <- p2%>%
  insert_left(normal_tree_heat,width = 0.8)
p3
ggsave("result/4.26_annotation/5.4_normaltree_full_dotplot.pdf",height = 15,width = 20)
saveRDS(p3,"result/4.26_annotation//5.4_normaltree_full_dotplot.Rds")
saveRDS(circular_tree_heat,"result/4.26_annotation/circular_hm.Rds")
saveRDS(full_seurat,"../important_processed_data/5.4_wtintegrate_full_seurat.Rds")



#== anno--------------------------

leaf_level = 5
leaf_level_column = "C49"

# full_seurat_bk <- full_seurat
# Idents(full_seurat) <- full_seurat$C19
# full_seurat <- full_seurat[,full_seurat$C19 != "C19-7"]
#
# #
#
# edgelist_bk <- full_seurat@misc$clustering_edgelist
#
# edgelist <- edgelist%>%dplyr::filter(!(from %in% c("C19-7","C36-13","C49-16","C90-31","C90-32","C90-33") | to %in%  c("C19-7","C36-13","C49-16","C90-31","C90-32","C90-33")))
# edgelist[26,1] <- "C7-5"
# annotation_df <- annotation_df%>%dplyr::filter(!(cluster_id %in%  c("C19-7","C36-13","C49-16","C90-31","C90-32","C90-33") ))


#== change order of column-------------------

full_seurat$Age.In.Detail. <-factor(full_seurat$Age.In.Detail.)

newOrder <- full_seurat@meta.data %>% dplyr::select(Project,Age.In.Detail.)%>%
  table()%>%as.data.frame()%>%
  filter(Freq!=0)%>%
  arrange(Age.In.Detail.)
dput(as.character(newOrder$Project))

heatmap_data <- heatmap_data[,c("C137",newLevel)]


heatmap_matrix = as.matrix(heatmap_data[,2:ncol(heatmap_data)])
rownames(heatmap_matrix) = heatmap_data[,leaf_level_column]
heatmap_matrix[is.na(heatmap_matrix)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix)),Dataset = colnames(heatmap_matrix))
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames.txt"),sep="\t")
colnames(heatmap_matrix) = colnames_overview$number

# make data for second heatmap with n cells
heatmap_data2 = full_seurat@meta.data %>% dplyr::select(cellid,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column)) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "ncells")  %>% dplyr::ungroup()  %>% dplyr::mutate(ncells_pct = ncells / sum(ncells)*100)  %>% as.data.frame()
heatmap_matrix2 = as.matrix(heatmap_data2[,"ncells_pct",drop=FALSE])
colnames(heatmap_matrix2) = "%"
rownames(heatmap_matrix2) = heatmap_data2[,leaf_level_column]
heatmap_matrix2[is.na(heatmap_matrix2)] = 0

# make data for third heatmap with regions
# heatmap_data3 = full_seurat@meta.data %>% dplyr::select(cellid,Region_summarized,!!sym(leaf_level_column)) %>%
#   dplyr::group_by(!!sym(leaf_level_column),Region_summarized) %>% dplyr::count() %>% dplyr::group_by(!!sym(leaf_level_column)) %>%
#   dplyr::top_n(n = 1,wt = n) %>% ungroup() %>%
#   dplyr::distinct(!!sym(leaf_level_column),Region_summarized,.keep_all=TRUE) %>% as.data.frame()
# heatmap_matrix3 = as.matrix(heatmap_data3[,"Region_summarized",drop=F])
# rownames(heatmap_matrix3) = heatmap_data3[,leaf_level_column]
# colnames(heatmap_matrix3) = "R"


#== age matrix--------------
heatmap_data3 = full_seurat@meta.data %>% dplyr::select(cellid,Age.In.Detail.,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Age.In.Detail.) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  tidyr::spread(key = Age.In.Detail.,value=presence) %>% as.data.frame()

heatmap_matrix_age = as.matrix(heatmap_data3[,2:ncol(heatmap_data3)])
rownames(heatmap_matrix_age) = heatmap_data3[,leaf_level_column]
heatmap_matrix_age[is.na(heatmap_matrix_age)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix_age)),Dataset = colnames(heatmap_matrix_age))
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames_age.txt"),sep="\t")
colnames(heatmap_matrix_age) = colnames_overview$number



#== tissue matrix--------------
heatmap_data4 = full_seurat@meta.data %>% dplyr::select(cellid,Tissue.Specific.,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Tissue.Specific.) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  tidyr::spread(key = Tissue.Specific.,value=presence) %>% as.data.frame()


heatmap_matrix_tissue = as.matrix(heatmap_data4[,2:ncol(heatmap_data4)])
rownames(heatmap_matrix_tissue) = heatmap_data4[,leaf_level_column]
heatmap_matrix_tissue[is.na(heatmap_matrix_tissue)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix_tissue)),Dataset = colnames(heatmap_matrix_tissue))
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames_tissue.txt"),sep="\t")
colnames(heatmap_matrix_tissue) = colnames_overview$number

#== Organ matrix--------------
heatmap_data5 = full_seurat@meta.data %>% dplyr::select(cellid,Organ,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Organ) %>% #dplyr::filter(predicted_Campbell!="NA")
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
  tidyr::spread(key = Organ,value=presence) %>% as.data.frame()


heatmap_matrix_organ = as.matrix(heatmap_data5[,2:ncol(heatmap_data5)])
rownames(heatmap_matrix_organ) = heatmap_data5[,leaf_level_column]
heatmap_matrix_organ[is.na(heatmap_matrix_organ)] = 0

colnames_overview = data.frame(number = 1:length(colnames(heatmap_matrix_organ)),Dataset = colnames(heatmap_matrix_organ))
data.table::fwrite(colnames_overview,paste0(results_path_figure2,"circular_tree_colnames_organ.txt"),sep="\t")
colnames(heatmap_matrix_organ) = colnames_overview$number


rgbList <- apply(heatmap_matrix_organ,1,function(x) rgb(x[[1]]/100,x[[2]]/100,x[[3]]/100))
heatmap_matrix3 = as.matrix(rgbList)
rownames(heatmap_matrix3) = rownames(heatmap_matrix_organ)
colnames(heatmap_matrix3) = "O"
rownames(heatmap_matrix3) = heatmap_data5[,leaf_level_column]

color_value_vector <- rgbList
names(color_value_vector) <- color_value_vector
circular_tree_heat = circular_tree_heat + ggplot2::scale_fill_manual(values = color_value_vector,na.value = "grey80")+
  guides(fill = FALSE)


## make annotation data frame
anno_df = annotation_df %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})

# plot cluster tree:

tree_color = "grey70"
circular_tree = plot_cluster_tree(edgelist = edgelist,
                                  leaf_level=6,
                                  anno_df = anno_df ,
                                  metadata=full_seurat@meta.data,
                                  label_size = 3,
                                  show_genes = TRUE,
                                  vjust_label = -0.5,
                                  edge_color = tree_color,
                                  node_color = tree_color)

circular_tree = rotate_tree(circular_tree, -90)
circular_tree
ggsave("result/4.26_annotation/5.4_annotree_circle_level5.pdf",width = 8,height = 8)
# circular_tree = plot_cluster_tree(edgelist = edgelist,
#                                   leaf_level=6,
#                                   anno_df = anno_df ,
#                                   metadata=full_seurat@meta.data,
#                                   label_size = 2,
#                                   show_genes = TRUE,
#                                   vjust_label = -0.25,
#                                   edge_color = tree_color,
#                                   node_color = tree_color)

# circular_tree
# circular_tree$data$tip.label <- circular_tree$data$cluster_name
circular_tree_heat = add_heatmap(circular_tree=circular_tree,
                                 heatmap_matrix = heatmap_matrix_age,
                                 heatmap_colors=c(bg_col,"darkred"),
                                 #scale_limits = c(0,100),
                                 heatmap_colnames =T,
                                 legend_title = "Pct Age",
                                 matrix_offset = 3.5,
                                 matrix_width =1.2,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=8,
                                 na_color = "white",
                                 heatmap_text_size=0)

circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix_tissue,
                                 heatmap_colors=c(bg_col,"darkorange"),
                                 #scale_limits = c(0,100),
                                 heatmap_colnames =T,
                                 legend_title = "Pct Tissue",
                                 matrix_offset = 9.5,
                                 matrix_width =0.7,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=20,
                                 na_color = "white",
                                 heatmap_text_size=0)
circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix_organ,
                                 heatmap_colors=c(bg_col,"dark green"),
                                 #scale_limits = c(0,100),
                                 heatmap_colnames =T,
                                 legend_title = "Pct Tissue",
                                 matrix_offset = 13.0,
                                 matrix_width =0.08,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=20,
                                 na_color = "white",
                                 heatmap_text_size=0)



heatmap_matrix3 = as.matrix(rgbList)
rownames(heatmap_matrix3) = rownames(heatmap_matrix_organ)
colnames(heatmap_matrix3) = " "
circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix3,
                                 heatmap_colors=c(bg_col,"darkred"),
                                 heatmap_colnames =TRUE,
                                 legend_title = "Pct Dataset",
                                 matrix_offset = 13.4,
                                 matrix_width =0.05,
                                 colnames_angle=0,
                                 legend_text_size = 7.2,
                                 hjust_colnames=12,
                                 na_color = "grey80",
                                 heatmap_text_size=3
)
circular_tree_heat + ggplot2::scale_fill_manual(values = color_value_vector,na.value = "grey80")+
  guides(fill = FALSE)
ggsave("result/4.26_annotation/5.4_annotree_circle_hm_level5.pdf",width = 15,height = 15)

# plot tree with heatmaps
circular_tree_heat = add_heatmap(circular_tree=circular_tree_heat,
                                 heatmap_matrix = heatmap_matrix2,
                                 heatmap_colors=c(bg_col,"gold"),
                                 scale_limits = c(0,2),
                                 heatmap_colnames =T,
                                 legend_title = "Cells",
                                 matrix_offset = 6.5,
                                 matrix_width =0.06,
                                 colnames_angle=0,
                                 legend_text_size = 8,
                                 hjust_colnames=10,
                                 na_color = "white",
                                 heatmap_text_size=0)
circular_tree_heat


#== dotplot----------------------


#
#
# circular_tree = plot_cluster_tree(edgelist = edgelist,
#                                   leaf_level=6,
#                                   anno_df = anno_df ,
#                                   metadata=full_seurat@meta.data,
#                                   label_size = 4,
#                                   show_genes = TRUE,
#                                   vjust_label = -0.5,
#                                   edge_color = tree_color,
#                                   node_color = tree_color)

# p1 <- normaltree
#
# library(aplot)
# p3 <- p2%>%
#   insert_left(p1,width = 1)
# p3
# library(cowplot)
# plot_grid(p1, p2, nrow = 1, rel_widths = c(0.5,2), align = 'h')
#

normaltree <- ggtree(circular_tree$data)+
  geom_nodelab(aes(x=branch, label=first_cluster_name), color="darkred",vjust=-0.25,size=3)+
  geom_tiplab(ggplot2::aes(x=branch, label=first_cluster_name),color="darkred",vjust=-0.25,size=3)
normal_tree_heat = add_heatmap(circular_tree=normaltree,
                               heatmap_matrix = heatmap_matrix_age,
                               heatmap_colors=c(bg_col,"darkred"),
                               #scale_limits = c(0,100),
                               heatmap_colnames =T,
                               legend_title = "Pct Age",
                               matrix_offset = 0.5,
                               matrix_width =0.6,
                               colnames_angle=0,
                               legend_text_size = 8,
                               hjust_colnames=0,
                               na_color = "white",scale_limits = c(0,40),
                               heatmap_text_size=0)

normal_tree_heat = add_heatmap(circular_tree=normal_tree_heat,
                               heatmap_matrix = heatmap_matrix_tissue,
                               heatmap_colors=c(bg_col,"darkorange"),
                               #scale_limits = c(0,100),
                               heatmap_colnames =T,
                               legend_title = "Pct Tissue",
                               matrix_offset = 3,
                               matrix_width =0.4,
                               colnames_angle=0,
                               legend_text_size = 8,
                               hjust_colnames=-1,
                               na_color = "white",
                               heatmap_text_size=0)
normal_tree_heat = add_heatmap(circular_tree=normal_tree_heat,
                               heatmap_matrix = heatmap_matrix_organ,
                               heatmap_colors=c(bg_col,"dark green"),
                               #scale_limits = orangec(0,100),
                               heatmap_colnames =T,
                               legend_title = "Pct Tissue",
                               matrix_offset = 5.0,
                               matrix_width =0.05,
                               colnames_angle=0,
                               legend_text_size = 8,
                               hjust_colnames=-0.5,
                               na_color = "white",
                               heatmap_text_size=0)
normal_tree_heat = add_heatmap(circular_tree=normal_tree_heat,
                               heatmap_matrix = heatmap_matrix3,
                               heatmap_colors=c(bg_col,"darkred"),
                               heatmap_colnames =TRUE,
                               legend_title = "Pct Dataset",
                               matrix_offset =5.3,
                               matrix_width =0.02,
                               colnames_angle=0,
                               legend_text_size = 8,
                               hjust_colnames=12,
                               na_color = "grey80",
                               heatmap_text_size=3
)
normal_tree_heat=normal_tree_heat + ggplot2::scale_fill_manual(values = color_value_vector,na.value = "grey80")+
  guides(fill = FALSE)
normal_tree_heat
ggsave("result/4.26_annotation/5.4_normaltree_hm.pdf",width = 12,height = 8)
Idents(full_seurat) <- full_seurat$C49
selectCluster <-levels(full_seurat)
index <- order(as.numeric(sub("C49-", "", selectCluster)))

# Use the index to get the sorted vector
selectCluster <- selectCluster[index]


clusterOrder <- normaltree$data%>%
  dplyr::filter(clusterlevel=="C49")%>%
  arrange(y)%>%
  select(label)%>%
  unlist()



comparisonGene <- comparisons_all_update%>%
  dplyr::filter(cluster_id%in%selectCluster&!(gene%in%all_exclusion_genes))%>%
  group_by(cluster_id)%>%
  arrange(desc(specificity))%>%
  top_n(1)%>%
  arrange(fct_relevel(cluster_id, clusterOrder))%>%
  ungroup()%>%
  select(gene)%>%
  unlist()
levels(full_seurat) <- clusterOrder

p2 <- DotPlot(full_seurat,
              features = unique(comparisonGene),assay = "originalexp",
              cluster.idents = T,cols = c("lightgrey", "tan1"))+
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 70,size = 10,face = "bold",hjust  = 1)
  )
p3 <- p2%>%
  insert_left(normal_tree_heat,width = 1.5)
p3
ggsave("result/5.4_normaltree_dotplot.pdf",width = 15,height = 8)
saveRDS(edgelist,"result/4.26_annotation/5.4update_edgelist.Rds")
