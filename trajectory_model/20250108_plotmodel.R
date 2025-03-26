combined <- readRDS("processed_data/trajectory/20250106_all_combined.Rds")

genelist <- list.files("results/trajectory/20250107_bone_tooth_conserved_plot/20250108_combined/",full.names = T,pattern = "*weight.csv")
geneName <- list.files("results/trajectory/20250107_bone_conserved/",pattern = "*weight.csv") %>% gsub("_fit_weight.csv","",.) %>% gsub("gene_","",.)
files <- lapply(genelist,read.csv, row.names = 1)
weightEstimate <- lapply(files, function(x) x$Estimate) %>% as.data.frame() %>% t()
rownames(weightEstimate) <- geneName
ComplexHeatmap::Heatmap(weightEstimate,show_row_names = F)
colnames(weightEstimate) <- files[[1]] %>% rownames() %>% gsub("b_shape_array","",.)
prefix = lapply(strsplit(colnames(weightEstimate),"_"), `[`, 1) %>% unlist()


n_clusters = 8
km <- kmeans(weightEstimate, centers = n_clusters)

# 2. Create heatmap with k-means clustering
ht <- ComplexHeatmap::Heatmap(
  weightEstimate,
  show_row_names = FALSE,
  row_split = km$cluster,  # split by k-means clusters
  cluster_rows = T,     # disable clustering
)
pdf("results/trajectory/20250107_bone_tooth_conserved_plot/20250107_heatmap_2nd.pdf",width = 6,height = 10)
draw(ht)
dev.off()
row_order <- row_order(ht)
row_order_vector <- unlist(row_order)
# Create ordered data frame
ordered_data <- data.frame(
  row_names = rownames(weightEstimate)[row_order_vector],
  cluster = km$cluster[row_order_vector],
  row_order = row_order_vector
)

# If you want to keep track of which list element (split) each row came from
ordered_data <- data.frame(
  row_names = rownames(weightEstimate)[row_order_vector],
  cluster = km$cluster[row_order_vector],
  row_order = row_order_vector,
  split = rep(names(row_order), sapply(row_order, length))
)
write.csv(ordered_data,"results/trajectory/20250107_bone_tooth_conserved_plot/20250107_vargeneHM_order2.csv")




pdf("results/trajectory/20250107_bone_tooth_conserved_plot/20250107_heatmap_2nd.pdf",width = 6,height = 10)
ht <- ComplexHeatmap::Heatmap(
  weightEstimate,
  show_row_names = FALSE,
  row_split = km$cluster,  # split by k-means clusters
  cluster_rows = T,     # disable clustering
  column_split = prefix,show_column_names = F
)
draw(ht)
dev.off()

cluster7Gene <- names(km$cluster)[km$cluster==7]


predictionFile <- list.files("results/trajectory/20250107_bone_tooth_conserved_plot/20250108_combined/",full.names = T,pattern = "*predictions.csv")
geneName2 <- list.files("results/trajectory/20250107_bone_tooth_conserved_plot/20250108_combined/",pattern = "*predictions.csv") %>% gsub("_predictions.csv","",.) %>% gsub("gene_","",.)
predfiles <- lapply(predictionFile,read.csv, row.names = 1)

predFit <- lapply(predfiles, function(x) x$fit) %>% as.data.frame() %>% t
rownames(predFit) <- geneName2
colnames(predFit) <- 1:100
predFitCluster7 <- predFit[cluster7Gene,]
predFitCluster7_scale <- t(scale(t(predFitCluster7)))
ComplexHeatmap::Heatmap(predFitCluster7_scale,cluster_columns = F,km = 8)


# Create and draw the heatmap
ht = Heatmap(predFitCluster7_scale,
             cluster_columns = FALSE,
             km = 8)
ht_drawn = draw(ht)

# Get row order numbers
row_order_nums = row_order(ht_drawn)

# Convert numbers to gene names for each cluster
gene_clusters = lapply(row_order_nums, function(x) {
  rownames(predFitCluster7_scale)[x]
})

# Create data frame with gene names and cluster assignments
cluster_df = data.frame(
  gene = unlist(gene_clusters),
  cluster = rep(1:length(gene_clusters), sapply(gene_clusters, length)),
  row.names = NULL
)
write.csv(cluster_df,"../202409_tooth/results/trajectory/20250107_bone_tooth_conserved_plot/20250108_cluster7_gene.csv")


Heatmap(predFitCluster7,
        cluster_columns = FALSE,
        km = 8)


row_maxes = apply(predFitCluster7, 1, max)

subsetGene = names(row_maxes)[row_maxes>1]

# Create and draw the heatmap
ht = Heatmap(predFitCluster7_scale[subsetGene,],
             cluster_columns = FALSE,
             km = 8)
ht_drawn = draw(ht)

# Get row order numbers
row_order_nums = row_order(ht_drawn)

# Convert numbers to gene names for each cluster
gene_clusters = lapply(row_order_nums, function(x) {
  rownames(predFitCluster7_scale)[x]
})

# Create data frame with gene names and cluster assignments
cluster_df = data.frame(
  gene = unlist(gene_clusters),
  cluster = rep(1:length(gene_clusters), sapply(gene_clusters, length)),
  row.names = NULL
)

FeaturePlot(mes,"Vim")


#== subset varGene---------
varGeneMes <- read.csv("processed_data/framework/hvg/20250108_hvg.csv",row.names = 1) %>% unlist
varGeneEpi <- read.csv("processed_data/framework/hvg/20250108_hvgEpi.csv",row.names = 1) %>% unlist
varGeneBone <- read.csv("processed_data/framework/hvg/20250108_hvgBone.csv",row.names = 1) %>% unlist
hvgGene <- reduce(list(varGeneMes,varGeneEpi,varGeneEpi),intersect)
hvgSub <- intersect(hvgGene,subsetGene)
write.csv(hvgSub,"processed_data/trajectory/20250108_cluster7_gene.csv")

#== load house keeping gene--------------------
# https://github.com/brianpenghe/Matlab-genomics/blob/master/He_2020_ENCODE3_RNA/GeneLists/Bulk%20Cluster%20Ubiquitous.txt
# https://github.com/Bidossessih/HRT_Atlas/blob/master/www/Housekeeping_GenesMouse.csv
housekeep1 <- read.csv("data/genes/Housekeeping_GenesMouse_.csv",sep = ";")
housekeep1 <- housekeep1$Gene
housekeep2 <- data.table::fread("data/genes/Bulk Cluster Ubiquitous (1).txt") %>% unlist
housekeep <- c(housekeep1,housekeep2) %>% unique()

#== mtGene---------
mtGene = grep("^mt-*", geneName2,value = T)
HbGene = grep("^Hb[a-b]", geneName2, value = T)
RbGene = grep("^Rp[ls]", geneName2, value = T)
RbMGene = grep("^Mrp[ls]", geneName2, value = T)

allRemoveGene <- c(housekeep,mtGene,HbGene,RbGene,RbMGene) %>% unique

# low expression genes------------------
highExpGene  = apply(predFit, 1, max)
highExpGene <- names(highExpGene)[highExpGene>0.3]


#== new genelist--------------------

keepGene <- setdiff(geneName2,allRemoveGene)
keepGene <- intersect(highExpGene,keepGene)
#weightEstimateSub <- weightEstimate[rownames(weightEstimate) %in% keepGene,]



#== cluster7-------------------
cluster7Gene <- names(km$cluster)[km$cluster==7]
cluster7GeneSub <- intersect(cluster7Gene,keepGene)
predSubCluster7 <- predFit[cluster7GeneSub,]
predSubCluster7_scale <- t(scale(t(predSubCluster7)))
ht7 = Heatmap(predSubCluster7_scale,
             cluster_columns = FALSE,
             km = 8)

ht_drawn7 = draw(ht7)

# Get row order numbers
row_order_nums = row_order(ht_drawn7)

# Convert numbers to gene names for each cluster
gene_clusters = lapply(row_order_nums, function(x) {
  rownames(predSubCluster7_scale)[x]
})

# Create data frame with gene names and cluster assignments
cluster_df = data.frame(
  gene = unlist(gene_clusters),
  cluster = rep(1:length(gene_clusters), sapply(gene_clusters, length)),
  row.names = NULL
)

get_heatmap_clusters <- function(expression_matrix,
                                 genes_to_use = NULL,
                                 n_clusters = 8,
                                 scale_data = TRUE,
                                 save_path = NULL,  # New parameter for save path
                                 fig_width = 6,     # New parameter for figure width
                                 fig_height = 10) { # New parameter for figure height

  # Subset matrix if genes provided
  if (!is.null(genes_to_use)) {
    expression_matrix <- expression_matrix[intersect(rownames(expression_matrix), genes_to_use), ]
  }

  # Scale data if requested
  if (scale_data) {
    expression_matrix <- t(scale(t(expression_matrix)))
  }

  # Create heatmap
  ht = Heatmap(expression_matrix,
               cluster_columns = FALSE,
               km = n_clusters)

  # Save and draw heatmap if save_path is provided
  if (!is.null(save_path)) {
    pdf(save_path, width = fig_width, height = fig_height)
    ht_drawn = draw(ht)
    dev.off()

    # Draw again to keep in memory
    ht_drawn
  } else {
    # Just draw if not saving
    ht_drawn = draw(ht)
  }

  # Get row order numbers
  row_order_nums = row_order(ht_drawn)

  # Convert numbers to gene names for each cluster
  gene_clusters = lapply(row_order_nums, function(x) {
    rownames(expression_matrix)[x]
  })

  # Create data frame with gene names and cluster assignments
  cluster_df = data.frame(
    gene = unlist(gene_clusters),
    cluster = rep(1:length(gene_clusters), sapply(gene_clusters, length)),
    row.names = NULL
  )

  # Return list with results
  return(list(
    heatmap = ht,
    drawn_heatmap = ht_drawn,
    cluster_assignments = cluster_df,
    scaled_matrix = expression_matrix
  ))
}

results7 <- get_heatmap_clusters(
  expression_matrix = predFit[cluster7Gene,],
  genes_to_use = keepGene,
  n_clusters = 8,
  save_path = "results/trajectory/20250107_bone_tooth_conserved_plot/hm_conserved_gene_7.pdf",
  fig_width = 6,
  fig_height = 10
)
cluster7 <- results7$cluster_assignments
candidateGene <- cluster7[cluster7$cluster%in%c(1,4,2,3),]
write.csv(candidateGene,"results/trajectory/20250107_bone_tooth_conserved_plot/candidate_gene/cluster7.csv")


cluster6Gene <- names(km$cluster)[km$cluster==6]
results6 <- get_heatmap_clusters(
  expression_matrix = predFit[cluster6Gene,],
  genes_to_use = keepGene,
  n_clusters = 8,
  save_path = "results/trajectory/20250107_bone_tooth_conserved_plot/hm_conserved_gene_6.pdf",
  fig_width = 6,
  fig_height = 10
)

cluster6 <- results6$cluster_assignments
candidateGene6 <- cluster6[cluster6$cluster%in%c(1,4,2,3,5),]
write.csv(candidateGene6,"results/trajectory/20250107_bone_tooth_conserved_plot/candidate_gene/cluster6.csv")


cluster4Gene <- names(km$cluster)[km$cluster==4]
results4 <- get_heatmap_clusters(
  expression_matrix = predFit[cluster4Gene,],
  genes_to_use = keepGene,
  n_clusters = 8,
  save_path = "results/trajectory/20250107_bone_tooth_conserved_plot/hm_conserved_gene_4.pdf",
  fig_width = 6,
  fig_height = 10
)

cluster4 <- results4$cluster_assignments
candidateGene4 <- cluster4[cluster4$cluster%in%c(1,4,2,3,5),]
write.csv(candidateGene4,"results/trajectory/20250107_bone_tooth_conserved_plot/candidate_gene/cluster4.csv")


cluster2Gene <- names(km$cluster)[km$cluster==2]
results2 <- get_heatmap_clusters(
  expression_matrix = predFit[cluster2Gene,],
  genes_to_use = keepGene,
  n_clusters = 8,
  save_path = "results/trajectory/20250107_bone_tooth_conserved_plot/hm_conserved_gene_2.pdf",
  fig_width = 6,
  fig_height = 10
)

cluster2 <- results2$cluster_assignments
candidateGene2 <- cluster2[cluster2$cluster%in%c(1,4,2,3,5),]
write.csv(candidateGene2,"results/trajectory/20250107_bone_tooth_conserved_plot/candidate_gene/cluster2.csv")

cluster7$prior <- 1
cluster6$prior <- 2
cluster4$prior <- 3
cluster2$prior <- 4

candidateList <- do.call(rbind,list(cluster7,cluster6,cluster4,cluster2))
write.csv(candidateList[c(1,3)],"results/trajectory/20250107_bone_tooth_conserved_plot/candidate_gene/2025_gene_combine.csv")

lapply(list(cluster7,cluster6,cluster4,cluster2), dim)
195
