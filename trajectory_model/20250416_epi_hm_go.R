fitHM <- read.csv("process/trajectory/20250414_epi_run_3/fitted_trajectories.csv",row.names = 1)
fitHM <- t(fitHM)
Heatmap(fitHM[c("Ambn","Amtn","Enam","Mcm3","Gpr45","Rpl7"),],cluster_rows = F,cluster_columns = F)
rownames(fitHM)%>%head
conserved <- read.csv("process/trajectory/20250414_epi_run_3/conservation_scores.csv",row.names = 1)

Heatmap(fitHM,km = 8,cluster_rows = F,cluster_columns = F)


# Assuming fitHM is your matrix/data frame for the heatmap
# Create the heatmap with k-means clustering (km=8)
set.seed(123)  # For reproducibility of clustering
ht <- Heatmap(fitHM,
              name = "Values",
              km = 8,  # k-means clustering with 8 clusters
              cluster_rows = F,
              cluster_columns = F,
              show_row_names = F,
              show_column_names = F,
              col=colorRamp2(c(-1, 0, 1), c("DeepSkyBlue3", "white", "red")))

# Draw the heatmap
ht_drawn <- draw(ht)

hmOrder <- row_order(ht_drawn)
hmGeneAll <- rownames(ht_drawn@ht_list$Values@matrix)
hmGeneListAll <- lapply(hmOrder, function(x) hmGeneAll[x])
library(stringi)
resAll <- as.data.frame((stri_list2matrix(hmGeneListAll)))
colnames(resAll) <- names(hmGeneListAll)
write.csv(resAll,"process/trajectory/20250414_epi_run_3/hm_cluster_4000.csv")

cluster1 <- hmGeneListAll[[1]]
cluster2 <- hmGeneListAll[[2]]

varGene <- read.csv("processed_data/framework/geneMeta/20250415_epi_vargene.csv",row.names = 1)%>%unlist
cluster1_var <- intersect(varGene,cluster1)
cluster2_var <- intersect(varGene,cluster2)
cluster_union <- c(cluster1_var,cluster2_var)


de_gene <- read.csv("process/trajectory/20250415_GO_downstream/20250415_epimarker_seurat_Filter.csv",row.names = 1)
de_gene_list <- de_gene%>%
  arrange(desc(avg_log2FC))%>%
  rownames()
de_gene_select = head(de_gene_list,length(cluster1_var))
compareList = list("DE" = de_gene_select,"DTW" = cluster1_var)
res <- clusterProfiler::compareCluster(compareList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(res,showCategory=20)

de_gene_select2 = head(de_gene_list,length(cluster_union))
compareUnion = list("DE" = de_gene_select2,"DTW" = cluster_union)
res2 <- clusterProfiler::compareCluster(compareUnion, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(res2,showCategory=30)

conservation_score = read.csv("process/trajectory/20250414_epi_run_3/conservation_scores.csv",row.names = 1)
conservation_score_filter <- conservation_score%>%
  filter(gene%in%varGene & gene %in%cluster_union)
DTW_gene_100 <- conservation_score_filter$gene%>% head(100)
de_gene_100 <- de_gene_list%>%head(100)
compareUnion3 = list("DE" = de_gene_100,"DTW" = DTW_gene_100)
res3 <- clusterProfiler::compareCluster(compareUnion3, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(res3,showCategory=30)


DTW_gene_50 <- conservation_score_filter$gene%>% head(50)
de_gene_50 <- de_gene_list%>%head(50)
compareUnion4 = list("DE" = de_gene_50,"DTW" = DTW_gene_50)
res4 <- clusterProfiler::compareCluster(compareUnion4, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(res4,showCategory=20)
compare_go_enrichment <- function(
                                  dtw_genes = conservation_score_filter$gene,
                                  de_genes = de_gene_list,
                                  top_n = 50,
                                  org_db = "org.Mm.eg.db",
                                  ontology = "ALL",
                                  key_type = "SYMBOL",
                                  ...) {

  # Check if required packages are installed
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is needed for this function to work. Please install it.")
  }

  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop(paste0("Package '", org_db, "' is needed for this function to work. Please install it."))
  }

  # Take top N genes from each list
  dtw_genes_top <- head(dtw_genes, top_n)
  de_genes_top <- head(de_genes, top_n)

  # Create named list for compareCluster
  gene_lists <- list("TradeSeq" = de_genes_top, "DTW" = dtw_genes_top)

  # Run compareCluster with enrichGO
  result <- clusterProfiler::compareCluster(
    gene_lists,
    fun = clusterProfiler::enrichGO,
    OrgDb = org_db,
    ont = ontology,
    keyType = key_type,
    ...
  )

  return(result)
}
res_150 <- compare_go_enrichment(top_n = 150)
dotplot(res_150,showCategory = 20)

res_20 <- compare_go_enrichment(top_n = 20)
dotplot(res_20,showCategory = 20)


res_200 <- compare_go_enrichment(top_n = 200)
dotplot(res_200,showCategory = 20)



gene_cluster <- read.csv("process/lzs_results/20250404_gene_cluster.csv",row.names = 1)
DiffGene <- gene_cluster$gene[gene_cluster$sigGene_ordered%in%c("GC8","GC7")]
tradeseq_diff <- intersect(DiffGene,varGene)
res_tradeseq <- compare_go_enrichment(top_n = 69,de_genes = tradeseq_diff)
dotplot(res_tradeseq,showCategory = 20)

