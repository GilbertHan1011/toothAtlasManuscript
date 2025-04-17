fitHM <- read.csv("process/trajectory/20250414_epi_run_4/fitted_trajectories.csv",row.names = 1)
fitHM <- t(fitHM)
Heatmap(fitHM[c("Ambn","Amtn","Enam","Mcm3","Gpr45","Rpl7"),],cluster_rows = F,cluster_columns = F)
rownames(fitHM)%>%head
conserved <- read.csv("process/trajectory/20250414_epi_run_4/conservation_scores.csv",row.names = 1)

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
write.csv(resAll,"process/trajectory/20250414_epi_run_4/hm_cluster_4000.csv")

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


res_50 <- compare_go_enrichment(top_n = 50,dtw_genes = conservation_score_filter$gene)
dotplot(res_50,showCategory=20)

res_100 <- compare_go_enrichment(top_n = 100,dtw_genes = conservation_score_filter$gene)
dotplot(res_100,showCategory=20)


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

res_tradeseq2 <- compare_go_enrichment(top_n = 69,de_genes = tradeseq_diff)
dotplot(res_tradeseq2,showCategory = 20)
conservation_score$rank <- 1:length(conservation_score$gene)
conservation_score_filter_raw <- conservation_score%>%
  filter(gene %in%c(cluster1,cluster2))
res_tradeseq4 <- compare_go_enrichment(top_n = length(DiffGene),de_genes = DiffGene,dtw_genes = conservation_score_filter_raw$gene)
dotplot(res_tradeseq4,showCategory = 10)


hist(conservation_score$raw_score)


create_elbow_plot <- function(scores, title = "Ranked Elbow Plot",
                              xlab = "Rank", ylab = "Score",
                              highlight_elbow = TRUE, top_n = NULL,
                              point_size = 1, line_size = 1) {

  # Sort values (use negative for ascending order if they're already negative)
  sorted_scores <- sort(scores)

  # Create ranks
  ranks <- seq_along(sorted_scores)

  # Create data frame for plotting
  plot_data <- data.frame(rank = ranks, score = sorted_scores)

  # Create basic plot
  library(ggplot2)
  p <- ggplot(plot_data, aes(x = rank, y = score)) +
    geom_line(size = line_size) +
    geom_point(size = point_size, alpha = 0.5) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal()

  # Try to detect elbow point if requested
  if (highlight_elbow) {
    # Simple method: find point of maximum curvature
    # Convert to cumulative sums first to make the elbow more pronounced
    y_norm <- cumsum(sorted_scores)
    y_norm <- y_norm / max(y_norm)
    x_norm <- ranks / max(ranks)

    # Calculate curvature (approximate)
    # This is a simple approach and may need refinement
    dx <- c(0, diff(x_norm))
    dy <- c(0, diff(y_norm))
    dx[dx == 0] <- .Machine$double.eps  # Avoid division by zero

    # Second derivative approximation
    d2y <- c(0, diff(dy/dx))

    # Find the point of maximum curvature
    elbow_idx <- which.max(abs(d2y))

    # Add elbow point to the plot
    p <- p +
      geom_vline(xintercept = elbow_idx, linetype = "dashed", color = "red") +
      annotate("text", x = elbow_idx, y = min(sorted_scores) + 0.1 * diff(range(sorted_scores)),
               label = paste("Elbow at rank", elbow_idx), hjust = -0.1, color = "red")
  }

  # Highlight top N if specified
  if (!is.null(top_n) && top_n < length(scores)) {
    threshold_value <- sorted_scores[top_n]
    p <- p +
      geom_hline(yintercept = threshold_value, linetype = "dashed", color = "blue") +
      annotate("text", x = max(ranks) * 0.8, y = threshold_value,
               label = paste("Top", top_n, "threshold"), vjust = -0.5, color = "blue")
  }

  return(p)
}

# Example usage with your data
# First, get all conservation scores (not just the first 6)
conservation_scores <- conservation_score$raw_score

# Create the elbow plot
elbow_plot <- create_elbow_plot(
  conservation_scores,
  title = "Conservation Score Elbow Plot",
  ylab = "Conservation Score (raw)",
  highlight_elbow = TRUE,
  top_n = 500  # Highlight top 500 genes
)
