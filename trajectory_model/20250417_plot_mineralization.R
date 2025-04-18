mine_conservation = read.csv("process/trajectory/20250417_mineralization/conservation_scores.csv",row.names = 1)
mine_conservation <- mine_conservation %>% dplyr::filter(n_valid_samples>=15)%>%
  head(110)

mineFitHM <- read.csv("process/trajectory/20250417_mineralization/fitted_trajectories_optimized.csv",row.names = 1)
mineFitHM <- mineFitHM[mine_conservation$gene]%>%t()
Heatmap(mineFitHM,cluster_rows = F,cluster_columns = F,km = 4,show_row_names = F,show_column_names = F)
# Perform k-means clustering
k <- 8  # Number of clusters
km <- kmeans(mineFitHM, centers = k, nstart = 25)


# Create a vector of cluster assignments for all rows
cluster_assignments <- km$cluster
# Create the heatmap with clustering results and column split
hmAll <- Heatmap(mineFitHM,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 row_split = cluster_assignments,
                 col = colorRamp2(c(-1, 0, 1), c("Deepskyblue3", "white", "red")),
                 na_col = "grey",
                 border = TRUE)
hmAll <- draw(hmAll)

cluster_assignments_df <- as.data.frame(cluster_assignments)
cluster_assignments_df$gene <- rownames(cluster_assignments_df)
write.csv(cluster_assignments_df,"process/trajectory/20250417_mineralization/20250417_cluster_mineralization.csv")
genedf = as.data.frame(cluster_assignments_df$gene)
colnames(genedf)= "Gene"
write.csv(genedf,"process/trajectory/20250417_mineralization/gene.csv",quote = F)

mineFitHM_ordered <- mineFitHM[order(apply(mineFitHM, 1, which.max)), ]
Heatmap(mineFitHM_ordered,
                 cluster_columns = FALSE,
                 cluster_rows = F,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 #row_split = cluster_assignments,
                 col = colorRamp2(c(-1, 0, 1), c("Deepskyblue3", "white", "red")),
                 na_col = "grey",
                 border = TRUE)
