fitHM_odonto <- read.csv("process/trajectory/20250415_odonto_run_2/fitted_trajectories_optimized.csv",row.names = 1)
fitHM_odonto <- t(fitHM_odonto)
Heatmap(fitHM_odonto[c("Sp7","Bglap3","Aspn","Sp6"),],cluster_rows = F,cluster_columns = F)

conserve_odonto <-  read.csv("process/trajectory/20250415_odonto_run_2/conservation_scores.csv",row.names = 1)
rownames(fitHM)%>%head

Heatmap(fitHM_odonto,km = 8,cluster_rows = F,cluster_columns = F)


# Assuming fitHM is your matrix/data frame for the heatmap
# Create the heatmap with k-means clustering (km=8)
set.seed(123)  # For reproducibility of clustering
ht_odonto <- Heatmap(fitHM_odonto,
              name = "Values",
              km = 8,  # k-means clustering with 8 clusters
              cluster_rows = F,
              cluster_columns = F,
              show_row_names = F,
              show_column_names = F,
              col=colorRamp2(c(-1, 0, 1), c("DeepSkyBlue3", "white", "red")))

# Draw the heatmap
ht_drawn <- draw(ht_odonto)

hmOrder <- row_order(ht_drawn)
hmGeneAll <- rownames(ht_drawn@ht_list$Values@matrix)
hmGeneListAll <- lapply(hmOrder, function(x) hmGeneAll[x])
library(stringi)
resAll <- as.data.frame((stri_list2matrix(hmGeneListAll)))
colnames(resAll) <- names(hmGeneListAll)
write.csv(resAll,"process/trajectory/20250415_odonto_run_2//hm_cluster_4000.csv")

# Reorder clusters according to specified order
geneOrder <- c("8","7","6","5","1","2","3","4")

# Create a flat vector of all indices in the desired order
ordered_indices <- unlist(hmOrder[geneOrder])

# Create a flat vector of all gene names in the desired order
ordered_genes <- hmGeneAll[ordered_indices]

# Create row split information for the new heatmap
# This maps each gene to its new cluster
cluster_map <- rep(NA, length(hmGeneAll))
for (i in seq_along(geneOrder)) {
  cluster_idx <- as.numeric(geneOrder[i])
  for (idx in hmOrder[[cluster_idx]]) {
    cluster_map[idx] <- i
  }
}

# Now create the reordered heatmap with row split
ht_reordered <- Heatmap(
  fitHM_odonto,
  name = "Values",
  cluster_rows = FALSE,  # No clustering since we're providing the order
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = colorRamp2(c(-1, 0, 1), c("DeepSkyBlue3", "white", "red")),
  row_order = ordered_indices,  # Use our ordered indices
  row_split = factor(cluster_map, levels = 1:8),  # Create the row split factor
  row_title = paste("GC", 1:8),
  row_gap = unit(2, "mm")
)

ht_reordered <- draw(ht_reordered)
pdf("results/trajectory/20250415_trajdtw_fit/20250417_odonto_hm_fit.pdf",width = 6,height = 8)
ht_reordered
dev.off()

genelist <- cbind(rownames(fitHM_odonto),cluster_map)
write.csv(genelist,"process/trajectory/20250415_odonto_run_2//20250417_genecluster.csv")


odontoGene <- c( hmGeneListAll[[2]], hmGeneListAll[[3]], hmGeneListAll[[4]])
conserve_odonto_filter <- conserve_odonto%>%filter(gene%in%odontoGene)
create_elbow_plot(na.omit(-conserve_odonto$raw_score))



