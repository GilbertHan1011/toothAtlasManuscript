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
odontoGene <- c( hmGeneListAll[[2]], hmGeneListAll[[3]], hmGeneListAll[[4]])
conserve_odonto_filter <- conserve_odonto%>%filter(gene%in%odontoGene)
create_elbow_plot(na.omit(-conserve_odonto$raw_score))
