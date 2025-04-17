fitHM <- read.csv("process/trajectory/20250414_epi_run_5//fitted_trajectories.csv",row.names = 1)
fitHM <- t(fitHM)
Heatmap(fitHM[c("Ambn","Amtn","Enam","Mcm3","Gpr45","Rpl7"),],cluster_rows = F,cluster_columns = F)
rownames(fitHM)%>%head
conserved <- read.csv("process/trajectory/20250414_epi_run_5/conservation_scores.csv",row.names = 1)
conserved2 <- read.csv("process/trajectory/20250414_epi_run_5/conservation_scores.csv",row.names = 1)

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
write.csv(resAll,"process/trajectory/20250414_epi_run_5/hm_cluster_5000.csv")

cluster1 <- hmGeneListAll[[1]]
cluster2 <- hmGeneListAll[[2]]
cluster3 <- hmGeneListAll[[3]]

conservation_score = read.csv("process/trajectory/20250414_epi_run_5/conservation_scores.csv",row.names = 1)
conservation_score_filter <- conservation_score%>%
  filter(gene %in%c(cluster1,cluster2,cluster3))

top500  <- conservation_score_filter$gene%>%head(500)

de_filter<- de_gene[!rownames(de_gene)%in%top500,]
de_gene <- de_gene%>%rownames_to_column("gene")
de_gene_join<- de_gene%>% left_join(conservation_score)

conserved_run4 <- read.csv("process/trajectory/20250414_epi_run_4/conservation_scores.csv",row.names = 1)
de_gene_join_run4 <- de_gene%>% left_join(conserved_run4)


plot(de_gene_join_run4$avg_log2FC, -de_gene_join_run4$raw_score)
de_gene_join_run4 <- de_gene_join_run4%>%filter()
# Create the scatter plot
ggplot(de_gene_join_run4, aes(x = avg_log2FC, y = -raw_score)) +
  geom_point(alpha = 0.7, size = 2, color = "DeepSkyBlue3") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  labs(
    title = "Relationship between DE and Conservation",
    x = "Average log2 Fold Change",
    y = "Conservation Score (-raw_score)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )


# First, let's calculate the counts of points in each region
# Filter points inside the box (x from 3 to 6, y from 0 to 100)
in_box <- de_gene_join_run4 %>%
  filter(avg_log2FC >= 1.5 & avg_log2FC <= 8 &
           -raw_score >= 0 & -raw_score <= 100)

line_thred = 29.5
# Count points above and below y = 20 within the box
above_line <- sum(-in_box$raw_score > line_thred)
below_line <- sum(-in_box$raw_score <= line_thred)

# Create a data frame for the count labels
count_labels <- data.frame(
  x = c(4.5, 4.5),
  y = c(60, 10),
  label = c(paste0("Above: ", above_line), paste0("Below: ", below_line))
)

# Now create the plot with all elements
ggplot(de_gene_join_run4, aes(x = avg_log2FC, y = -raw_score)) +
  # Original points
  geom_point(alpha = 0.7, size = 2, color = "DeepSkyBlue3") +

  # Add red box from x=(3,6) and y=(0,100)
  geom_rect(aes(xmin = 2, xmax = 8, ymin = 0, ymax = 100),
            fill = NA, color = "red", linewidth = 1) +

  # Add red line at y = 20
  geom_hline(yintercept = 20, color = "red", linewidth = 1) +

  # Original reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +

  # Add count labels
  geom_text(data = count_labels, aes(x = x, y = y, label = label),
            color = "red", fontface = "bold", size = 5) +

  # Rest of the original styling
  labs(
    title = "Relationship between DE and Conservation",
    subtitle = paste("Points in box:", nrow(in_box),
                     "| Above y=20:", above_line,
                     "| Below y=20:", below_line),
    x = "Average log2 Fold Change",
    y = "Conservation Score (-raw_score)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "darkred"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  # Set axis limits to ensure the entire box is visible
  coord_cartesian(xlim = c(min(de_gene_join_run4$avg_log2FC, 0),
                           max(de_gene_join_run4$avg_log2FC, 7)),
                  ylim = c(min(-de_gene_join_run4$raw_score, 0),
                           max(-de_gene_join_run4$raw_score, 105)))
ggsave("results/trajectory/20250415_trajdtw_fit/20250416_epi_benchmarking.pdf",width = 8,height = 6)

gene_above = in_box$gene[-in_box$raw_score > line_thred]
gene_below = in_box$gene[-in_box$raw_score <= line_thred]

bm_list <- list(conservate = gene_below, unconservate = gene_above)
bm_ck_res <- clusterProfiler::compareCluster(bm_list, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(bm_ck_res)





# First, let's calculate the counts of points in each region
# Filter points inside the box (x from 3 to 6, y from 0 to 100)
in_box <- de_gene_join_run4 %>%
  filter(avg_log2FC >= 1 & avg_log2FC <= 8 &
           -raw_score >= 0 & -raw_score <= 100)

line_thred = 33
# Count points above and below y = 20 within the box
above_line <- sum(-in_box$raw_score > line_thred)
below_line <- sum(-in_box$raw_score <= line_thred)

# Create a data frame for the count labels
count_labels <- data.frame(
  x = c(4.5, 4.5),
  y = c(60, 10),
  label = c(paste0("Above: ", above_line), paste0("Below: ", below_line))
)

# Now create the plot with all elements
ggplot(de_gene_join_run4, aes(x = avg_log2FC, y = -raw_score)) +
  # Original points
  geom_point(alpha = 0.7, size = 2, color = "DeepSkyBlue3") +

  # Add red box from x=(3,6) and y=(0,100)
  geom_rect(aes(xmin = 2, xmax = 8, ymin = 0, ymax = 100),
            fill = NA, color = "red", linewidth = 1) +

  # Add red line at y = 20
  geom_hline(yintercept = 20, color = "red", linewidth = 1) +

  # Original reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +

  # Add count labels
  geom_text(data = count_labels, aes(x = x, y = y, label = label),
            color = "red", fontface = "bold", size = 5) +

  # Rest of the original styling
  labs(
    title = "Relationship between DE and Conservation",
    subtitle = paste("Points in box:", nrow(in_box),
                     "| Above y=20:", above_line,
                     "| Below y=20:", below_line),
    x = "Average log2 Fold Change",
    y = "Conservation Score (-raw_score)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "darkred"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  # Set axis limits to ensure the entire box is visible
  coord_cartesian(xlim = c(min(de_gene_join_run4$avg_log2FC, 0),
                           max(de_gene_join_run4$avg_log2FC, 7)),
                  ylim = c(min(-de_gene_join_run4$raw_score, 0),
                           max(-de_gene_join_run4$raw_score, 105)))
ggsave("results/trajectory/20250415_trajdtw_fit/20250416_epi_benchmarking.pdf",width = 8,height = 6)

gene_above = in_box$gene[-in_box$raw_score > line_thred]
gene_below = in_box$gene[-in_box$raw_score <= line_thred]

bm_list <- list(conservate = gene_below, unconservate = gene_above)
bm_ck_res2 <- clusterProfiler::compareCluster(bm_list, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(bm_ck_res2)
