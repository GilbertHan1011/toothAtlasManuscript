combined <- readRDS("processed_data/trajectory/20250106_all_combined.Rds")

genelist <- list.files("results/trajectory/20250107_bone_conserved/",full.names = T,pattern = "*.csv")
geneName <- list.files("results/trajectory/20250107_bone_conserved/",pattern = "*.csv") %>% gsub("_fit_weight.csv","",.) %>% gsub("gene_","",.)
files <- lapply(genelist,read.csv, row.names = 1)
weightEstimate <- lapply(files, function(x) x$Estimate) %>% as.data.frame() %>% t()
rownames(weightEstimate) <- geneName
ComplexHeatmap::Heatmap(weightEstimate,show_row_names = F)



n_clusters = 8
km <- kmeans(weightEstimate, centers = n_clusters)

# 2. Create heatmap with k-means clustering
ht <- ComplexHeatmap::Heatmap(
  weightEstimate,
  show_row_names = FALSE,
  row_split = km$cluster,  # split by k-means clusters
  cluster_rows = T,     # disable clustering
)
pdf("results/trajectory/20250107_bone_tooth_conserved_plot/20250107_heatmap_1st.pdf",width = 6,height = 10)
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
write.csv(ordered_data,"results/trajectory/20250107_bone_tooth_conserved_plot/20250107_vargeneHM_order.csv")
