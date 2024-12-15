genelist <- list.files("processed_data/trajectory/model_fits/",full.names = T)
geneName <- list.files("processed_data/trajectory/model_fits/") %>% gsub("_fit.Rds","",.) %>% gsub("gene_","",.)
files <- lapply(genelist,readRDS)
weightMatList <- lapply(files, function(x) x$array_weight)
weightEstimate <- lapply(weightMatList, function(x) x$Estimate) %>% as.data.frame() %>% t()
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

# 3. Get row order

pdf("results/trajectory/20241129_conserved_trajectory_model/20241202_allgene_heatmap.pdf",width = 6,height = 6)
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
write.csv(ordered_data,"processed_data/trajectory/20241202_vargeneHM_order.csv")


#== visualize every cluster--------------------
plot_results_brms(files[[which(geneName%in%"Arl6ip1")]])
Heatmap(inputDataMat[,,which(genes%in%"Arl6ip1")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
plot_results_brms(files[[which(geneName%in%"Bmp4")]])
Heatmap(inputDataMat[,,which(genes%in%"Bmp4")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

#== cluster 4
plot_results_brms(files[[which(geneName%in%"Kcnc3")]])
Heatmap(inputDataMat[,,which(genes%in%"Kcnc3")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

plot_results_brms(files[[which(geneName%in%"Kif20a")]])
Heatmap(inputDataMat[,,which(genes%in%"Kif20a")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

plot_results_brms(files[[which(geneName%in%"Itih5")]])
Heatmap(inputDataMat[,,which(genes%in%"Itih5")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

plot_results_brms(files[[which(geneName%in%"Il34")]])
Heatmap(inputDataMat[,,which(genes%in%"Il34")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)


# cluster 8

plot_results_brms(files[[which(geneName%in%"Wnt11")]])
Heatmap(inputDataMat[,,which(genes%in%"Wnt11")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

plot_results_brms(files[[which(geneName%in%"Xpr1")]])
Heatmap(inputDataMat[,,which(genes%in%"Xpr1")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

# cluster 3


plot_results_brms(files[[which(geneName%in%"Ghr")]])
Heatmap(inputDataMat[,,which(genes%in%"Ghr")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

plot_results_brms(files[[which(geneName%in%"Gng11")]])
Heatmap(inputDataMat[,,which(genes%in%"Gng11")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)


gene_means <- apply(inputDataMat, 3, mean, na.rm = TRUE)
gene_var <- apply(inputDataMat, 3, sd, na.rm = TRUE)
weightMean <- weightEstimate %>% rowMeans()
plot(gene_means,weightMean)
plot(gene_var,weightMean)



plot(gene_means, weightMean,
     xlab = "Gene Mean Expression",
     ylab = "Weight Mean",
     main = "Gene Mean Expression vs Weights",
     pch = 16,
     col = rgb(0,0,0,0.5))
abline(lm(weightMean ~ gene_means), col = "red", lwd = 2)

# Add correlation coefficient
cor_mean <- cor(gene_means, weightMean, use = "complete.obs")
legend("topright",
       legend = sprintf("r = %.3f", cor_mean),
       bty = "n")

# For variance vs weights
plot(gene_var, weightMean,
     xlab = "Gene Standard Deviation",
     ylab = "Weight Mean",
     main = "Gene SD vs Weights",
     pch = 16,
     col = rgb(0,0,0,0.5))
abline(lm(weightMean ~ gene_var), col = "red", lwd = 2)

# Add correlation coefficient
cor_var <- cor(gene_var, weightMean, use = "complete.obs")
legend("topright",
       legend = sprintf("r = %.3f", cor_var),
       bty = "n")

# Create data frame
plot_data <- data.frame(
  mean = gene_means,
  sd = gene_var,
  weight = weightMean
)

# Mean vs weights plot
p1 <- ggplot(plot_data, aes(x = mean, y = weight)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  theme_minimal() +
  labs(
    title = "Gene Mean Expression vs Weights",
    x = "Gene Mean Expression",
    y = "Weight Mean"
  ) +
  annotate("text",
           x = max(plot_data$mean, na.rm = TRUE),
           y = max(plot_data$weight, na.rm = TRUE),
           label = sprintf("r = %.3f", cor_mean),
           hjust = 1)

# SD vs weights plot
p2 <- ggplot(plot_data, aes(x = sd, y = weight)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  theme_minimal() +
  labs(
    title = "Gene Standard Deviation vs Weights",
    x = "Gene Standard Deviation",
    y = "Weight Mean"
  ) +
  annotate("text",
           x = max(plot_data$sd, na.rm = TRUE),
           y = max(plot_data$weight, na.rm = TRUE),
           label = sprintf("r = %.3f", cor_var),
           hjust = 1)

# Arrange plots side by side
library(gridExtra)
grid.arrange(p1, p2, ncol = 2)
ggsave("results/trajectory/20241129_conserved_trajectory_model/20241201_gene_mean_vs_weights.pdf", p1, width = 6, height = 5)
ggsave("results/trajectory/20241129_conserved_trajectory_model/20241201_gene_sd_vs_weights.pdf", p2, width = 6, height = 5)

#== get residue-----------------------------------

# Function to calculate residuals
get_model_residuals <- function(fit, data) {
  # Get fitted values (posterior mean)
  fitted_values <- fitted(fit, summary = TRUE)[,"Estimate"]

  # Calculate residuals
  residuals_df <- data.frame(
    # Raw residuals
    raw = data$y - fitted_values,

    # Pearson residuals
    pearson = (data$y - fitted_values) / sqrt(fitted_values),

    # Response residuals
    response = residuals(fit, type = "ordinary"),

    # Additional information
    fitted = fitted_values,
    observed = data$y,
    array = data$array,
    x = data$x
  )

  return(residuals_df)
}


#== explore the relationship between residual and shape
resList <- list()
for (i in 1:length(genes)){
  preparedData <- prepare_data_for_gam(inputDataMat[,,i])
  x <- preparedData$x
  y <- preparedData$y
  array_idx <- preparedData$array_idx
  df <- data.frame(
    y = round(y[!is.na(y)]),  # ensure integers
    x = x[!is.na(y)],
    array = factor(array_idx[!is.na(y)])
  )
  fit = files[[i]]$fit
  resDf <- get_model_residuals(fit = fit,data = df)
  resList[[i]] <- resDf
}
fitted_values <- fitted(fit, summary = TRUE)[,"Estimate"]


# Function to calculate mean raw residuals by array
calculate_mean_raw_by_array <- function(residuals_df) {
  residuals_df %>%
    group_by(array) %>%
    summarise(mean_raw = mean(abs(pearson), na.rm = TRUE))
}

# Apply the function to each element in resList
mean_raw_by_array_list <- lapply(resList, calculate_mean_raw_by_array)
meanResMatrix <- lapply(mean_raw_by_array_list,function(x) x$mean_raw)%>%do.call(rbind, .)
meanResMatrix <- t(scale(t(meanResMatrix)))
weightEstimate_scale <- t(scale(t(weightEstimate)))
row_correlations <- numeric(nrow(weightEstimate))
for(i in 1:nrow(weightEstimate)) {
  row_correlations[i] <- cor(weightEstimate_scale[i,], meanResMatrix[i,],method = "spearman")
}

hist(row_correlations,breaks = 20)
plot(row_correlations)

weightEstimate%>%dim

which.min(row_correlations)

plot_results_brms(files[[which(geneName%in%"Adgrg6")]])
Heatmap(inputDataMat[,,which(genes%in%"Adgrg6")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

meanResMatrix[60,]


plot_results_brms(files[[which(geneName%in%"Gm26917")]])
Heatmap(inputDataMat[,,which(genes%in%"Gm26917")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
meanResMatrix[777,]

which.max(row_correlations)

plot_results_brms(files[[which(geneName%in%"Aqp1")]])
Heatmap(inputDataMat[,,which(genes%in%"Aqp1")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

plot(meanResMatrix[116,])

which(row_correlations> 0.9)

plot_results_brms(files[[which(geneName%in%"Abi3bp")]])
Heatmap(inputDataMat[,,which(genes%in%"Abi3bp")],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

plot(meanResMatrix[30,], type = "b")



arrayId <- resList[[1]]$array
rawResidue <- lapply(resList, function(x) x$raw)
rawResidue <- do.call(rbind,rawResidue)

weightEstimateDf <- as.data.frame(weightEstimate)
colnames(weightEstimateDf) <- files[[1]]$diagnostics$array

write.csv(weightEstimateDf,"processed_data/trajectory/20241203_run1_estimate.csv")
