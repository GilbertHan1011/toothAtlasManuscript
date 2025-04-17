# Simple script to plot metrics from noise_benchmark_results
# ----------------------------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create output directory
dir.create("results/trajectory", recursive = TRUE, showWarnings = FALSE)

# Reshape the data for plotting
plot_data <- noise_benchmark_results$summary %>%
  select(n_noise_arrays, true_noise_ratio, 
         mean_bayes_weight, mean_pearson_cor, mean_dtw_sim) %>%
  pivot_longer(
    cols = c(mean_bayes_weight, mean_pearson_cor, mean_dtw_sim),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_name = case_when(
      metric == "mean_bayes_weight" ~ "Bayes Weight",
      metric == "mean_pearson_cor" ~ "Pearson Correlation",
      metric == "mean_dtw_sim" ~ "DTW Similarity"
    )
  )

# Create the plot
p <- ggplot(plot_data, aes(x = true_noise_ratio, y = value, color = metric_name)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(
    title = "Method Performance by True/Noise Ratio",
    subtitle = "Comparing Bayes Weight, Pearson Correlation, and DTW Similarity",
    x = "True/Noise Array Ratio (higher = more signal)",
    y = "Metric Value",
    color = "Metric"
  ) +
  scale_x_log10() +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  scale_color_manual(values = c(
    "Bayes Weight" = "darkblue",
    "Pearson Correlation" = "darkorange",
    "DTW Similarity" = "darkgreen"
  ))

# Create a multi-axis plot to handle the different scales better
# This creates three separate plots with aligned x-axes

# Metrics by ratio - separate y-axes
p1 <- ggplot(filter(plot_data, metric == "mean_bayes_weight"), 
            aes(x = true_noise_ratio, y = value)) +
  geom_line(size = 1.2, color = "darkblue") +
  geom_point(size = 3, color = "darkblue") +
  labs(
    title = "Bayes Weight",
    x = NULL,
    y = "Weight Value"
  ) +
  scale_x_log10() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    plot.title = element_text(face = "bold", color = "darkblue")
  )

p2 <- ggplot(filter(plot_data, metric == "mean_pearson_cor"), 
            aes(x = true_noise_ratio, y = value)) +
  geom_line(size = 1.2, color = "darkorange") +
  geom_point(size = 3, color = "darkorange") +
  labs(
    title = "Pearson Correlation",
    x = NULL,
    y = "Correlation Value"
  ) +
  scale_x_log10() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    plot.title = element_text(face = "bold", color = "darkorange")
  )

p3 <- ggplot(filter(plot_data, metric == "mean_dtw_sim"), 
            aes(x = true_noise_ratio, y = value)) +
  geom_line(size = 1.2, color = "darkgreen") +
  geom_point(size = 3, color = "darkgreen") +
  labs(
    title = "DTW Similarity",
    x = "True/Noise Array Ratio (log scale)",
    y = "Similarity Value"
  ) +
  scale_x_log10() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    plot.title = element_text(face = "bold", color = "darkgreen")
  )

# Combine the plots using patchwork if available
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  multi_plot <- p1 / p2 / p3 + 
    plot_layout(heights = c(1, 1, 1)) +
    plot_annotation(
      title = "Method Performance by True/Noise Ratio",
      subtitle = "Each metric shown on its own scale",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Save the multi-plot
  ggsave("results/trajectory/metrics_separate_scales.pdf", multi_plot, width = 10, height = 12)
}

# Save the single plot
ggsave("results/trajectory/metrics_line_plot.pdf", p, width = 10, height = 7)

# Also create a normalized version for direct comparison
norm_data <- noise_benchmark_results$summary %>%
  select(n_noise_arrays, true_noise_ratio, 
         mean_bayes_weight, mean_pearson_cor, mean_dtw_sim) %>%
  mutate(
    # Scale each metric to 0-1 range
    norm_bayes = (mean_bayes_weight - min(mean_bayes_weight)) / 
                (max(mean_bayes_weight) - min(mean_bayes_weight)),
    norm_pearson = (mean_pearson_cor - min(mean_pearson_cor)) / 
                   (max(mean_pearson_cor) - min(mean_pearson_cor)),
    norm_dtw = (mean_dtw_sim - min(mean_dtw_sim)) / 
              (max(mean_dtw_sim) - min(mean_dtw_sim))
  ) %>%
  pivot_longer(
    cols = c(norm_bayes, norm_pearson, norm_dtw),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_name = case_when(
      metric == "norm_bayes" ~ "Bayes Weight",
      metric == "norm_pearson" ~ "Pearson Correlation",
      metric == "norm_dtw" ~ "DTW Similarity"
    )
  )

# Create normalized plot
p_norm <- ggplot(norm_data, aes(x = true_noise_ratio, y = value, color = metric_name)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(
    title = "Normalized Method Performance by True/Noise Ratio",
    subtitle = "All metrics scaled to 0-1 range for direct comparison",
    x = "True/Noise Array Ratio (higher = more signal)",
    y = "Normalized Value (0-1)",
    color = "Metric"
  ) +
  scale_x_log10() +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  scale_color_manual(values = c(
    "Bayes Weight" = "darkblue",
    "Pearson Correlation" = "darkorange",
    "DTW Similarity" = "darkgreen"
  ))

# Save the normalized plot
ggsave("results/trajectory/metrics_normalized_plot.pdf", p_norm, width = 10, height = 7)

cat("Plots created successfully! Check the results/trajectory directory.\n") 