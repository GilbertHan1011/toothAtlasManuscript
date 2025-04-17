# Script to plot metrics against true_noise_ratio
# --------------------------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create output directory
dir.create("results/trajectory", recursive = TRUE, showWarnings = FALSE)

# Function to create the plot
plot_metrics_by_noise_ratio <- function(data, output_file = NULL) {
  # Reshape the data from wide to long format for plotting
  long_data <- data %>%
    select(n_noise_arrays, true_noise_ratio, 
           mean_bayes_weight, mean_pearson_cor, mean_dtw_sim) %>%
    pivot_longer(
      cols = c(mean_bayes_weight, mean_pearson_cor, mean_dtw_sim),
      names_to = "metric",
      values_to = "value"
    ) %>%
    # Make metric names more readable
    mutate(
      metric_name = case_when(
        metric == "mean_bayes_weight" ~ "Bayes Weight",
        metric == "mean_pearson_cor" ~ "Pearson Correlation",
        metric == "mean_dtw_sim" ~ "DTW Similarity"
      )
    )
  
  # Create the line plot
  p <- ggplot(long_data, aes(x = true_noise_ratio, y = value, color = metric_name)) +
    # Add lines
    geom_line(size = 1.2) +
    # Add points at each data point
    geom_point(size = 3) +
    # Add labels
    labs(
      title = "Metrics by True/Noise Array Ratio",
      subtitle = "Relationship between metrics and signal-to-noise ratio",
      x = "True/Noise Array Ratio (higher = more signal)",
      y = "Metric Value",
      color = "Metric"
    ) +
    # Use a logarithmic scale for x-axis since ratios can span wide range
    scale_x_log10() +
    # Custom theme
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    # Set custom colors
    scale_color_manual(values = c(
      "Bayes Weight" = "darkblue",
      "Pearson Correlation" = "darkorange",
      "DTW Similarity" = "darkgreen"
    ))
  
  # Because the metrics have very different scales, we might want to create
  # a normalized version where each metric is scaled to 0-1 range
  
  # Normalize the data
  normalized_data <- data %>%
    select(n_noise_arrays, true_noise_ratio, 
           mean_bayes_weight, mean_pearson_cor, mean_dtw_sim) %>%
    mutate(
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
      values_to = "normalized_value"
    ) %>%
    # Make metric names more readable
    mutate(
      metric_name = case_when(
        metric == "norm_bayes" ~ "Bayes Weight",
        metric == "norm_pearson" ~ "Pearson Correlation",
        metric == "norm_dtw" ~ "DTW Similarity"
      )
    )
  
  # Create the normalized plot
  p_norm <- ggplot(normalized_data, 
                  aes(x = true_noise_ratio, y = normalized_value, color = metric_name)) +
    # Add lines
    geom_line(size = 1.2) +
    # Add points at each data point
    geom_point(size = 3) +
    # Add labels
    labs(
      title = "Normalized Metrics by True/Noise Array Ratio",
      subtitle = "Each metric normalized to 0-1 scale for easier comparison",
      x = "True/Noise Array Ratio (higher = more signal)",
      y = "Normalized Value (0-1)",
      color = "Metric"
    ) +
    # Use a logarithmic scale for x-axis
    scale_x_log10() +
    # Custom theme
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    # Set custom colors
    scale_color_manual(values = c(
      "Bayes Weight" = "darkblue",
      "Pearson Correlation" = "darkorange",
      "DTW Similarity" = "darkgreen"
    ))
  
  # Save the plots if output file is specified
  if (!is.null(output_file)) {
    # Save the raw values plot
    ggsave(
      paste0(dirname(output_file), "/raw_", basename(output_file)),
      p, width = 10, height = 7
    )
    
    # Save the normalized plot
    ggsave(
      paste0(dirname(output_file), "/norm_", basename(output_file)),
      p_norm, width = 10, height = 7
    )
    
    message(paste("Plots saved to", dirname(output_file)))
  }
  
  # Return both plots
  return(list(raw = p, normalized = p_norm))
}

# Example usage:
# If you have the noise_benchmark_results loaded in your environment:
# plots <- plot_metrics_by_noise_ratio(
#   noise_benchmark_results$summary,
#   "results/trajectory/metrics_by_noise_ratio.pdf"
# )

# Usage with the sample data provided
sample_data <- data.frame(
  n_noise_arrays = c(1, 2, 5, 8, 10, 15),
  true_noise_ratio = c(10.0000000, 5.0000000, 2.0000000, 1.2500000, 1.0000000, 0.6666667),
  mean_bayes_weight = c(3.190578, 3.194538, 3.197551, 3.188779, 3.182674, 3.182810),
  mean_pearson_cor = c(0.23109440, 0.19626046, 0.12391574, 0.08671810, 0.07013007, 0.04726964),
  mean_dtw_sim = c(0.6038562, 0.5264506, 0.4097408, 0.3523004, 0.3299454, 0.3056311),
  weight_true_arrays = c(3.192397, 3.202184, 3.209877, 3.194804, 3.186401, 3.188658),
  weight_noise_arrays = c(3.172386, 3.156308, 3.172899, 3.181248, 3.178947, 3.178911),
  bayes_precision = c(0.9333333, 0.9000000, 0.7000000, 0.6333333, 0.5666667, 0.4333333),
  pearson_truth_cor = c(0.5000000, 0.6477503, 0.8183171, 0.8619942, 0.8671100, 0.8492078),
  dtw_truth_cor = c(0.5000000, 0.6477503, 0.8183171, 0.8619942, 0.8671100, 0.8492078)
)

# Create the plots with the sample data
plots <- plot_metrics_by_noise_ratio(
  sample_data,
  "results/trajectory/metrics_by_noise_ratio.pdf"
)

# Display the plots (works in RStudio or interactive R)
if (interactive()) {
  print(plots$raw)
  print(plots$normalized)
}

cat("Plots created successfully! Check the results/trajectory directory.\n") 