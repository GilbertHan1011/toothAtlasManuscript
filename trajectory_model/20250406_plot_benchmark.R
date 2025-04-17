source("script/trajectory_model/20250406_benchmarking_utils.R")
source("script/utils/trajectory_model_util.R")
library(brms)
library(dtw)
# Function to format dataset for Bayesian model
formDataset <- function(dataMat, x) {
  x_values = rep(x, nrow(dataMat))
  y_values = dataMat %>% t() %>% as.vector()
  array_idx = rep(1:nrow(dataMat), each=ncol(dataMat))
  return(list(x_values, y_values, array_idx))
}

test_data <- generate_simulation_data(
  n_true_arrays = 7,
  n_noise_arrays = 3,
  n_positions = 50,
  base_expression = 1.5,
  signal_strength = 10,
  signal_noise_sd = 1,
  noise_sd = 7,
  seed = 42
)
test_data2 <- generate_simulation_data(
  n_true_arrays = 7,
  n_noise_arrays = 3,
  n_positions = 50,
  base_expression = 1.5,
  signal_strength = 10,
  signal_noise_sd = 1,
  noise_sd = 7,
  seed = 42
)
viz_results <- simple_visualize_signal(
  sim_data = test_data,show_plot = T,
  output_file = "results/trajectory/benchmark/simulation_visualization.pdf"
)
viz_results2 <- simple_visualize_signal(
  sim_data = test_data2,show_plot = T,
  output_file = "results/trajectory/benchmark/simulation_visualization_2.pdf"
)


benchmark_results_static1 <- tryCatch({
  benchmark_algorithms(
    n_simulations = 3,  # Reduced for demonstration
    n_true_arrays = 15,  # Smaller arrays
    n_noise_arrays = 5,
    n_positions = 100,
    base_expression = 1.5,
    signal_strength = 7,
    signal_noise_sd_values = c(0.1,0.2,0.4,0.6,1,1.5,2,3,4,8,10),
    noise_sd_values = c(7),
    correlation_methods = c("pearson", "dtw"),  # Include DTW in the test
    smooth_options = c(FALSE)          # Only test without smoothing for speed
  )
}, error = function(e) {
  cat("Benchmark failed:", e$message, "\n")
  return(NULL)
})

benchmark_results_static1$summary %>% View

library(tidyverse)

write.csv(benchmark_results_static1$summary,"results/trajectory/20250406_noise_ratio_benchmraking.csv")
# Function to prepare and plot correlation by noise ratio
plot_correlation_by_noise <- function(benchmark_results, output_file = NULL) {
  # Extract the summary data
  summary_data <- benchmark_results$summary

  # Calculate noise ratio (signal_noise_sd / noise_sd)
  # This represents the ratio of signal noise to overall noise
  processed_data <- summary_data %>%
    mutate(noise_ratio = signal_noise_sd / noise_sd) %>%
    # Create a column to identify correlation type for Bayes vs other methods
    mutate(
      bayes_cor = bayes_correlation,
      method_cor = cor_correlation
    ) %>%
    # Reshape to long format to make plotting easier
    pivot_longer(
      cols = c(bayes_cor, method_cor),
      names_to = "correlation_type",
      values_to = "correlation_value"
    ) %>%
    # Create clearer method names
    mutate(
      correlation_method = case_when(
        correlation_type == "bayes_cor" ~ "Bayes",
        correlation_type == "method_cor" ~ as.character(method)
      )
    ) %>%
    # Remove duplicate Bayes entries
    filter(!(correlation_method == "Bayes" & method != "pearson"))

  # Calculate mean correlation by method, noise ratio and signal noise
  summary_stats <- processed_data %>%
    group_by(signal_noise_sd, noise_ratio, correlation_method) %>%
    summarize(
      mean_correlation = mean(correlation_value, na.rm = TRUE),
      sd_correlation = sd(correlation_value, na.rm = TRUE),
      n = n(),
      se_correlation = sd_correlation / sqrt(n),
      .groups = "drop"
    )

  # Create the line plot
  p <- ggplot(summary_stats,
              aes(x = noise_ratio, y = mean_correlation,
                  color = correlation_method, group = correlation_method)) +
    # Add confidence interval ribbon
    geom_ribbon(
      aes(ymin = mean_correlation - se_correlation,
          ymax = mean_correlation + se_correlation,
          fill = correlation_method),
      alpha = 0.2,
      color = NA
    ) +
    # Add lines connecting points
    geom_line(size = 1) +
    # Add points at each noise ratio
    geom_point(size = 3) +
    # Set labels
    labs(
      title = "Correlation Performance by Noise Ratio",
      subtitle = "Comparing Bayes, Pearson, and DTW methods",
      x = "Noise Ratio (signal_noise_sd / noise_sd)",
      y = "Mean Correlation",
      color = "Method",
      fill = "Method"
    ) +
    # Use a minimal theme
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    # Set custom colors
    scale_color_manual(values = c("Bayes" = "darkblue", "pearson" = "darkorange", "dtw" = "darkgreen")) +
    scale_fill_manual(values = c("Bayes" = "darkblue", "pearson" = "darkorange", "dtw" = "darkgreen"))

  # Save the plot if output file is specified
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(output_file, p, width = 10, height = 7)
    message(paste("Correlation plot saved to", output_file))
  }

  # Return the plot object
  return(p)
}
plot_correlation_by_noise(benchmark_results_static1,output_file = "results/trajectory/benchmark/20250406_noise_level.pdf")



noise_benchmark_results <- tryCatch({
  benchmark_noise_response(
    n_simulations = 3,            # Small number of simulations for each noise level
    n_true_arrays = 10,            # Smaller arrays
    noise_array_values = c(1,2, 5, 8,10, 15, 20,40),  # Different numbers of noise arrays
    n_positions = 100,             # Fewer positions
    base_expression = 1.5,
    signal_strength = 7,          # Low signal strength to make noise impact more visible
    signal_noise_sd = 1,        # Fixed signal noise
    noise_sd = 5.0                # Fixed noise level
  )
}, error = function(e) {
  cat("Noise benchmark failed:", e$message, "\n")
  return(NULL)
})


# Function to create the plot
plot_metrics_by_noise_ratio <- function(data, output_file = NULL) {
  # Reshape the data from wide to long format for plotting
  long_data <- data %>%
    select(n_noise_arrays, true_noise_ratio,
           mean_bayes_weight, mean_pearson_cor, mean_dtw_sim) %>%
    mutate(true_noise_ratio = 1/true_noise_ratio)%>%
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
      title = "Normalized Metrics by Noise/True Array Ratio",
      subtitle = "Each metric normalized to 0-1 scale for easier comparison",
      x = "Noise/True Array Ratio (higher = low signal)",
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
#noise_benchmark_results$summary$true_noise_ratio <- 1/noise_benchmark_results$summary$true_noise_ratio
plot_metrics_by_noise_ratio(noise_benchmark_results$summary,output_file = "results/trajectory/benchmark/20250406_truenoise_weight.pdf")


ggsave("", width = 10, height = 7)
