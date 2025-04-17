# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(brms)
library(mgcv)
library(rstan)
library(dtw)  # Add DTW package

# Enable parallel execution for Stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Create directory for results if it doesn't exist
dir.create("results/trajectory", recursive = TRUE, showWarnings = FALSE)

# Source the utility functions - fix the paths
tryCatch({
  source("utils/trajectory_model_util.R")
  source("trajectory_model/20250406_benchmarking_utils.R")
}, error = function(e) {
  # Try alternative paths
  source("script/utils/trajectory_model_util.R")
  source("script/trajectory_model/20250406_benchmarking_utils.R")
})

# Function to format dataset for Bayesian model
formDataset <- function(dataMat, x) {
  x_values = rep(x, nrow(dataMat))
  y_values = dataMat %>% t() %>% as.vector()
  array_idx = rep(1:nrow(dataMat), each=ncol(dataMat))
  return(list(x_values, y_values, array_idx))
}

#---------------------------------------------------------
# Noise Response Benchmark - Test behavior with varying numbers of noise arrays
#---------------------------------------------------------

benchmark_noise_response <- function(
    n_simulations = 5,
    n_true_arrays = 10,
    noise_array_values = c(2, 5, 10, 20),  # Different numbers of noise arrays
    n_positions = 50,
    base_expression = 1.5,
    signal_strength = 2,
    signal_noise_sd = 0.3,
    noise_sd = 1.0  # Fixed noise level
) {
  # Initialize results storage
  noise_results <- list()

  # Run simulations with different numbers of noise arrays
  for (n_noise_arrays in noise_array_values) {
    cat(sprintf("\nRunning benchmark with %d noise arrays\n", n_noise_arrays))

    # Store results for this noise array count
    level_results <- list()

    for (sim in 1:n_simulations) {
      cat(sprintf("  Simulation %d/%d\n", sim, n_simulations))

      # Generate simulation data with fixed noise but varying number of noise arrays
      sim_data <- generate_simulation_data(
        n_true_arrays = n_true_arrays,
        n_noise_arrays = n_noise_arrays,
        n_positions = n_positions,
        base_expression = base_expression,
        signal_strength = signal_strength,
        signal_noise_sd = signal_noise_sd,
        noise_sd = noise_sd,
        seed = 1000 + sim
      )

      # Fit Bayesian model
      bayes_results <- fit_bayesian_model(sim_data, n_samples = 1000)

      # Calculate Pearson correlations
      pearson_results <- tryCatch({
        calculate_trajectory_correlations(
          trajectories = sim_data$data,
          method = "pearson",
          smooth_first = FALSE,
          benchmark = TRUE,
          true_signal = sim_data$true_signal,
          true_arrays = sim_data$true_arrays
        )
      }, error = function(e) {
        message("Error in calculating Pearson correlations: ", e$message)
        return(list(
          correlation_matrix = matrix(0, n_true_arrays + n_noise_arrays, n_true_arrays + n_noise_arrays),
          benchmark = list(
            correlation_with_truth = 0,
            mean_correlation_true_arrays = 0,
            mean_correlation_overall = 0
          )
        ))
      })

      # Calculate DTW similarities
      dtw_results <- tryCatch({
        calculate_trajectory_correlations(
          trajectories = sim_data$data,
          method = "dtw",
          smooth_first = FALSE,
          benchmark = TRUE,
          true_signal = sim_data$true_signal,
          true_arrays = sim_data$true_arrays
        )
      }, error = function(e) {
        message("Error in calculating DTW similarities: ", e$message)
        return(list(
          correlation_matrix = matrix(0, n_true_arrays + n_noise_arrays, n_true_arrays + n_noise_arrays),
          benchmark = list(
            correlation_with_truth = 0,
            mean_correlation_true_arrays = 0,
            mean_correlation_overall = 0
          )
        ))
      })

      # Get mean weight from Bayes and mean correlation from Pearson
      mean_bayes_weight <- mean(bayes_results$weights$Estimate)
      weight_for_true_arrays <- mean(bayes_results$weights$Estimate[sim_data$true_arrays])
      weight_for_noise_arrays <- mean(bayes_results$weights$Estimate[-sim_data$true_arrays])

      # Mean correlation values
      mean_pearson_cor <- mean(pearson_results$correlation_matrix[upper.tri(pearson_results$correlation_matrix)])
      mean_dtw_sim <- mean(dtw_results$correlation_matrix[upper.tri(dtw_results$correlation_matrix)])

      # Store results for this simulation
      level_results[[sim]] <- list(
        mean_bayes_weight = mean_bayes_weight,
        weight_for_true_arrays = weight_for_true_arrays,
        weight_for_noise_arrays = weight_for_noise_arrays,
        mean_pearson_cor = mean_pearson_cor,
        mean_dtw_sim = mean_dtw_sim,
        bayes_precision = bayes_results$precision,
        pearson_truth_cor = pearson_results$benchmark$correlation_with_truth,
        dtw_truth_cor = dtw_results$benchmark$correlation_with_truth,
        true_noise_ratio = n_true_arrays / n_noise_arrays
      )
    }

    # Average results across simulations for this noise array count
    avg_results <- list(
      n_noise_arrays = n_noise_arrays,
      true_noise_ratio = n_true_arrays / n_noise_arrays,
      mean_bayes_weight = mean(sapply(level_results, function(x) x$mean_bayes_weight)),
      weight_for_true_arrays = mean(sapply(level_results, function(x) x$weight_for_true_arrays)),
      weight_for_noise_arrays = mean(sapply(level_results, function(x) x$weight_for_noise_arrays)),
      mean_pearson_cor = mean(sapply(level_results, function(x) x$mean_pearson_cor)),
      mean_dtw_sim = mean(sapply(level_results, function(x) x$mean_dtw_sim)),
      bayes_precision = mean(sapply(level_results, function(x) x$bayes_precision)),
      pearson_truth_cor = mean(sapply(level_results, function(x) x$pearson_truth_cor)),
      dtw_truth_cor = mean(sapply(level_results, function(x) x$dtw_truth_cor)),
      simulations = level_results
    )

    noise_results[[as.character(n_noise_arrays)]] <- avg_results
  }

  # Create summary dataframe
  summary_df <- do.call(rbind, lapply(noise_results, function(result) {
    data.frame(
      n_noise_arrays = result$n_noise_arrays,
      true_noise_ratio = result$true_noise_ratio,
      mean_bayes_weight = result$mean_bayes_weight,
      weight_true_arrays = result$weight_for_true_arrays,
      weight_noise_arrays = result$weight_for_noise_arrays,
      mean_pearson_cor = result$mean_pearson_cor,
      mean_dtw_sim = result$mean_dtw_sim,
      bayes_precision = result$bayes_precision,
      pearson_truth_cor = result$pearson_truth_cor,
      dtw_truth_cor = result$dtw_truth_cor
    )
  }))

  # Create visualization of noise array response
  p_weights <- ggplot(summary_df, aes(x = n_noise_arrays)) +
    geom_line(aes(y = weight_true_arrays, color = "True Arrays"), size = 1) +
    geom_line(aes(y = weight_noise_arrays, color = "Noise Arrays"), size = 1) +
    scale_x_log10() +
    labs(title = "Bayesian Weight Response to Number of Noise Arrays",
         x = "Number of Noise Arrays (log scale)",
         y = "Mean Array Weight",
         color = "Array Type") +
    theme_minimal()

  p_correlations <- ggplot(summary_df, aes(x = n_noise_arrays)) +
    geom_line(aes(y = mean_pearson_cor, color = "Pearson"), size = 1) +
    geom_line(aes(y = mean_dtw_sim, color = "DTW"), size = 1) +
    scale_x_log10() +
    labs(title = "Similarity Measures Response to Number of Noise Arrays",
         x = "Number of Noise Arrays (log scale)",
         y = "Mean Similarity",
         color = "Method") +
    theme_minimal()

  p_combined <- ggplot(summary_df, aes(x = n_noise_arrays)) +
    geom_line(aes(y = bayes_precision, color = "Bayes Precision"), size = 1) +
    geom_line(aes(y = pearson_truth_cor, color = "Pearson-Truth Correlation"), size = 1) +
    geom_line(aes(y = dtw_truth_cor, color = "DTW-Truth Correlation"), size = 1) +
    scale_x_log10() +
    labs(title = "Performance Metrics vs Number of Noise Arrays",
         x = "Number of Noise Arrays (log scale)",
         y = "Metric Value",
         color = "Metric") +
    theme_minimal()

  # Create a plot showing performance vs true/noise ratio
  p_ratio <- ggplot(summary_df, aes(x = true_noise_ratio)) +
    geom_line(aes(y = bayes_precision, color = "Bayes Precision"), size = 1) +
    geom_line(aes(y = pearson_truth_cor, color = "Pearson-Truth Correlation"), size = 1) +
    geom_line(aes(y = dtw_truth_cor, color = "DTW-Truth Correlation"), size = 1) +
    labs(title = "Performance Metrics vs True/Noise Array Ratio",
         x = "True/Noise Array Ratio",
         y = "Metric Value",
         color = "Metric") +
    theme_minimal() +
    scale_x_continuous(trans = "log10")

  # Save plots
  ggsave("results/trajectory/noise_arrays_weights.pdf", p_weights, width = 8, height = 6)
  ggsave("results/trajectory/noise_arrays_correlations.pdf", p_correlations, width = 8, height = 6)
  ggsave("results/trajectory/noise_arrays_combined.pdf", p_combined, width = 8, height = 6)
  ggsave("results/trajectory/noise_arrays_ratio.pdf", p_ratio, width = 8, height = 6)

  # Combined visualization
  combined_plot <- (p_weights / p_correlations) | (p_combined / p_ratio)
  ggsave("results/trajectory/noise_arrays_summary.pdf", combined_plot, width = 16, height = 12)

  return(list(
    results = noise_results,
    summary = summary_df,
    plots = list(
      weights = p_weights,
      correlations = p_correlations,
      combined = p_combined,
      ratio = p_ratio
    )
  ))
}

#---------------------------------------------------------
# Test with a small simulation before running the full benchmark
#---------------------------------------------------------

# Test a single, small simulation first
cat("Running test simulation...\n")
test_simulation <- tryCatch({
  test_data <- generate_simulation_data(
    n_true_arrays = 10,
    n_noise_arrays = 3,
    n_positions = 60,
    base_expression = 1.5,
    signal_strength = 10,
    signal_noise_sd = 0.5,
    noise_sd = 1,
    seed = 42
  )

  # Create a simple plot to test
  pdf("results/trajectory/test_signal.pdf", width = 8, height = 4)
  plot(test_data$x, test_data$true_signal, type = "l", col = "red", lwd = 2,
       xlab = "Position", ylab = "Signal", main = "Test Signal Pattern")
  dev.off()

  cat("Test simulation successful!\n")
  TRUE
}, error = function(e) {
  cat("Test simulation failed:", e$message, "\n")
  FALSE
})

# Only run the benchmark if the test was successful
if (test_simulation) {
  # Run the benchmarking with reduced parameters
  cat("Running benchmark with reduced parameters...\n")

  benchmark_results <- tryCatch({
    benchmark_algorithms(
      n_simulations = 2,  # Reduced for demonstration
      n_true_arrays = 7,  # Smaller arrays
      n_noise_arrays = 3,
      n_positions = 50,    # Fewer positions
      base_expression = 1.5,
      signal_strength = 7,
      signal_noise_sd_values = c(0.7),
      noise_sd_values = c(1.5,3,5,7),
      correlation_methods = c("pearson", "dtw"),  # Include DTW in the test
      smooth_options = c(FALSE)          # Only test without smoothing for speed
    )
  }, error = function(e) {
    cat("Benchmark failed:", e$message, "\n")
    return(NULL)
  })

  # Save the results if benchmark was successful
  if (!is.null(benchmark_results)) {
    saveRDS(benchmark_results, "results/trajectory/benchmark_results.rds")

    #---------------------------------------------------------
    # Create summary visualizations
    #---------------------------------------------------------

    # Check if we have results to plot
    if (nrow(benchmark_results$summary) > 0) {
      # Average performance across all simulations
      summary_plot_data <- benchmark_results$summary %>%
        group_by(method, smooth, signal_noise_sd, noise_sd) %>%
        summarize(
          bayes_precision_mean = mean(bayes_precision),
          bayes_correlation_mean = mean(bayes_correlation),
          bayes_runtime_mean = mean(bayes_runtime),
          cor_precision_mean = mean(cor_precision),
          cor_recall_mean = mean(cor_recall),
          cor_f1_score_mean = mean(cor_f1_score),
          cor_correlation_mean = mean(cor_correlation),
          cor_runtime_mean = mean(cor_runtime),
          .groups = "drop"
        )

      # Create plots

      # 1. F1 Score by method and noise level
      p1 <- ggplot(summary_plot_data, aes(x = interaction(method, smooth), y = cor_f1_score_mean,
                                         fill = factor(signal_noise_sd))) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~noise_sd, labeller = label_both) +
        labs(title = "F1 Score by Correlation Method and Noise Level",
             x = "Method & Smoothing",
             y = "F1 Score",
             fill = "Signal Noise SD") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      # 2. Comparison of Bayesian vs Correlation precision
      p2 <- summary_plot_data %>%
        mutate(method_smooth = paste(method, smooth)) %>%
        ggplot(aes(x = bayes_precision_mean, y = cor_precision_mean, color = method_smooth)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        facet_grid(signal_noise_sd ~ noise_sd, labeller = label_both) +
        labs(title = "Bayesian vs Correlation Method Precision",
             x = "Bayesian Model Precision",
             y = "Correlation Method Precision",
             color = "Method & Smooth") +
        theme_minimal() +
        theme(legend.position = "bottom")

      # 3. Runtime comparison
      p3 <- summary_plot_data %>%
        ggplot(aes(x = interaction(method, smooth), y = cor_runtime_mean)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_hline(aes(yintercept = mean(bayes_runtime_mean)), color = "red", linetype = "dashed") +
        labs(title = "Runtime Comparison (Bayesian model in red)",
             x = "Method & Smoothing",
             y = "Runtime (seconds)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      # Combine plots
      combined_plot <- (p1 / p2) | p3

      # Save plots
      ggsave("results/trajectory/benchmark_f1_score.pdf", p1, width = 8, height = 6)
      ggsave("results/trajectory/benchmark_precision.pdf", p2, width = 10, height = 8)
      ggsave("results/trajectory/benchmark_runtime.pdf", p3, width = 8, height = 6)
      ggsave("results/trajectory/benchmark_summary.pdf", combined_plot, width = 15, height = 10)

      cat("Benchmark results visualized and saved to results/trajectory/\n")
    } else {
      cat("No summary data available to plot.\n")
    }

    # Run noise response benchmark if main benchmark was successful
    if (!is.null(benchmark_results)) {
      cat("\n\nRunning noise response benchmark...\n")

      noise_benchmark_results <- tryCatch({
        benchmark_noise_response(
          n_simulations = 3,            # Small number of simulations for each noise level
          n_true_arrays = 8,            # Smaller arrays
          noise_array_values = c(2, 5, 10, 20, 40),  # Different numbers of noise arrays
          n_positions = 50,             # Fewer positions
          base_expression = 1.5,
          signal_strength = 5,          # Low signal strength to make noise impact more visible
          signal_noise_sd = 0.3,        # Fixed signal noise
          noise_sd = 1.0                # Fixed noise level
        )
      }, error = function(e) {
        cat("Noise benchmark failed:", e$message, "\n")
        return(NULL)
      })

      # Save noise benchmark results
      if (!is.null(noise_benchmark_results)) {
        saveRDS(noise_benchmark_results, "results/trajectory/noise_benchmark_results.rds")
        cat("Noise benchmark complete. Results saved to results/trajectory/noise_benchmark_results.rds\n")

        # Print a summary of findings
        cat("\nSummary of noise response:\n")
        print(noise_benchmark_results$summary[, c("n_noise_arrays", "bayes_precision", "pearson_truth_cor", "dtw_truth_cor")])

        # Print a message about DTW method
        cat("\nDTW method has been added to the benchmarking. This method may perform better with\n")
        cat("time-shifted patterns where traditional correlation methods struggle.\n")
        cat("The results now include DTW-based similarity metrics for comparison.\n")

        # Print message about the noise array benchmark
        cat("\nThe benchmark now tests how methods perform with increasing numbers of noise arrays\n")
        cat("while keeping a fixed number of true signal arrays. This simulates scenarios where\n")
        cat("the signal-to-noise ratio decreases as more noisy samples are added to the dataset.\n")
        cat("The true/noise ratio plot shows how performance correlates with the proportion of\n")
        cat("true signal arrays in the dataset.\n")
      }
    }
  }
} else {
  cat("Skipping benchmark due to test failure.\n")
}

