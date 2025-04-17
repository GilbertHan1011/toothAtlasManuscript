# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(mgcv)
library(rstan)
library(brms)

# Enable parallel execution for Stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Source the utility functions
source("utils/trajectory_model_util.R")
source("trajectory_model/20250406_benchmarking_utils.R")

# Create directory for results if it doesn't exist
dir.create("results/trajectory/benchmark", recursive = TRUE, showWarnings = FALSE)

# Function to format dataset for Bayesian model
formDataset <- function(dataMat, x) {
  x_values = rep(x, nrow(dataMat))
  y_values = dataMat %>% t() %>% as.vector()
  array_idx = rep(1:nrow(dataMat), each=ncol(dataMat))
  return(list(x_values, y_values, array_idx))
}

#---------------------------------------------------------
# Run a single simulation for testing
#---------------------------------------------------------

test_single_simulation <- function() {
  # Generate simulation data
  cat("Generating test simulation data...\n")
  sim_data <- generate_simulation_data(
    n_true_arrays = 10,  # Reduced for testing
    n_noise_arrays = 3,  # Reduced for testing
    n_positions = 50,    # Reduced for testing
    base_expression = 1.5,
    signal_strength = 3,
    signal_noise_sd = 0.5,
    noise_sd = 1,
    seed = 42
  )
  
  # Plot the true signal
  pdf("results/trajectory/benchmark/test_true_signal.pdf", width = 8, height = 4)
  plot(sim_data$x, sim_data$true_signal, type = "l", col = "red", lwd = 2,
       xlab = "Position", ylab = "Signal", main = "True Signal Pattern")
  dev.off()
  
  # Create a heatmap of the data
  pdf("results/trajectory/benchmark/test_data_heatmap.pdf", width = 8, height = 6)
  trueCol <- c(rep("true", length(sim_data$true_arrays)), 
               rep("false", length(sim_data$noise_arrays)))
  Heatmap(sim_data$data, cluster_rows = FALSE, cluster_columns = FALSE,
          show_column_names = FALSE, show_row_names = TRUE,
          col = colorRamp2(c(0, 2, 4), c("blue", "white", "red")),
          split = trueCol, border = TRUE,
          column_title = "Simulated Data Patterns")
  dev.off()
  
  # Test the correlation method
  cat("Testing correlation method...\n")
  cor_results <- calculate_trajectory_correlations(
    trajectories = sim_data$data,
    method = "pearson",
    smooth_first = TRUE,
    benchmark = TRUE,
    true_signal = sim_data$true_signal,
    true_arrays = sim_data$true_arrays
  )
  
  # Visualize correlation matrix
  pdf("results/trajectory/benchmark/test_correlation_matrix.pdf", width = 8, height = 7)
  Heatmap(
    cor_results$correlation_matrix,
    name = "Correlation",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_split = factor(trueCol, levels = c("true", "false")),
    column_split = factor(trueCol, levels = c("true", "false")),
    column_title = "Pairwise Correlations (Pearson)"
  )
  dev.off()
  
  # Test the Bayesian model
  cat("Testing Bayesian model...\n")
  
  # Prepare data for Bayesian GAM
  input_data <- formDataset(sim_data$data, sim_data$x)
  
  # Use a try-catch to handle potential errors
  bayes_result <- tryCatch({
    # Fit Bayesian model with reduced samples for testing
    fit <- bayesian_gam_regression_nb_shape(
      x = input_data[[2]],
      y = input_data[[1]],
      array_idx = input_data[[3]],
      n_samples = 500  # Reduced for testing
    )
    
    # Extract array weights
    weights <- fit$array_weights
    
    # Plot the weights
    pdf("results/trajectory/benchmark/test_bayes_weights.pdf", width = 8, height = 4)
    plot(weights$Estimate, type = "h", lwd = 3, col = ifelse(1:nrow(weights) %in% sim_data$true_arrays, "blue", "red"),
         xlab = "Array", ylab = "Weight Estimate", main = "Bayesian Model Array Weights")
    legend("topright", legend = c("True Arrays", "Noise Arrays"), 
           col = c("blue", "red"), lwd = 3)
    dev.off()
    
    return(fit)
  }, error = function(e) {
    cat("Error in Bayesian model:", e$message, "\n")
    return(NULL)
  })
  
  return(list(
    sim_data = sim_data,
    correlation_results = cor_results,
    bayesian_results = bayes_result
  ))
}

#---------------------------------------------------------
# Run the benchmarking with reduced parameters
#---------------------------------------------------------

run_benchmark <- function() {
  cat("Starting benchmarking process...\n")
  
  # Use smaller parameters for faster execution
  benchmark_results <- benchmark_algorithms(
    n_simulations = 2,  # Just 2 simulations per configuration
    n_true_arrays = 10,  # Reduced number of arrays
    n_noise_arrays = 3,
    n_positions = 50,    # Fewer positions
    base_expression = 1.5,
    signal_strength = 3,
    signal_noise_sd_values = c(0.3, 0.7),  # Just two noise levels
    noise_sd_values = c(0.5, 1.5),        # Just two noise levels
    correlation_methods = c("pearson", "weighted"),
    smooth_options = c(TRUE, FALSE)
  )
  
  # Save the results
  saveRDS(benchmark_results, "results/trajectory/benchmark/benchmark_results.rds")
  
  # Create summary plots
  create_summary_plots(benchmark_results)
  
  return(benchmark_results)
}

#---------------------------------------------------------
# Create summary visualizations
#---------------------------------------------------------

create_summary_plots <- function(benchmark_results) {
  # Extract summary data
  summary_data <- benchmark_results$summary
  
  # 1. F1 Score by method and noise level
  p1 <- ggplot(summary_data, aes(x = interaction(method, smooth), y = cor_f1_score, 
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
  p2 <- summary_data %>%
    mutate(method_smooth = paste(method, smooth)) %>%
    ggplot(aes(x = bayes_precision, y = cor_precision, color = method_smooth)) +
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
  p3 <- summary_data %>%
    group_by(method, smooth) %>%
    summarize(
      mean_cor_runtime = mean(cor_runtime),
      mean_bayes_runtime = mean(bayes_runtime),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = interaction(method, smooth), y = mean_cor_runtime)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(aes(yintercept = mean(mean_bayes_runtime)), color = "red", linetype = "dashed") +
    labs(title = "Runtime Comparison (Bayesian model in red)",
         x = "Method & Smoothing",
         y = "Runtime (seconds)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine plots
  combined_plot <- (p1 / p2) | p3
  
  # Save the combined plot
  ggsave("results/trajectory/benchmark/summary_plots.pdf", combined_plot, width = 15, height = 10)
  
  # Individual plots
  ggsave("results/trajectory/benchmark/f1_score_plot.pdf", p1, width = 8, height = 6)
  ggsave("results/trajectory/benchmark/precision_comparison_plot.pdf", p2, width = 10, height = 8)
  ggsave("results/trajectory/benchmark/runtime_comparison_plot.pdf", p3, width = 8, height = 6)
}

#---------------------------------------------------------
# Run the benchmarking
#---------------------------------------------------------

cat("Running test simulation...\n")
test_results <- test_single_simulation()

# Check if the test was successful
if (!is.null(test_results)) {
  cat("Test successful! Running full benchmark...\n")
  
  # Here you can choose whether to run the full benchmark
  # This is commented out by default as it may take a long time
  # benchmark_results <- run_benchmark()
  
  cat("You can run the full benchmark by uncommenting the run_benchmark() line in the script.\n")
} else {
  cat("Test failed. Please review the error messages.\n")
} 