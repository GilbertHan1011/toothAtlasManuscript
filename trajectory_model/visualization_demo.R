# Demonstration script for simulation and benchmark visualization
# --------------------------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(RColorBrewer)

# Source the utility functions
source("trajectory_model/20250406_benchmarking_utils.R")

# 1. Generate simulation data with known signal and noise
# --------------------------------------------------------
cat("Generating simulation data...\n")
sim_data <- generate_simulation_data(
  n_true_arrays = 12,      # Number of arrays containing the true signal
  n_noise_arrays = 8,      # Number of arrays containing only noise
  n_positions = 80,        # Number of positions in each array
  base_expression = 2,     # Base expression level
  signal_strength = 3,     # Strength of the true signal
  signal_noise_sd = 0.5,   # Amount of noise in true signal arrays
  noise_sd = 1.0,          # Amount of noise in noise-only arrays
  seed = 123               # For reproducibility
)

# 2. Visualize the simulation data
# --------------------------------------------------------
cat("Creating static visualization...\n")
viz_results <- visualize_simulation(
  sim_data = sim_data,
  output_file = "results/trajectory/simulation_visualization.pdf"
)

# 3. Run correlation analysis on the data
# --------------------------------------------------------
cat("Calculating correlations with different methods...\n")

# Pearson correlation
pearson_results <- calculate_trajectory_correlations(
  trajectories = sim_data$data,
  method = "pearson",
  smooth_first = FALSE,
  benchmark = TRUE,
  true_signal = sim_data$true_signal,
  true_arrays = sim_data$true_arrays
)

# Try DTW correlation (will be skipped if package not available)
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
  cat("DTW calculation failed:", e$message, "\n")
  NULL
})

# 4. Visualize simulation with correlation results
# --------------------------------------------------------
cat("Creating visualization with correlation results...\n")
viz_with_results <- visualize_simulation(
  sim_data = sim_data,
  results = list(
    correlation_matrix = pearson_results$correlation_matrix,
    benchmark = pearson_results$benchmark,
    bayes = list(
      # Fake Bayesian weights for demonstration
      weights = data.frame(
        Estimate = c(rep(0.8, length(sim_data$true_arrays)), 
                     rep(0.2, length(sim_data$noise_arrays))),
        weight_norm = rep(1, nrow(sim_data$data))
      ),
      precision = 0.9
    )
  ),
  output_file = "results/trajectory/simulation_with_correlations.pdf"
)

# 5. Try interactive visualization if plotly is available
# --------------------------------------------------------
tryCatch({
  cat("Creating interactive visualization...\n")
  interactive_viz <- create_interactive_viz(
    sim_data = sim_data,
    output_file = "results/trajectory/interactive_simulation.html"
  )
  cat("Interactive visualization created! Open in a web browser.\n")
}, error = function(e) {
  cat("Interactive visualization failed:", e$message, "\n")
})

# 6. Simulate benchmark results and visualize them
# --------------------------------------------------------
cat("Creating mock benchmark visualization...\n")

# Create mock benchmark summary data
mock_summary <- data.frame(
  simulation = rep(1:3, each = 8),
  signal_noise_sd = rep(c(0.5, 1.0), each = 4, times = 3),
  noise_sd = rep(c(1.0, 2.0), each = 2, times = 6),
  method = rep(c("pearson", "dtw"), times = 12),
  smooth = rep(c(FALSE, TRUE), times = 12),
  bayes_precision = runif(24, 0.6, 0.9),
  bayes_correlation = runif(24, 0.5, 0.8),
  bayes_runtime = runif(24, 5, 15),
  cor_precision = runif(24, 0.5, 0.95),
  cor_recall = runif(24, 0.5, 0.95),
  cor_f1_score = runif(24, 0.5, 0.95),
  cor_correlation = runif(24, 0.4, 0.9),
  cor_runtime = runif(24, 0.2, 2)
)

# Create mock benchmark results
mock_benchmark <- list(
  summary = mock_summary,
  detailed_results = list()  # Not needed for visualization
)

# Visualize benchmark results
benchmark_viz <- visualize_benchmark(
  benchmark_results = mock_benchmark,
  output_dir = "results/trajectory/mock_benchmark"
)

cat("\nVisualization demo complete! Check the results in the 'results/trajectory' directory.\n") 