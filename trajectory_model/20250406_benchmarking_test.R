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


cor_results <-
  calculate_trajectory_correlations(
    trajectories = test_data$data,
    method = "pearson",
    weights = NULL,  # No longer using array weights since the function now creates its own
    smooth_first = F,
    benchmark = TRUE,
    true_signal = test_data$true_signal,
    true_arrays = test_data$true_arrays
  )

Heatmap(cor_results$correlation_matrix,cluster_rows = F,cluster_columns = F)


benchmark_results <- tryCatch({
  benchmark_algorithms(
    n_simulations = 1,  # Reduced for demonstration
    n_true_arrays = 7,  # Smaller arrays
    n_noise_arrays = 3,
    n_positions = 50,    # Fewer positions
    base_expression = 1.5,
    signal_strength = 3,
    signal_noise_sd_values = c(1,3,5,7),
    noise_sd_values = c(7),
    correlation_methods = c("pearson", "dtw"),  # Include DTW in the test
    smooth_options = c(FALSE)          # Only test without smoothing for speed
  )
}, error = function(e) {
  cat("Benchmark failed:", e$message, "\n")
  return(NULL)
})


benchmark_results2 <- tryCatch({
  benchmark_algorithms(
    n_simulations = 1,  # Reduced for demonstration
    n_true_arrays = 15,  # Smaller arrays
    n_noise_arrays = 5,
    n_positions = 50,    # Fewer positions
    base_expression = 1.5,
    signal_strength = 3,
    signal_noise_sd_values = c(1,3,5,7),
    noise_sd_values = c(7),
    correlation_methods = c("pearson", "dtw"),  # Include DTW in the test
    smooth_options = c(FALSE)          # Only test without smoothing for speed
  )
}, error = function(e) {
  cat("Benchmark failed:", e$message, "\n")
  return(NULL)
})

benchmark_results2 <- tryCatch({
  benchmark_algorithms(
    n_simulations = 1,  # Reduced for demonstration
    n_true_arrays = 15,  # Smaller arrays
    n_noise_arrays = 5,
    n_positions = 50,    # Fewer positions
    base_expression = 1.5,
    signal_strength = 6,
    signal_noise_sd_values = c(1,3,5,7),
    noise_sd_values = c(7),
    correlation_methods = c("pearson", "dtw"),  # Include DTW in the test
    smooth_options = c(FALSE)          # Only test without smoothing for speed
  )
}, error = function(e) {
  cat("Benchmark failed:", e$message, "\n")
  return(NULL)
})


noise_benchmark_results <- tryCatch({
  benchmark_noise_response(
    n_simulations = 1,            # Small number of simulations for each noise level
    n_true_arrays = 7,            # Smaller arrays
    noise_array_values = c(2, 5, 10, 20),  # Different numbers of noise arrays
    n_positions = 50,             # Fewer positions
    base_expression = 1.5,
    signal_strength = 7,          # Low signal strength to make noise impact more visible
    signal_noise_sd = 0.3,        # Fixed signal noise
    noise_sd = 3.0                # Fixed noise level
  )
}, error = function(e) {
  cat("Noise benchmark failed:", e$message, "\n")
  return(NULL)
})


benchmark_results3 <- tryCatch({
  benchmark_algorithms(
    n_simulations = 1,  # Reduced for demonstration
    n_true_arrays = 15,  # Smaller arrays
    n_noise_arrays = 5,
    n_positions = 50,    # Fewer positions
    base_expression = 1.5,
    signal_strength = 15,
    signal_noise_sd_values = c(1),
    noise_sd_values = c(7),
    correlation_methods = c("pearson", "dtw"),  # Include DTW in the test
    smooth_options = c(FALSE)          # Only test without smoothing for speed
  )
}, error = function(e) {
  cat("Benchmark failed:", e$message, "\n")
  return(NULL)
})
