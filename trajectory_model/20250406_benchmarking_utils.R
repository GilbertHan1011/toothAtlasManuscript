calculate_trajectory_correlations <- function(
    trajectories,                # Can be the simulation data matrix
    method = c("pearson", "spearman", "kendall", "weighted", "dtw"),
    weights = NULL,              # Use array weights from Bayesian model if weighted
    smooth_first = FALSE,        # Option to smooth trajectories before correlation
    n_knots = 5,                 # Number of knots if smoothing
    benchmark = FALSE,           # Whether to run benchmarking
    true_signal = NULL,          # True signal for benchmarking
    true_arrays = NULL           # Indices of true arrays for benchmarking
) {

  # Fix: Don't use match.arg() which requires exact matches in the default values
  # Instead, validate the method manually to allow "dtw" as a valid method
  valid_methods <- c("pearson", "spearman", "kendall", "weighted", "dtw")
  if (!method %in% valid_methods) {
    stop(sprintf("'arg' should be one of %s", paste(shQuote(valid_methods), collapse = ", ")))
  }

  n_trajectories <- nrow(trajectories)

  # Initialize correlation matrix
  cor_matrix <- matrix(
    NA,
    nrow = n_trajectories,
    ncol = n_trajectories,
    dimnames = list(rownames(trajectories), rownames(trajectories))
  )

  # Smooth trajectories if requested
  if (smooth_first) {
    smoothed_trajectories <- matrix(
      NA,
      nrow = n_trajectories,
      ncol = ncol(trajectories),
      dimnames = dimnames(trajectories)
    )

    for (i in 1:n_trajectories) {
      # Use GAM to smooth each trajectory
      x_seq <- seq(0, 1, length.out = ncol(trajectories))
      gam_fit <- mgcv::gam(trajectories[i,] ~ s(x_seq, bs = "cr", k = n_knots))
      smoothed_trajectories[i,] <- predict(gam_fit)
    }
    working_trajectories <- smoothed_trajectories
  } else {
    working_trajectories <- trajectories
  }

  # Check if weights are provided for weighted correlation
  if (method == "weighted" && !is.null(weights)) {
    # Ensure weights match the number of trajectories
    if (length(weights) != n_trajectories) {
      warning("Number of weights does not match number of trajectories. Using equal weights.")
      weights <- rep(1, n_trajectories)
    }
  }

  # For DTW, use the dedicated function
  if (method == "dtw") {
    # Check if dtw package is available
    if (!requireNamespace("dtw", quietly = TRUE)) {
      stop("Package 'dtw' is required for DTW similarity. Please install it.")
    }

    # Use the DTW similarity function
    cor_matrix <- calculate_dtw_similarity(working_trajectories)
  } else {
    # Calculate pairwise correlations for other methods
    for (i in 1:n_trajectories) {
      for (j in i:n_trajectories) {
        if (i == j) {
          cor_matrix[i, j] <- 1
        } else {
          if (method == "weighted" && !is.null(weights)) {
            # For weighted correlation, use weights from Bayesian model
            # Create weights for this specific correlation calculation
            cor_weights <- rep(1, ncol(working_trajectories))

            cor_val <- calculate_weighted_correlation(
              working_trajectories[i,],
              working_trajectories[j,],
              weights = cor_weights  # Use a simple weight vector matching the length of the data
            )
          } else {
            # Standard correlation methods
            cor_val <- cor(
              working_trajectories[i,],
              working_trajectories[j,],
              method = method,
              use = "pairwise.complete.obs"
            )
          }
          cor_matrix[i, j] <- cor_matrix[j, i] <- cor_val
        }
      }
    }
  }

  # Benchmarking if requested
  if (benchmark && !is.null(true_signal) && !is.null(true_arrays)) {
    benchmark_results <- run_benchmarking(
      cor_matrix = cor_matrix,
      trajectories = working_trajectories,
      true_signal = true_signal,
      true_arrays = true_arrays
    )

    return(list(
      correlation_matrix = cor_matrix,
      benchmark = benchmark_results
    ))
  }

  return(cor_matrix)
}

# Helper function for weighted correlation
calculate_weighted_correlation <- function(x, y, weights) {
  # Make sure none of the inputs are NULL
  if (is.null(x) || is.null(y) || is.null(weights)) {
    return(NA)
  }

  # Handle NA values
  valid_indices <- !is.na(x) & !is.na(y)
  if (sum(valid_indices) < 2) {
    return(NA)  # Can't calculate correlation with less than 2 points
  }

  x <- x[valid_indices]
  y <- y[valid_indices]

  # Make sure weights match the length of x and y after NA removal
  if (length(weights) != length(x)) {
    # If weights are per array (e.g., one per row), repeat weight values
    if (length(weights) == 1) {
      weights <- rep(weights, length(x))
    } else {
      # Try to adjust weights to match the data
      weights <- weights[1:min(length(weights), length(valid_indices))]
      weights <- weights[valid_indices]

      # If still not matching, handle the error condition
      if (length(weights) != length(x)) {
        warning("Weights length doesn't match vector length after NA removal. Using equal weights.")
        weights <- rep(1, length(x))
      }
    }
  }

  # Ensure weights are positive
  weights[weights <= 0] <- min(weights[weights > 0], na.rm = TRUE)
  if (all(is.na(weights))) weights <- rep(1, length(x))

  # Center the variables
  x_centered <- x - weighted.mean(x, weights, na.rm = TRUE)
  y_centered <- y - weighted.mean(y, weights, na.rm = TRUE)

  # Calculate weighted covariance
  weighted_cov <- sum(weights * x_centered * y_centered, na.rm = TRUE) / sum(weights, na.rm = TRUE)

  # Calculate weighted standard deviations
  sd_x <- sqrt(sum(weights * x_centered^2, na.rm = TRUE) / sum(weights, na.rm = TRUE))
  sd_y <- sqrt(sum(weights * y_centered^2, na.rm = TRUE) / sum(weights, na.rm = TRUE))

  # Check for division by zero
  if (sd_x <= .Machine$double.eps || sd_y <= .Machine$double.eps) {
    return(NA)
  }

  # Return weighted correlation
  return(weighted_cov / (sd_x * sd_y))
}

# Function to benchmark the correlation approach
run_benchmarking <- function(cor_matrix, trajectories, true_signal, true_arrays) {
  n_trajectories <- nrow(trajectories)

  # Calculate correlation of each trajectory with true signal
  true_cors <- sapply(1:n_trajectories, function(i) {
    cor(trajectories[i,], true_signal, method = "pearson")
  })

  # Calculate metrics
  # 1. Accuracy of identifying true arrays
  # Sort trajectories by average correlation with other trajectories
  avg_cors <- rowMeans(cor_matrix, na.rm = TRUE)
  ranked_trajectories <- order(avg_cors, decreasing = TRUE)

  precision <- sum(ranked_trajectories[1:length(true_arrays)] %in% true_arrays) / length(true_arrays)
  recall <- sum(ranked_trajectories[1:length(true_arrays)] %in% true_arrays) / length(true_arrays)
  f1_score <- 2 * (precision * recall) / (precision + recall)

  # 2. Correlation between true/noise classification and our ranking
  true_false_vector <- rep(0, n_trajectories)
  true_false_vector[true_arrays] <- 1

  cor_with_truth <- cor(avg_cors, true_false_vector, method = "spearman")

  # 3. Mean correlation among true arrays vs. overall
  true_array_pairs <- expand.grid(true_arrays, true_arrays)
  true_array_pairs <- true_array_pairs[true_array_pairs$Var1 != true_array_pairs$Var2,]

  mean_cor_true <- mean(sapply(1:nrow(true_array_pairs), function(i) {
    cor_matrix[true_array_pairs$Var1[i], true_array_pairs$Var2[i]]
  }))

  mean_cor_overall <- mean(cor_matrix[upper.tri(cor_matrix)])

  # Return benchmark results
  return(list(
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    correlation_with_truth = cor_with_truth,
    mean_correlation_true_arrays = mean_cor_true,
    mean_correlation_overall = mean_cor_overall,
    true_signal_correlations = true_cors
  ))
}

# benchmark_algorithms <- function(
#     n_simulations = 3,
#     n_true_arrays = 20,
#     n_noise_arrays = 5,
#     n_positions = 100,
#     base_expression = 3,
#     signal_strength = 2,
#     signal_noise_sd_values = c(0.3, 0.5, 1.0),
#     noise_sd_values = c(0.5, 1.0, 2.0),
#     correlation_methods = c("pearson", "spearman", "kendall", "weighted"),
#     smooth_options = c(TRUE, FALSE)
# ) {
#   # Initialize results storage
#   all_results <- list()
#
#   # Run simulations with different noise levels
#   simulation_counter <- 1
#
#   for (signal_noise in signal_noise_sd_values) {
#     for (noise_sd in noise_sd_values) {
#       for (sim in 1:n_simulations) {
#         cat(sprintf("Running simulation %d/%d with signal noise %.2f and noise SD %.2f\n",
#                    sim, n_simulations, signal_noise, noise_sd))
#
#         # Generate simulation data
#         sim_data <- generate_simulation_data(
#           n_true_arrays = n_true_arrays,
#           n_noise_arrays = n_noise_arrays,
#           n_positions = n_positions,
#           base_expression = base_expression,
#           signal_strength = signal_strength,
#           signal_noise_sd = signal_noise,
#           noise_sd = noise_sd,
#           seed = simulation_counter
#         )
#
#         simulation_counter <- simulation_counter + 1
#
#         #-------------------------------------------------------
#         # Method 1: Bayesian GAM Regression Model
#         #-------------------------------------------------------
#
#         # Prepare data for Bayesian GAM
#         input_data <- formDataset(sim_data$data, sim_data$x)
#
#         # Time the fitting process
#         bayes_start_time <- Sys.time()
#
#         # Fit Bayesian model - handle potential errors
#         bayes_fit <- tryCatch({
#           bayesian_gam_regression_nb_shape(
#             x = input_data[[2]],
#             y = input_data[[1]],
#             array_idx = input_data[[3]],
#             n_samples = 1000  # Reduced for speed in benchmarking
#           )
#         }, error = function(e) {
#           message("Error in Bayesian fit: ", e$message)
#           # Return a placeholder with default weights
#           return(list(
#             array_weights = data.frame(
#               Estimate = rep(0, n_true_arrays + n_noise_arrays),
#               weight_norm = rep(1, n_true_arrays + n_noise_arrays)
#             )
#           ))
#         })
#
#         bayes_end_time <- Sys.time()
#         bayes_runtime <- as.numeric(difftime(bayes_end_time, bayes_start_time, units = "secs"))
#
#         # Extract array weights - ensure they exist and have proper form
#         if (!is.null(bayes_fit$array_weights)) {
#           bayes_weights <- bayes_fit$array_weights
#           if (!"weight_norm" %in% names(bayes_weights)) {
#             bayes_weights$weight_norm <- rep(1, n_true_arrays + n_noise_arrays)
#           }
#         } else {
#           # Default weights if we couldn't get them
#           bayes_weights <- data.frame(
#             Estimate = rep(0, n_true_arrays + n_noise_arrays),
#             weight_norm = rep(1, n_true_arrays + n_noise_arrays)
#           )
#         }
#
#         # Calculate accuracy of Bayesian model in identifying true arrays
#         if ("Estimate" %in% names(bayes_weights)) {
#           weights_ordered <- order(bayes_weights$Estimate, decreasing = TRUE)
#           bayes_precision <- sum(weights_ordered[1:n_true_arrays] %in% sim_data$true_arrays) / n_true_arrays
#
#           # Create true/false vector
#           true_false_vector <- rep(0, n_true_arrays + n_noise_arrays)
#           true_false_vector[sim_data$true_arrays] <- 1
#
#           # Correlation between weights and true signal classification
#           bayes_correlation <- cor(bayes_weights$Estimate, true_false_vector, method = "spearman")
#         } else {
#           bayes_precision <- 0
#           bayes_correlation <- 0
#         }
#
#         #-------------------------------------------------------
#         # Method 2: Correlation-based Methods
#         #-------------------------------------------------------
#
#         # For each correlation method and smooth option
#         for (method in correlation_methods) {
#           for (smooth in smooth_options) {
#             # Time the correlation calculation
#             cor_start_time <- Sys.time()
#
#             # Calculate correlations - handle potential errors
#             cor_results <- tryCatch({
#               calculate_trajectory_correlations(
#                 trajectories = sim_data$data,
#                 method = method,
#                 weights = NULL,  # No longer using array weights since the function now creates its own
#                 smooth_first = smooth,
#                 benchmark = TRUE,
#                 true_signal = sim_data$true_signal,
#                 true_arrays = sim_data$true_arrays
#               )
#             }, error = function(e) {
#               message("Error in calculating correlations: ", e$message)
#               # Return empty results
#               return(list(
#                 correlation_matrix = matrix(NA, n_true_arrays + n_noise_arrays, n_true_arrays + n_noise_arrays),
#                 benchmark = list(
#                   precision = 0,
#                   recall = 0,
#                   f1_score = 0,
#                   correlation_with_truth = 0,
#                   mean_correlation_true_arrays = 0,
#                   mean_correlation_overall = 0
#                 )
#               ))
#             })
#
#             cor_end_time <- Sys.time()
#             cor_runtime <- as.numeric(difftime(cor_end_time, cor_start_time, units = "secs"))
#
#             # Store results
#             result_name <- sprintf("sim_%d_signal_noise_%.2f_noise_sd_%.2f_method_%s_smooth_%s",
#                                   sim, signal_noise, noise_sd, method, smooth)
#
#             all_results[[result_name]] <- list(
#               params = list(
#                 sim = sim,
#                 signal_noise_sd = signal_noise,
#                 noise_sd = noise_sd,
#                 method = method,
#                 smooth = smooth
#               ),
#               bayes = list(
#                 precision = bayes_precision,
#                 correlation = bayes_correlation,
#                 runtime = bayes_runtime,
#                 weights = bayes_weights
#               ),
#               correlation = list(
#                 precision = cor_results$benchmark$precision,
#                 recall = cor_results$benchmark$recall,
#                 f1_score = cor_results$benchmark$f1_score,
#                 correlation_with_truth = cor_results$benchmark$correlation_with_truth,
#                 mean_correlation_true_arrays = cor_results$benchmark$mean_correlation_true_arrays,
#                 mean_correlation_overall = cor_results$benchmark$mean_correlation_overall,
#                 runtime = cor_runtime
#               ),
#               correlation_matrix = cor_results$correlation_matrix
#             )
#           }
#         }
#       }
#     }
#   }
#
#   # Create summary data frame
#   summary_df <- do.call(rbind, lapply(names(all_results), function(name) {
#     result <- all_results[[name]]
#     data.frame(
#       simulation = result$params$sim,
#       signal_noise_sd = result$params$signal_noise_sd,
#       noise_sd = result$params$noise_sd,
#       method = result$params$method,
#       smooth = result$params$smooth,
#       bayes_precision = result$bayes$precision,
#       bayes_correlation = result$bayes$correlation,
#       bayes_runtime = result$bayes$runtime,
#       cor_precision = result$correlation$precision,
#       cor_recall = result$correlation$recall,
#       cor_f1_score = result$correlation$f1_score,
#       cor_correlation = result$correlation$correlation_with_truth,
#       cor_runtime = result$correlation$runtime
#     )
#   }))
#
#   return(list(
#     detailed_results = all_results,
#     summary = summary_df
#   ))
# }

generate_simulation_data <- function(
    n_true_arrays = 20,
    n_noise_arrays = 5,
    n_positions = 100,
    base_expression = 3,
    signal_strength = 2,
    signal_noise_sd = 0.5,
    overdispersion =4 ,  # different overdispersion levels
    noise_sd = 1,  # standard deviation for Gaussian noise
    seed = 123
) {
  set.seed(seed)

  n_total_arrays <- n_true_arrays + n_noise_arrays

  # Generate position vector
  x <- seq(0, 1, length.out = n_positions)

  # Create polynomial signal pattern
  true_signal <- base_expression +
    signal_strength * (-x^3 + 0.5*x^2 + 0.3*x)  # cubic polynomial

  # Ensure signal is positive
  true_signal <- true_signal - min(true_signal) + 1

  # Create empty matrix
  sim_data <- matrix(
    NA,
    nrow = n_total_arrays,
    ncol = n_positions,
    dimnames = list(
      paste0("Array", 1:n_total_arrays),
      paste0("Pos", 1:n_positions)
    )
  )

  # Generate true signal arrays (1:20)
  for(i in 1:n_true_arrays) {
    # Add random noise to the true signal
    noisy_signal <- true_signal + rnorm(n_positions, mean = 0, sd = signal_noise_sd)
    # Ensure positivity
    noisy_signal <- pmax(noisy_signal, 0)
    sim_data[i, ] <- noisy_signal
  }

  # Generate pure noise arrays (21:25)
  for(i in (n_true_arrays + 1):n_total_arrays) {
    # Pure random Gaussian noise around base expression
    sim_data[i, ] <- pmax(rnorm(n_positions,
                           mean = base_expression,
                           sd = noise_sd),0)
  }

  return(list(
    data = sim_data,
    x = x,
    true_signal = true_signal,
    true_arrays = 1:n_true_arrays,
    noise_arrays = (n_true_arrays + 1):n_total_arrays
  ))
}

# Function for Dynamic Time Warping similarity
calculate_dtw_similarity <- function(trajectories, smooth_first = FALSE) {
  # Check if dtw package is available
  if (!requireNamespace("dtw", quietly = TRUE)) {
    stop("Package 'dtw' is required for DTW similarity. Please install it.")
  }

  n_trajectories <- nrow(trajectories)

  # Smooth trajectories if requested (should be handled by parent function)
  if (smooth_first) {
    warning("Smoothing should be handled by the parent function. Ignoring here.")
  }

  # Initialize similarity matrix
  sim_matrix <- matrix(0,
                      nrow = n_trajectories,
                      ncol = n_trajectories,
                      dimnames = list(rownames(trajectories), rownames(trajectories)))

  # Track the maximum distance for normalization
  max_distance <- 0

  # Calculate DTW distances
  for (i in 1:n_trajectories) {
    for (j in i:n_trajectories) {
      if (i == j) {
        # Distance to self is 0, similarity is 1
        sim_matrix[i, j] <- 1
      } else {
        # Calculate DTW distance
        tryCatch({
          alignment <- dtw::dtw(trajectories[i,], trajectories[j,],
                               window.type = "sakoechiba", window.size = 10)

          # Get normalized distance
          distance <- alignment$normalizedDistance

          # Update max distance if needed
          if (distance > max_distance && is.finite(distance)) {
            max_distance <- distance
          }

          # Store distance temporarily
          sim_matrix[i, j] <- sim_matrix[j, i] <- distance
        }, error = function(e) {
          warning(paste("DTW calculation failed for arrays", i, "and", j, ":", e$message))
          sim_matrix[i, j] <- sim_matrix[j, i] <- NA
        })
      }
    }
  }

  # Convert non-diagonal entries from distances to similarities
  for (i in 1:n_trajectories) {
    for (j in 1:n_trajectories) {
      if (i != j) {
        if (is.na(sim_matrix[i, j]) || max_distance == 0) {
          sim_matrix[i, j] <- 0
        } else {
          # Convert distance to similarity (1 - normalized distance)
          sim_matrix[i, j] <- 1 - (sim_matrix[i, j] / max_distance)
        }
      }
    }
  }

  return(sim_matrix)
}


#---------------------------------------------------------
# Benchmarking script to compare algorithms
#---------------------------------------------------------

benchmark_algorithms <- function(
    n_simulations = 5,
    n_true_arrays = 20,
    n_noise_arrays = 5,
    n_positions = 50,
    base_expression = 3,
    signal_strength = 5,
    signal_noise_sd_values = c( 0.5, 1.0,2,4,6),
    noise_sd_values = c(1.0),
    correlation_methods = c("pearson", "spearman", "weighted", "dtw"),  # Add DTW as an option
    smooth_options = c(TRUE, FALSE)
) {
  # Initialize results storage
  all_results <- list()

  # Run simulations with different noise levels
  simulation_counter <- 1

  for (signal_noise in signal_noise_sd_values) {
    for (noise_sd in noise_sd_values) {
      for (sim in 1:n_simulations) {
        cat(sprintf("Running simulation %d/%d with signal noise %.2f and noise SD %.2f\n",
                    sim, n_simulations, signal_noise, noise_sd))

        # Generate simulation data
        sim_data <- generate_simulation_data(
          n_true_arrays = n_true_arrays,
          n_noise_arrays = n_noise_arrays,
          n_positions = n_positions,
          base_expression = base_expression,
          signal_strength = signal_strength,
          signal_noise_sd = signal_noise,
          noise_sd = noise_sd,
          seed = simulation_counter
        )

        simulation_counter <- simulation_counter + 1

        #-------------------------------------------------------
        # Method 1: Bayesian GAM Regression Model
        #-------------------------------------------------------

        # Use the Bayesian model function to fit and evaluate
        bayes_results <- fit_bayesian_model(sim_data, n_samples = 2000)

        # Extract the necessary results
        bayes_fit <- bayes_results$fit
        bayes_weights <- bayes_results$weights
        bayes_precision <- bayes_results$precision
        bayes_correlation <- bayes_results$correlation
        bayes_runtime <- bayes_results$runtime

        #-------------------------------------------------------
        # Method 2: Correlation-based Methods
        #-------------------------------------------------------

        # For each correlation method and smooth option
        for (method in correlation_methods) {
          for (smooth in smooth_options) {
            # Time the correlation calculation
            cor_start_time <- Sys.time()

            # Calculate correlations - handle potential errors
            cor_results <- tryCatch({
              calculate_trajectory_correlations(
                trajectories = sim_data$data,
                method = method,
                weights = NULL,  # No longer using array weights since the function now creates its own
                smooth_first = smooth,
                benchmark = TRUE,
                true_signal = sim_data$true_signal,
                true_arrays = sim_data$true_arrays
              )
            }, error = function(e) {
              message("Error in calculating correlations: ", e$message)
              # Return empty results
              return(list(
                correlation_matrix = matrix(NA, n_true_arrays + n_noise_arrays, n_true_arrays + n_noise_arrays),
                benchmark = list(
                  precision = 0,
                  recall = 0,
                  f1_score = 0,
                  correlation_with_truth = 0,
                  mean_correlation_true_arrays = 0,
                  mean_correlation_overall = 0
                )
              ))
            })

            cor_end_time <- Sys.time()
            cor_runtime <- as.numeric(difftime(cor_end_time, cor_start_time, units = "secs"))

            # Store results
            result_name <- sprintf("sim_%d_signal_noise_%.2f_noise_sd_%.2f_method_%s_smooth_%s",
                                   sim, signal_noise, noise_sd, method, smooth)

            all_results[[result_name]] <- list(
              params = list(
                sim = sim,
                signal_noise_sd = signal_noise,
                noise_sd = noise_sd,
                method = method,
                smooth = smooth
              ),
              bayes = list(
                precision = bayes_precision,
                correlation = bayes_correlation,
                runtime = bayes_runtime,
                weights = bayes_weights
              ),
              correlation = list(
                precision = cor_results$benchmark$precision,
                recall = cor_results$benchmark$recall,
                f1_score = cor_results$benchmark$f1_score,
                correlation_with_truth = cor_results$benchmark$correlation_with_truth,
                mean_correlation_true_arrays = cor_results$benchmark$mean_correlation_true_arrays,
                mean_correlation_overall = cor_results$benchmark$mean_correlation_overall,
                runtime = cor_runtime
              ),
              correlation_matrix = cor_results$correlation_matrix
            )
          }
        }
      }
    }
  }

  # Create summary data frame - add error handling
  summary_df <- tryCatch({
    do.call(rbind, lapply(names(all_results), function(name) {
      result <- all_results[[name]]
      data.frame(
        simulation = result$params$sim,
        signal_noise_sd = result$params$signal_noise_sd,
        noise_sd = result$params$noise_sd,
        method = result$params$method,
        smooth = result$params$smooth,
        bayes_precision = result$bayes$precision,
        bayes_correlation = result$bayes$correlation,
        bayes_runtime = result$bayes$runtime,
        cor_precision = result$correlation$precision,
        cor_recall = result$correlation$recall,
        cor_f1_score = result$correlation$f1_score,
        cor_correlation = result$correlation$correlation_with_truth,
        cor_runtime = result$correlation$runtime
      )
    }))
  }, error = function(e) {
    message("Error creating summary data frame: ", e$message)
    return(data.frame())
  })

  return(list(
    detailed_results = all_results,
    summary = summary_df
  ))
}

#---------------------------------------------------------
# Function to fit and evaluate Bayesian model
#---------------------------------------------------------

fit_bayesian_model <- function(sim_data, n_samples = 2000) {
  # Get array dimensions
  n_true_arrays <- length(sim_data$true_arrays)
  n_noise_arrays <- length(sim_data$noise_arrays)
  n_total_arrays <- n_true_arrays + n_noise_arrays

  # Prepare data for Bayesian GAM
  input_data <- formDataset(sim_data$data, sim_data$x)

  # Time the fitting process
  bayes_start_time <- Sys.time()

  # Fit Bayesian model - handle potential errors
  bayes_fit <- tryCatch({
    bayesian_gam_regression_nb_shape(
      x = input_data[[2]],
      y = input_data[[1]],
      array_idx = input_data[[3]],
      n_samples = n_samples  # Reduced for speed in benchmarking
    )
  }, error = function(e) {
    message("Error in Bayesian fit: ", e$message)
    # Return a placeholder with default weights
    return(list(
      array_weights = data.frame(
        Estimate = rep(0, n_total_arrays),
        weight_norm = rep(1, n_total_arrays)
      )
    ))
  })

  bayes_end_time <- Sys.time()
  bayes_runtime <- as.numeric(difftime(bayes_end_time, bayes_start_time, units = "secs"))

  # Extract array weights - ensure they exist and have proper form
  if (!is.null(bayes_fit$array_weights)) {
    bayes_weights <- bayes_fit$array_weights
    if (!"weight_norm" %in% names(bayes_weights)) {
      bayes_weights$weight_norm <- rep(1, n_total_arrays)
    }
  } else {
    # Default weights if we couldn't get them
    bayes_weights <- data.frame(
      Estimate = rep(0, n_total_arrays),
      weight_norm = rep(1, n_total_arrays)
    )
  }

  # Calculate accuracy of Bayesian model in identifying true arrays
  if ("Estimate" %in% names(bayes_weights)) {
    weights_ordered <- order(bayes_weights$Estimate, decreasing = TRUE)
    bayes_precision <- sum(weights_ordered[1:n_true_arrays] %in% sim_data$true_arrays) / n_true_arrays

    # Create true/false vector
    true_false_vector <- rep(0, n_total_arrays)
    true_false_vector[sim_data$true_arrays] <- 1

    # Correlation between weights and true signal classification
    bayes_correlation <- cor(bayes_weights$Estimate, true_false_vector, method = "spearman")
  } else {
    bayes_precision <- 0
    bayes_correlation <- 0
  }

  # Return results
  return(list(
    fit = bayes_fit,
    weights = bayes_weights,
    precision = bayes_precision,
    correlation = bayes_correlation,
    runtime = bayes_runtime
  ))
}


#---------------------------------------------------------
# Visualization functions for benchmark data
#---------------------------------------------------------

visualize_simulation <- function(sim_data, results = NULL, output_file = NULL) {
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for visualization. Please install it.")
  }

  # Load required visualization packages
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  require(patchwork)
  require(viridis)  # For better color scales

  # Extract data from sim_data
  trajectories <- sim_data$data
  x <- sim_data$x
  true_signal <- sim_data$true_signal
  true_arrays <- sim_data$true_arrays
  noise_arrays <- sim_data$noise_arrays

  n_true <- length(true_arrays)
  n_noise <- length(noise_arrays)

  # Prepare data for plotting

  # 1. True signal for reference
  true_signal_df <- data.frame(
    x = x,
    y = true_signal,
    type = "True Signal"
  )

  # 2. Sample of true arrays (max 3 for clarity)
  sample_true <- sample(true_arrays, min(3, n_true))
  true_arrays_df <- do.call(rbind, lapply(sample_true, function(i) {
    data.frame(
      x = x,
      y = trajectories[i,],
      type = paste0("True Array ", i)
    )
  }))

  # 3. Sample of noise arrays (max 3 for clarity)
  sample_noise <- sample(noise_arrays, min(3, n_noise))
  noise_arrays_df <- do.call(rbind, lapply(sample_noise, function(i) {
    data.frame(
      x = x,
      y = trajectories[i,],
      type = paste0("Noise Array ", i)
    )
  }))

  # Combine all trajectory data
  trajectory_df <- rbind(
    true_signal_df,
    true_arrays_df,
    noise_arrays_df
  )

  # Create trajectory plot
  p_trajectories <- ggplot(trajectory_df, aes(x = x, y = y, color = type, linetype = type)) +
    geom_line(size = 1) +
    scale_color_manual(
      values = c("True Signal" = "black",
                 setNames(viridis::viridis(length(sample_true)), paste0("True Array ", sample_true)),
                 setNames(RColorBrewer::brewer.pal(min(length(sample_noise), 8), "Set2"), paste0("Noise Array ", sample_noise))),
      name = "Trajectory Type"
    ) +
    scale_linetype_manual(
      values = c("True Signal" = 1,
                 rep(1, length(sample_true)),
                 rep(2, length(sample_noise))),
      name = "Trajectory Type"
    ) +
    labs(title = "Simulation Trajectories",
         subtitle = paste0(n_true, " true arrays and ", n_noise, " noise arrays"),
         x = "Position", y = "Expression") +
    theme_minimal()

  # Calculate correlation matrix if not provided
  if (is.null(results)) {
    # Calculate Pearson correlations
    pearson_matrix <- calculate_trajectory_correlations(
      trajectories = trajectories,
      method = "pearson",
      smooth_first = FALSE
    )

    # Try to calculate DTW similarities if dtw package is available
    dtw_matrix <- tryCatch({
      calculate_trajectory_correlations(
        trajectories = trajectories,
        method = "dtw",
        smooth_first = FALSE
      )
    }, error = function(e) {
      message("DTW calculation failed, using only Pearson correlations.")
      NULL
    })
  } else {
    # Extract correlation matrices from results
    pearson_matrix <- results$correlation_matrix
    dtw_matrix <- NULL  # Can be extended if DTW results are provided
  }

  # Create heatmap data for Pearson correlations
  heatmap_data <- expand.grid(
    Array1 = 1:nrow(trajectories),
    Array2 = 1:nrow(trajectories)
  )
  heatmap_data$Correlation <- NA
  for (i in 1:nrow(heatmap_data)) {
    row <- heatmap_data$Array1[i]
    col <- heatmap_data$Array2[i]
    heatmap_data$Correlation[i] <- pearson_matrix[row, col]
  }

  # Add array type information
  heatmap_data$Type1 <- ifelse(heatmap_data$Array1 %in% true_arrays, "True", "Noise")
  heatmap_data$Type2 <- ifelse(heatmap_data$Array2 %in% true_arrays, "True", "Noise")
  heatmap_data$PairType <- paste(heatmap_data$Type1, "-", heatmap_data$Type2)

  # Create correlation heatmap
  p_heatmap <- ggplot(heatmap_data, aes(x = Array1, y = Array2, fill = Correlation)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(-1, 1), option = "plasma") +
    labs(title = "Pearson Correlation Matrix",
         x = "Array", y = "Array",
         fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8)) +
    # Highlight true/noise arrays
    geom_vline(xintercept = c(min(true_arrays) - 0.5, max(true_arrays) + 0.5),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = c(min(true_arrays) - 0.5, max(true_arrays) + 0.5),
               linetype = "dashed", color = "black")

  # Create boxplot of correlation distributions by array pair type
  p_boxplot <- ggplot(subset(heatmap_data, Array1 != Array2), aes(x = PairType, y = Correlation, fill = PairType)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = "Distribution of Correlations by Array Type",
         x = "Pair Type", y = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # If array weights are provided in results, visualize them
  p_weights <- NULL
  if (!is.null(results) && !is.null(results$bayes) && !is.null(results$bayes$weights)) {
    weights_df <- data.frame(
      Array = 1:nrow(trajectories),
      Weight = results$bayes$weights$Estimate,
      Type = ifelse(1:nrow(trajectories) %in% true_arrays, "True", "Noise")
    )

    p_weights <- ggplot(weights_df, aes(x = Array, y = Weight, fill = Type)) +
      geom_col() +
      scale_fill_manual(values = c("True" = "darkgreen", "Noise" = "darkred")) +
      labs(title = "Bayesian Model Array Weights",
           subtitle = paste("Precision:", round(results$bayes$precision, 2)),
           x = "Array", y = "Weight") +
      theme_minimal()
  }

  # Combine plots
  if (!is.null(p_weights)) {
    combined_plot <- (p_trajectories / p_weights) | (p_heatmap / p_boxplot)
  } else {
    combined_plot <- p_trajectories | (p_heatmap / p_boxplot)
  }

  # Save plot if output file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, combined_plot, width = 15, height = 10)
    message(paste("Visualization saved to", output_file))
  }

  # Return plots for further customization
  invisible(list(
    trajectory_plot = p_trajectories,
    correlation_heatmap = p_heatmap,
    correlation_boxplot = p_boxplot,
    weights_plot = p_weights,
    combined_plot = combined_plot
  ))
}

# Function to visualize benchmark results
visualize_benchmark <- function(benchmark_results, output_dir = "results/trajectory") {
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for visualization. Please install it.")
  }

  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Load required packages
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(patchwork)

  # Extract summary data
  summary_df <- benchmark_results$summary

  # Check if we have data to visualize
  if (nrow(summary_df) == 0) {
    message("No benchmark data to visualize.")
    return(invisible(NULL))
  }

  # Summarize data by groups
  summary_plot_data <- summary_df %>%
    group_by(method, smooth, signal_noise_sd, noise_sd) %>%
    summarize(
      bayes_precision_mean = mean(bayes_precision, na.rm = TRUE),
      bayes_correlation_mean = mean(bayes_correlation, na.rm = TRUE),
      bayes_runtime_mean = mean(bayes_runtime, na.rm = TRUE),
      cor_precision_mean = mean(cor_precision, na.rm = TRUE),
      cor_recall_mean = mean(cor_recall, na.rm = TRUE),
      cor_f1_score_mean = mean(cor_f1_score, na.rm = TRUE),
      cor_correlation_mean = mean(cor_correlation, na.rm = TRUE),
      cor_runtime_mean = mean(cor_runtime, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # 1. Create method comparison plot by noise level
  p_method <- ggplot(summary_plot_data, aes(x = signal_noise_sd, y = cor_f1_score_mean,
                                          color = method, shape = factor(smooth))) +
    geom_point(size = 3) +
    geom_line(aes(linetype = factor(smooth))) +
    facet_wrap(~noise_sd, labeller = labeller(noise_sd = function(x) paste("Noise SD:", x))) +
    labs(title = "Method Performance by Noise Level",
         x = "Signal Noise SD",
         y = "F1 Score",
         color = "Method",
         shape = "Smoothing",
         linetype = "Smoothing") +
    theme_minimal() +
    scale_color_viridis_d(option = "plasma", end = 0.9)

  # 2. Bayes vs Correlation methods comparison
  p_comparison <- summary_plot_data %>%
    mutate(method_smooth = paste(method, ifelse(smooth, "Smoothed", "Raw"))) %>%
    ggplot(aes(x = bayes_precision_mean, y = cor_precision_mean,
              color = method, shape = factor(smooth))) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    facet_grid(signal_noise_sd ~ noise_sd,
              labeller = labeller(signal_noise_sd = function(x) paste("Signal Noise:", x),
                                 noise_sd = function(x) paste("Noise SD:", x))) +
    labs(title = "Bayesian vs Correlation Methods",
         subtitle = "Points above line indicate correlation method outperforms Bayesian model",
         x = "Bayesian Model Precision",
         y = "Correlation Method Precision",
         color = "Method",
         shape = "Smoothing") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_viridis_d(option = "plasma", end = 0.9)

  # 3. Runtime comparison
  p_runtime <- summary_plot_data %>%
    mutate(method_smooth = paste(method, ifelse(smooth, "Smoothed", "Raw"))) %>%
    ggplot(aes(x = method_smooth, y = cor_runtime_mean, fill = method)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = mean(bayes_runtime_mean)),
              color = "red", linetype = "dashed", size = 1) +
    labs(title = "Runtime Comparison",
         subtitle = "Red line shows mean Bayesian model runtime",
         x = "Method & Smoothing",
         y = "Runtime (seconds)",
         fill = "Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_viridis_d(option = "plasma", end = 0.9)

  # 4. Precision vs Noise level for all methods
  p_noise_impact <- ggplot() +
    # Correlation method lines
    geom_line(data = summary_plot_data,
             aes(x = signal_noise_sd, y = cor_precision_mean,
                 color = method, linetype = factor(smooth)),
             size = 1) +
    # Bayesian model reference line
    geom_line(data = summary_plot_data %>%
                 group_by(signal_noise_sd, noise_sd) %>%
                 summarize(bayes_mean = mean(bayes_precision_mean), .groups = "drop"),
               aes(x = signal_noise_sd, y = bayes_mean),
               color = "black", size = 1.5, alpha = 0.7) +
    facet_wrap(~noise_sd, labeller = labeller(noise_sd = function(x) paste("Noise SD:", x))) +
    labs(title = "Method Performance vs Signal Noise",
         subtitle = "Black line shows Bayesian model performance",
         x = "Signal Noise SD",
         y = "Precision",
         color = "Method",
         linetype = "Smoothing") +
    theme_minimal() +
    scale_color_viridis_d(option = "plasma", end = 0.9)

  # Combine plots into a grid
  p_grid1 <- p_method / p_noise_impact
  p_grid2 <- p_comparison / p_runtime
  combined_plot <- p_grid1 | p_grid2

  # Save individual plots
  ggsave(paste0(output_dir, "/benchmark_method_comparison.pdf"), p_method, width = 10, height = 8)
  ggsave(paste0(output_dir, "/benchmark_bayes_vs_correlation.pdf"), p_comparison, width = 12, height = 10)
  ggsave(paste0(output_dir, "/benchmark_runtime.pdf"), p_runtime, width = 8, height = 6)
  ggsave(paste0(output_dir, "/benchmark_noise_impact.pdf"), p_noise_impact, width = 10, height = 8)

  # Save combined plot
  ggsave(paste0(output_dir, "/benchmark_comprehensive.pdf"), combined_plot, width = 18, height = 12)

  message("Benchmark visualization complete. Results saved to ", output_dir)

  # Return plots for further customization
  invisible(list(
    method_comparison = p_method,
    bayes_vs_correlation = p_comparison,
    runtime_comparison = p_runtime,
    noise_impact = p_noise_impact,
    combined_plot = combined_plot
  ))
}

# Simple function to create an interactive visualization using plotly
create_interactive_viz <- function(sim_data, output_file = NULL) {
  # Check if plotly is available
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive visualization. Please install it.")
  }

  require(plotly)
  require(tidyr)
  require(dplyr)

  # Extract data from sim_data
  trajectories <- sim_data$data
  x <- sim_data$x
  true_signal <- sim_data$true_signal
  true_arrays <- sim_data$true_arrays
  noise_arrays <- sim_data$noise_arrays

  # Create long format data
  trajectory_data <- as.data.frame(t(trajectories))
  colnames(trajectory_data) <- paste0("Array", 1:ncol(trajectory_data))
  trajectory_data$Position <- x
  trajectory_data$TrueSignal <- true_signal

  # Convert to long format for plotting
  plot_data <- trajectory_data %>%
    pivot_longer(cols = -c(Position, TrueSignal),
                 names_to = "Array",
                 values_to = "Expression") %>%
    mutate(
      ArrayNum = as.numeric(gsub("Array", "", Array)),
      Type = ifelse(ArrayNum %in% true_arrays, "True Array", "Noise Array")
    )

  # Create plotly figure
  fig <- plot_ly()

  # Add true signal
  fig <- fig %>% add_trace(
    data = trajectory_data,
    x = ~Position,
    y = ~TrueSignal,
    type = 'scatter',
    mode = 'lines',
    name = 'True Signal',
    line = list(color = 'black', width = 3)
  )

  # Add true arrays
  for (i in true_arrays) {
    array_data <- subset(plot_data, ArrayNum == i)
    fig <- fig %>% add_trace(
      x = array_data$Position,
      y = array_data$Expression,
      type = 'scatter',
      mode = 'lines',
      name = paste0('True Array ', i),
      line = list(dash = 'solid'),
      visible = ifelse(i <= 3, TRUE, 'legendonly')  # Show only first 3 by default
    )
  }

  # Add noise arrays
  for (i in noise_arrays) {
    array_data <- subset(plot_data, ArrayNum == i)
    fig <- fig %>% add_trace(
      x = array_data$Position,
      y = array_data$Expression,
      type = 'scatter',
      mode = 'lines',
      name = paste0('Noise Array ', i),
      line = list(dash = 'dot'),
      visible = ifelse(i <= min(noise_arrays) + 2, TRUE, 'legendonly')  # Show only first 3 by default
    )
  }

  # Update layout
  fig <- fig %>% layout(
    title = "Interactive Simulation Trajectories",
    xaxis = list(title = "Position"),
    yaxis = list(title = "Expression"),
    legend = list(title = list(text = "Arrays")),
    hovermode = "closest"
  )

  # Save as HTML if output file specified
  if (!is.null(output_file)) {
    htmlwidgets::saveWidget(fig, output_file)
    message(paste("Interactive visualization saved to", output_file))
  }

  return(fig)
}

simple_visualize_signal <- function(sim_data, output_file = NULL, show_plot = TRUE) {
  # Make sure we have ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for visualization. Please install it.")
  }

  # Load basic visualization packages
  library(ggplot2)

  # Extract the data
  x <- sim_data$x  # positions vector
  true_signal <- sim_data$true_signal  # The true underlying signal
  trajectories <- sim_data$data  # All trajectories matrix
  true_arrays <- sim_data$true_arrays  # Which arrays contain signal
  noise_arrays <- sim_data$noise_arrays  # Which arrays are just noise

  # Create a data frame for the true signal
  df_true_signal <- data.frame(x = x, y = true_signal, type = "True Signal", category = "True Signal")

  # Create a data frame for true arrays
  df_true_arrays <- data.frame()
  for (i in true_arrays[1:min(3, length(true_arrays))]) {
    temp_df <- data.frame(
      x = x,
      y = trajectories[i, ],
      type = paste0("True Array ", i),
      category = "True Arrays"  # All true arrays share the same category
    )
    df_true_arrays <- rbind(df_true_arrays, temp_df)
  }

  # Create a data frame for noise arrays
  df_noise_arrays <- data.frame()
  for (i in noise_arrays[1:min(3, length(noise_arrays))]) {
    temp_df <- data.frame(
      x = x,
      y = trajectories[i, ],
      type = paste0("Noise Array ", i),
      category = "Noise Arrays"  # All noise arrays share the same category
    )
    df_noise_arrays <- rbind(df_noise_arrays, temp_df)
  }

  # Combine all data
  plot_data <- rbind(df_true_signal, df_true_arrays, df_noise_arrays)

  # Create the plot with separate facets for the three categories
  p <- ggplot(plot_data, aes(x = x, y = y, color = type, linetype = type)) +
    geom_point(size = 2) +  # Adjust point size for visibility
    geom_smooth(se = FALSE, method = "loess") +  # Use loess for smoothing
    labs(
      title = "Trajectory Visualization",
      x = "Position",
      y = "Expression",
      color = "Trajectory Type",
      linetype = "Trajectory Type"
    ) +
    scale_color_manual(values = c("True Signal" = "black",
                                   "True Array 1" = "blue",
                                   "True Array 2" = "royalblue",
                                   "True Array 3" = "steelblue",
                                   "Noise Array 1" = "red",
                                   "Noise Array 2" = "darkred",
                                   "Noise Array 3" = "firebrick")) +  # Customize colors
    scale_linetype_manual(values = c("True Signal" = "solid",
                                     "True Array 1" = "solid",
                                     "True Array 2" = "solid",
                                     "True Array 3" = "solid",
                                     "Noise Array 1" = "dashed",
                                     "Noise Array 2" = "dashed",
                                     "Noise Array 3" = "dashed")) +  # Customize linetypes
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA)
    ) +
    facet_wrap(~ category, scales = "free_y", ncol = 1)  # Create only three facets

  # Save the plot if output file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 10)  # Taller to accommodate three facets
    message(paste("Simple visualization saved to", output_file))
  }

  # Show the plot if requested
  if (show_plot) {
    print(p)
  }

  # Return the plot and data for further customization
  return(list(
    plot = p,
    data = plot_data
  ))
}
