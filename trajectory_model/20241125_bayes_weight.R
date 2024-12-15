bayesian_poly_regression <- function(x, y, array_idx, degree = 3, n_samples = 2000) {
  # Remove NA values and keep track of which array each point belongs to
  valid_idx <- !is.na(y)
  x_clean <- x[valid_idx]
  y_clean <- y[valid_idx]
  array_clean <- array_idx[valid_idx]

  # Normalize x to [0,1]
  x_norm <- (x_clean - min(x_clean)) / (max(x_clean) - min(x_clean))

  # Normalize y for better fitting
  y_mean <- mean(y_clean)
  y_sd <- sd(y_clean)
  y_norm <- (y_clean - y_mean) / y_sd

  # Create polynomial features
  poly_features <- matrix(1, nrow = length(x_norm), ncol = degree + 1)
  for (i in 1:degree) {
    poly_features[, i + 1] <- x_norm^i
  }

  # Manual scaling of features
  feature_means <- colMeans(poly_features)
  feature_sds <- apply(poly_features, 2, sd)
  feature_sds[feature_sds == 0] <- 1

  poly_features_scaled <- sweep(poly_features, 2, feature_means, "-")
  poly_features_scaled <- sweep(poly_features_scaled, 2, feature_sds, "/")

  # Prepare data for Stan
  stan_data <- list(
    N = length(x_norm),
    K = degree + 1,
    M = length(unique(array_clean)),
    X = poly_features_scaled,
    y = y_norm,  # Using normalized y
    array_idx = array_clean
  )

  # Improved Stan model
  stan_code <- "
  data {
    int<lower=0> N;
    int<lower=0> K;
    int<lower=0> M;
    matrix[N, K] X;
    vector[N] y;
    array[N] int<lower=1, upper=M> array_idx;
  }

  parameters {
    vector[K] beta;
    real<lower=0> sigma;
    vector<lower=0>[M] array_weight;  // Changed prior
    real<lower=0> alpha;  // Shape parameter for array weights
  }

  model {
    // Improved priors
    beta ~ student_t(3, 0, 2);  // More robust prior
    sigma ~ student_t(3, 0, 2);
    alpha ~ gamma(2, 0.1);
    array_weight ~ gamma(alpha, alpha);  // More flexible prior

    // Likelihood with improved weighting
    for (n in 1:N) {
      real mu = dot_product(X[n], beta);
      y[n] ~ student_t(3, mu, sigma / sqrt(array_weight[array_idx[n]]));
    }
  }

  generated quantities {
    vector<lower=0, upper=1>[M] array_probs;
    for (m in 1:M) {
      array_probs[m] = array_weight[m] / max(array_weight);
    }
  }
  "

  # Fit the model with more iterations and adaptation
  fit <- stan(
    model_code = stan_code,
    data = stan_data,
    iter = n_samples,
    chains = 4,
    cores = 4,
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 12
    )
  )

  # Extract results
  array_probs <- summary(fit, pars = "array_probs")$summary[, "mean"]
  coefficients <- summary(fit, pars = "beta")$summary[, "mean"]

  return(list(
    coefficients = coefficients,
    array_probs = array_probs,
    fit = fit,
    feature_means = feature_means,
    feature_sds = feature_sds,
    x_range = range(x_clean),
    y_mean = y_mean,
    y_sd = y_sd
  ))
}

# Improved plotting function
plot_results <- function(x, y, array_idx, results) {
  # Plot array probabilities
  p1 <- ggplot(data.frame(index = 1:length(results$array_probs),
                          prob = results$array_probs),
               aes(x = index, y = prob)) +
    geom_point(size = 2) +
    geom_line() +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    labs(title = "Array Signal Probability",
         x = "Array Index",
         y = "Probability") +
    theme_minimal() +
    ylim(0, 1)

  # Create prediction data
  x_seq <- seq(min(x), max(x), length.out = 100)
  x_norm <- (x_seq - min(x)) / (max(x) - min(x))

  # Create polynomial features for prediction
  poly_pred <- matrix(1, nrow = length(x_norm), ncol = length(results$coefficients))
  for (i in 1:(length(results$coefficients)-1)) {
    poly_pred[, i + 1] <- x_norm^i
  }

  # Scale features and predict
  poly_pred_scaled <- sweep(poly_pred, 2, results$feature_means, "-")
  poly_pred_scaled <- sweep(poly_pred_scaled, 2, results$feature_sds, "/")

  # Transform prediction back to original scale
  y_pred <- (as.vector(poly_pred_scaled %*% results$coefficients) *
               results$y_sd + results$y_mean)

  # Create plot data
  df_plot <- data.frame(
    x = x,
    y = y,
    array = factor(array_idx),
    weight = rep(results$array_probs, each = length(x)/length(results$array_probs))
  )

  df_pred <- data.frame(
    x = x_seq,
    y = y_pred
  )

  p2 <- ggplot() +
    geom_point(data = df_plot, aes(x = x, y = y, color = array, alpha = weight)) +
    geom_line(data = df_pred, aes(x = x, y = y),
              color = "red", size = 1) +
    labs(title = "Data and Fitted Polynomial",
         x = "X",
         y = "Y") +
    theme_minimal() +
    theme(legend.position = "none")

  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

# Fit model with more iterations
results <- bayesian_poly_regression(x, y, array_idx, n_samples = 3000)
pdf("results/trajectory/20241129_conserved_trajectory_model/20241125_bayes_test1.pdf",width = 8,height = 4)
plot_results(x, y, array_idx, results)
dev.off()
#ggsave("results/trajectory/20241129_conserved_trajectory_model/20241125_bayes_test1.pdf",width = 8,height = 4)
# Check specific arrays
par(mfrow=c(2,2))
plot(x[1:100], y[1:100], main="Array 1 (Signal)")
plot(x[2401:2500], y[2401:2500], main="Array 24 (Noise)")
plot(x[2301:2400], y[2301:2400], main="Array 23 (Noise)")
plot(x[101:200], y[101:200], main="Array 2 (Signal)")

