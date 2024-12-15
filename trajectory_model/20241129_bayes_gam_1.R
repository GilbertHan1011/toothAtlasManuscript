library(brms)
bayesian_gam_regression <- function(x, y,n_knots = 5, array_idx, n_samples = 1000) {
  # Create data frame
  df <- data.frame(
    y = y[!is.na(y)],
    x = x[!is.na(y)],
    array = factor(array_idx[!is.na(y)])
  )
  init_fit <- mgcv::gam(
    y ~ s(x, bs = "cr") + factor(array_idx),
    family = gaussian()
  )
  init_values <- coef(init_fit)
  # Define formula with array-specific smooths
  formula <- bf(
    y ~ s(x, bs = "cr", k = 5) + (1|array),
    sigma ~ 0 + array
  )

  # Set correct priors
  prior <- c(
    # Prior for smoothing parameters
    prior(student_t(3, 0, 2), class = "b"),
    # Prior for array-specific sigmas
    prior(student_t(3, 0, 2), class = "b", dpar = "sigma"),
    # Prior for random effects SD
    prior(gamma(2, 0.1), class = "sd")
  )

  # Fit the model
  fit <- brm(
    formula = formula,
    data = df,
    family = gaussian,
    prior = prior,
    chains = 4,
    cores = 4,
    iter = n_samples,
    init = init_values,
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 12
    )
  )

  array_weights <- posterior_summary(fit, pars = "b_sigma") %>%
    as.data.frame() %>%
    mutate(
      array = 1:n_distinct(array_idx),
      weight = 1 / Estimate,
      weight_norm = weight / max(weight)
    )

  return(list(
    fit = fit,
    array_weights = array_weights,
    data = df
  ))
}

plot_results_brms <- function(fit) {
  # Extract data
  df <- fit$data

  # Plot 1: Array weights
  p1 <- ggplot(fit$array_weights,
               aes(x = array, y = Estimate)) +
    geom_point(size = 3) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Array Weights",
         x = "Array",
         y = "Normalized Weight") +
    theme_minimal() +
    ylim(-3, 1)

  # Plot 2: Data and fitted curve
  # Generate prediction data
  x_seq <- seq(min(df$x), max(df$x), length.out = 100)
  pred_data <- data.frame(
    x = x_seq,
    # Use the first array for predictions
    array = factor(rep(levels(df$array)[1], 100))
  )

  # Get predictions
  predictions <- predict(fit$fit, newdata = pred_data)

  pred_df <- data.frame(
    x = x_seq,
    y = predictions[,"Estimate"],
    lower = predictions[,"Q2.5"],
    upper = predictions[,"Q97.5"]
  )

  # Create the second plot
  p2 <- ggplot() +
    # Original data points
    geom_point(data = df,
               aes(x = x, y = y, color = array,
                   alpha = fit$array_weights$weight_norm[as.numeric(array)])) +
    # Fitted line
    geom_line(data = pred_df,
              aes(x = x, y = y),
              color = "red",
              size = 1) +
    # Confidence interval
    geom_ribbon(data = pred_df,
                aes(x = x, ymin = lower, ymax = upper),
                alpha = 0.2,
                fill = "red") +
    labs(title = "Data and Fitted GAM",
         x = "X",
         y = "Y") +
    theme_minimal() +
    theme(legend.position = "right")

  # Arrange plots
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

get_array_errors <- function(model) {
  # Extract posterior samples for sigma parameters
  sigma_samples <- posterior_samples(model$fit, pars = "b_sigma")

  # Calculate summary statistics
  error_summary <- data.frame(
    array = levels(model$data$array),
    error = 1 / colMeans(sigma_samples),  # Convert to weights
    error_lower = 1 / apply(sigma_samples, 2, quantile, probs = 0.975),
    error_upper = 1 / apply(sigma_samples, 2, quantile, probs = 0.025)
  )

  return(error_summary)
}

# Function to make predictions
predict_regression <- function(model, newdata) {
  predictions <- predict(model$fit, newdata = newdata)
  return(predictions[, "Estimate"])
}

# Example usage:
# fit <- bayesian_poly_regression(x, y, array_idx, n_knots = 20)
#
# # Make predictions
# new_data <- data.frame(
#   x = seq(min(x), max(x), length.out = 100)
# )
# predictions <- predict_regression(fit, new_data)
fit <- bayesian_gam_regression(x, y, array_idx, n_knots = 10)
pdf("results/trajectory/20241129_conserved_trajectory_model/20241125_bayes_brms_gam_test2.pdf",width = 8,height = 4)
plot_results_brms(fit)
dev.off()
