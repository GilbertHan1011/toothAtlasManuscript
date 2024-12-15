bayesian_gam_regression_fast <- function(x, y, n_knots = 5, array_idx, n_samples = 2000) {
  # Create data frame
  df <- data.frame(
    y = round(y[!is.na(y)]),  # ensure integers for Poisson
    x = x[!is.na(y)],
    array = factor(array_idx[!is.na(y)])
  )

  # Initial fit with mgcv using Poisson
  init_fit <- mgcv::gam(
    y ~ s(x, bs = "cr") + factor(array_idx),
    family = poisson(link = "log")  # Changed to Poisson
  )
  init_values <- coef(init_fit)

  # Define formula with array-specific smooths
  # Note: For Poisson, we model the log of the mean
  formula <- bf(
    y ~ s(x, bs = "cr", k = 5) + array, # Instead use array-specific smooth I use just array
    # For Poisson, we can still model overdispersion through array-specific effects
    phi ~ 0 + array  # phi parameter for overdispersion
  )

  # Set correct priors
  # Minimal priors
  prior <- c(
    prior(student_t(3, 0, 2), class = "b", dpar = "sigma"),
    prior(gamma(2, 0.1), class = "sd")
  )

  # Fit the model with Poisson family
  fit <- brm(
    formula = formula,
    data = df,
    family = poisson(),
    prior = prior,
    chains = 4,
    cores = 4,
    iter = n_samples,
    init = init_values,
    backend = "cmdstanr",
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 12
    )
  )

  # Extract array-specific weights based on overdispersion parameter
  array_weights <- posterior_summary(fit, pars = "b_phi") %>%
    as.data.frame() %>%
    mutate(
      array = 1:n_distinct(array_idx),
      weight = Estimate,  # For Poisson, larger phi means more reliable
      weight_norm = weight / max(weight)
    )

  return(list(
    fit = fit,
    array_weights = array_weights,
    data = df
  ))
}


geneData1 <- reshapedData[,,1]
prepare_data_for_gam <- function(matrix_data) {
  # Convert matrix to long format
  df <- as.data.frame(matrix_data)
  df$array <- rownames(matrix_data)

  # Reshape to long format
  long_df <- tidyr::pivot_longer(
    df,
    cols = -array,
    names_to = "x",
    values_to = "y"
  )

  # Convert x to numeric
  long_df$x <- as.numeric(long_df$x)  # Changed from long_df$name to long_df$x

  # Return vectors in required format
  list(
    x = long_df$x,
    y = long_df$y,
    array_idx = long_df$array
  )
}


preparedData <- prepare_data_for_gam(geneData1)
fit_count <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)

plot_results_brms(fit_count)
pdf("results/trajectory/20241129_conserved_trajectory_model/20241125_bayes_brms_gam_test_sulf.pdf",width = 8,height = 4)
plot_results_brms(fit_count)
dev.off()


geneData2 <- reshaped_data[,,1]
preparedData <- prepare_data_for_gam(geneData2)
fit_count_data <- bayesian_gam_regression(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
pdf("results/trajectory/20241129_conserved_trajectory_model/20241125_bayes_brms_gam_test_data_sulf_fit.pdf",width = 12,height = 4)
plot_results_brms(fit_count_data)
dev.off()

plot(geneData2[1,])
plot(geneData2[2,])
plot(geneData2[3,])
plot(geneData2[4,])
plot(geneData2[5,])
plot(geneData2[6,])
plot(geneData2[7,])
plot(geneData2[8,])
plot(geneData2[10,])
plot(geneData2[11,])

x_seq <- 1:100
pred_data <- data.frame(
  x = x_seq,
  # Use the first array for predictions
  array = factor(rep("Epi_Chiba", 100))
)

# Get predictions
predictions <- predict(fit_count_data$fit, newdata = pred_data)
estimate <- predictions[,1]
geneData2 <- rbind(estimate,geneData2)
geneData2_scale <- t(scale(t(geneData2)))
Heatmap(geneData2_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)



