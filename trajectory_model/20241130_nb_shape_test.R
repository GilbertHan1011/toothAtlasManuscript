bayesian_gam_regression_nb_shape <- function(x, y, n_knots = 5, array_idx, n_samples = 2000) {
  # Create data frame
  df <- data.frame(
    y = round(y[!is.na(y)]),  # ensure integers
    x = x[!is.na(y)],
    array = factor(array_idx[!is.na(y)])
  )

  # # Initial fit with mgcv using negative binomial
  # init_fit <- mgcv::gam(
  #   y ~ s(x, bs = "cr") + factor(array_idx),
  #   family = nb()  # Changed to negative binomial
  # )
  # init_values <- coef(init_fit)

  # Define formula with array-specific shape parameters
  formula <- bf(
    y ~ s(x, bs = "cr", k = 5) + array,
    shape ~ 0 + array  # shape parameter varies by array
  )

  # Set priors
  prior <- c(
    prior(normal(0, 5), class = "b"),  # prior for fixed effects
    prior(normal(0, 2), class = "b", dpar = "shape"),  # prior for shape parameters
    prior(normal(0, 2), class = "sds")  # prior for smooth terms
  )

  # Fit the model with negative binomial family
  fit <- brm(
    formula = formula,
    data = df,
    family = negbinomial(),
    prior = prior,
    chains = 4,
    cores = 4,
    iter = n_samples,
    #init = init_values,
    backend = "cmdstanr",
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 12
    )
  )

  # Extract array-specific weights based on shape parameter
  # Higher shape = less overdispersion = more reliable
  array_weights <- posterior_summary(fit, pars = "b_shape") %>%
    as.data.frame() %>%
    mutate(
      array = 1:n_distinct(array_idx),
      shape = exp(Estimate),  # Convert from log scale
      weight = shape,  # Higher shape means more reliable
      weight_norm = weight / max(weight)
    )

  # Calculate additional diagnostics
  diagnostics <- df %>%
    group_by(array) %>%
    summarise(
      mean = mean(y),
      variance = var(y),
      overdispersion = variance/mean,
      n_obs = n()
    )

  return(list(
    fit = fit,
    array_weights = array_weights,
    diagnostics = diagnostics,
    data = df
  ))
}


#fit_count_2_Prex2 <- bayesian_gam_regression_fast(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
inputDataMat <- readRDS("processed_data/trajectory/20241130_data_varmat.Rds")
preparedData <- prepare_data_for_gam(inputDataMat[,,431])
fit_count_3_Bglap2 <- bayesian_gam_regression_nb_shape(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)


plot_results_brms(fit_count_3_Bglap2)


preparedData <- prepare_data_for_gam(inputDataMat[,,1])
fit_count_3_gene1 <- bayesian_gam_regression_nb_shape(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)
plot_results_brms(fit_count_3_gene1)

preparedData <- prepare_data_for_gam(inputDataMat[,,2])
fit_count_3_gene2 <- bayesian_gam_regression_nb_shape(preparedData$x, preparedData$y, preparedData$array_idx, n_knots = 5)

plot_results_brms(fit_count_3_gene2)

Heatmap(inputDataMat[,,2],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)

plot_results_brms(fit_count_2_Bglap2)
geneData3_scale <- t(scale(t(inputDataMat[,,431])))
Heatmap(geneData3_scale,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
Heatmap(inputDataMat[,,431],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
