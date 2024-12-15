bayesian_gam_regression_nb_shape <- function(x, y, n_knots = 5, array_idx, n_samples = 2000) {
  # Create data frame
  df <- data.frame(
    y = round(y[!is.na(y)]),  # ensure integers
    x = x[!is.na(y)],
    array = factor(array_idx[!is.na(y)])
  )

  # Initial fit with mgcv using negative binomial
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



scRNA_2_mat <- function(mes, assay ,slot ,pseudo_col, project_col, thred = 0.1, batch_thred = 0.3, n_bin = 100) {
  # Extract expression matrix
  mesDf <- GetAssayData(mes,assay = assay,slot = slot) %>% as.matrix()
  pseudotime <- mes@meta.data[[pseudo_col]]
  batch <- mes@meta.data[[project_col]]

  # Function to bin values
  bin_to_100 <- function(x, n_bins = n_bin) {
    bins <- cut(x,
                breaks = seq(min(x), max(x), length.out = n_bins + 1),
                labels = FALSE,
                include.lowest = TRUE)
    return(bins)
  }

  # Bin pseudotime
  pseudotime_binned <- bin_to_100(pseudotime)
  metaDf <- data.frame(batch, pseudotime_binned)
  metaDf$bin <- paste0(metaDf$batch, "_", metaDf$pseudotime_binned)

  # Calculate bin means
  calculate_bin_means_fast1 <- function(expression_matrix, bin_labels) {
    bin_factors <- factor(bin_labels, levels = sort(unique(bin_labels)))
    result_matrix <- t(apply(expression_matrix, 1, function(x) {
      tapply(x, bin_factors, mean)
    }))
    return(result_matrix)
  }

  binned_means <- calculate_bin_means_fast1(mesDf, metaDf$bin)

  # Filter genes and batches
  geneNum <- round(thred * ncol(binned_means))
  filteredGene <- rownames(binned_means)[rowSums(binned_means > 0) > geneNum]

  #parts <- strsplit(colnames(binned_means), "_")
  #prefixes <- sapply(parts, function(x) paste(x[1:2], collapse = "_"))
  prefixes <- sapply(strsplit(colnames(binned_means), "_"),
                     function(x) paste(x[-length(x)], collapse = "_"))
  numbers <- sapply(strsplit(colnames(binned_means), "_"),
                    function(x) paste(x[length(x)], collapse = "_")) %>% as.numeric()


  bath_thred_real <- batch_thred * n_bin
  batchName <- names(table(prefixes) > bath_thred_real)[table(prefixes) > bath_thred_real]
  binned_means_filter <- binned_means[filteredGene, prefixes %in% batchName]

  # Reshape to 3D array
  reshape_to_3d <- function(matrix_data, prefixes, numbers, n_bins = n_bin) {
    unique_prefixes <- unique(prefixes)
    result <- array(NA,
                    dim = c(length(unique_prefixes), n_bins, nrow(matrix_data)),
                    dimnames = list(unique_prefixes,
                                    1:n_bins,
                                    rownames(matrix_data)))

    for(i in seq_along(prefixes)) {
      prefix <- prefixes[i]
      number <- numbers[i]
      if(number <= n_bins) {
        result[prefix, number, ] <- matrix_data[, i]
      }
    }
    return(result)
  }

  # Process final data
  prefixes <- sapply(strsplit(colnames(binned_means_filter), "_"),
                     function(x) paste(x[-length(x)], collapse = "_"))
  numbers <- sapply(strsplit(colnames(binned_means_filter), "_"),
                    function(x) paste(x[length(x)], collapse = "_")) %>% as.numeric()
  reshaped_data <- reshape_to_3d(binned_means_filter, prefixes, numbers)

  # Return results
  return(list(
    reshaped_data = reshaped_data,
    binned_means = binned_means,
    binned_means_filter = binned_means_filter,
    filtered_genes = filteredGene,
    batch_names = batchName,
    metadata = metaDf
  ))
}

# Example usage:
# result <- process_scRNA_data(
#   mes = mes,
#   lightGBM = pseudo$lightGBM,
#   Project = mes@meta.data$Project,
#   thred = 0.1,
#   batch_thred = 0.3,
#   n_bin = 100
# )


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
    ylim(min(fit$array_weights$Estimate), max(fit$array_weights$Estimate))

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
    theme(legend.position = "none")

  # Arrange plots
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}


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
