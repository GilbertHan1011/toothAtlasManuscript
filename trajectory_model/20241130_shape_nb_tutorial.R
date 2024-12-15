library(ggplot2)
library(MASS)
library(brms)
library(tidyverse)

# 1. Generate example data with different shape parameters
set.seed(123)

# Function to generate data with different shape parameters
generate_nb_data <- function(n = 1000, mean_value = 10) {
  # Different shape parameters
  shapes <- c(0.5, 2, 10, 50,100)  # from high to low overdispersion

  data <- lapply(shapes, function(shape) {
    rnbinom(n, mu = mean_value, size = shape)
  })

  data.frame(
    value = unlist(data),
    shape = rep(paste("Shape =", shapes), each = n)
  ) %>%
    mutate(shape = factor(shape, levels = paste("Shape =", shapes)))
}

# Generate and plot data
data <- generate_nb_data()

# Plot distributions
p1 <- ggplot(data, aes(x = value, fill = shape)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  facet_wrap(~shape, ncol = 1) +
  theme_minimal() +
  labs(title = "Negative Binomial Distributions with Different Shape Parameters",
       subtitle = "Higher shape = Less overdispersion (closer to Poisson)",
       x = "Count",
       y = "Frequency")

# Calculate summary statistics
summary_stats <- data %>%
  group_by(shape) %>%
  summarise(
    mean = mean(value),
    variance = var(value),
    variance_mean_ratio = var(value)/mean(value)
  )

print(summary_stats)

# 2. Time series example with different shapes
set.seed(123)
n_samples <- 100
x <- seq(1, n_samples)

# Generate counts with different shapes
generate_nb_ts_data <- function(x, shape) {
  # Base trend (sine wave)
  mu <- 10 * sin(x/10) + 20

  # Generate NB data with specified shape
  counts <- rnbinom(length(x), mu = mu, size = shape)
  counts
}

# Generate data with different shapes
shapes <- c(0.5, 2, 10, 50)
ts_data <- data.frame(
  x = rep(x, length(shapes)),
  counts = unlist(lapply(shapes, function(s) generate_nb_ts_data(x, s))),
  shape = rep(paste("Shape =", shapes), each = n_samples)
) %>%
  mutate(shape = factor(shape, levels = paste("Shape =", shapes)))

# Plot time series
p2 <- ggplot(ts_data, aes(x = x, y = counts)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~shape, ncol = 1) +
  theme_minimal() +
  labs(title = "Time Series with Different Shape Parameters",
       subtitle = "Higher shape = Less variation around trend",
       x = "Position",
       y = "Counts")

# 3. Fit models with different shapes
# Take one set of data and fit with different shape specifications
test_data <- ts_data[ts_data$shape == "Shape = 2", ]

# Fit model with estimated shape
model_free <- brm(
  counts ~ s(x, bs = "cr"),
  data = test_data,
  family = negbinomial(),
  chains = 2,
  iter = 1000
)

# Fit model with fixed shape
model_fixed <- brm(
  bf(counts ~ s(x, bs = "cr"),
     shape ~ 1),  # Fixed shape
  data = test_data,
  family = negbinomial(),
  chains = 2,
  iter = 1000
)

# Print results
print("Model with estimated shape:")
summary(model_free)
print("\nModel with fixed shape:")
summary(model_fixed)

# Save plots
ggsave("results/trajectory/20241129_conserved_trajectory_model/shape_distributions.pdf", p1, height = 10, width = 8)
ggsave("results/trajectory/20241129_conserved_trajectory_model/shape_timeseries.pdf", p2, height = 10, width = 8)
