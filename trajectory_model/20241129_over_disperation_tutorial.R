library(ggplot2)
library(MASS)
library(brms)

# 1. Generate example data to demonstrate overdispersion
set.seed(123)

# Function to generate data with different levels of dispersion
generate_count_data <- function(n = 1000, mean_value = 10) {
  # Poisson data (no overdispersion)
  poisson_data <- rpois(n, lambda = mean_value)

  # Negative Binomial data (overdispersed)
  # smaller size parameter = more overdispersion
  nb_data_high <- rnbinom(n, mu = mean_value, size = 1)   # high overdispersion
  nb_data_low <- rnbinom(n, mu = mean_value, size = 5)    # low overdispersion

  # Combine into data frame
  data.frame(
    value = c(poisson_data, nb_data_low, nb_data_high),
    type = rep(c("Poisson", "NB (Low OD)", "NB (High OD)"), each = n)
  )
}

# Generate and plot data
data <- generate_count_data()

# Plot distributions
ggplot(data, aes(x = value, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  facet_wrap(~type, ncol = 1) +
  theme_minimal() +
  labs(title = "Comparison of Count Distributions",
       subtitle = "Showing different levels of overdispersion",
       x = "Count",
       y = "Frequency")
ggsave("results/trajectory/20241129_conserved_trajectory_model/20241129_overdisperion_show.pdf")
# Calculate mean and variance for each distribution
summary_stats <- data %>%
  group_by(type) %>%
  summarise(
    mean = mean(value),
    variance = var(value),
    variance_mean_ratio = var(value)/mean(value)
  )

print(summary_stats)

# 2. Real-world example with RNA-seq like data
set.seed(123)

# Generate example gene expression data
n_samples <- 100
x <- seq(1, n_samples)

# Generate counts with different dispersion levels
generate_expression_data <- function(x, dispersion = "low") {
  # Base trend
  mu <- 10 * sin(x/10) + 20

  if(dispersion == "none") {
    # Poisson (no overdispersion)
    counts <- rpois(length(x), lambda = mu)
  } else if(dispersion == "low") {
    # Low overdispersion
    counts <- rnbinom(length(x), mu = mu, size = 5)
  } else {
    # High overdispersion
    counts <- rnbinom(length(x), mu = mu, size = 1)
  }
  counts
}

# Generate data
expression_data <- data.frame(
  x = rep(x, 3),
  counts = c(
    generate_expression_data(x, "none"),
    generate_expression_data(x, "low"),
    generate_expression_data(x, "high")
  ),
  dispersion = rep(c("No OD (Poisson)", "Low OD", "High OD"), each = n_samples)
)

# Plot the data
ggplot(expression_data, aes(x = x, y = counts)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~dispersion, ncol = 1) +
  theme_minimal() +
  labs(title = "Gene Expression-like Data with Different Dispersion Levels",
       x = "Position",
       y = "Counts")

# 3. Fit models to compare Poisson vs Negative Binomial
# Take high overdispersion data
high_od_data <- expression_data[expression_data$dispersion == "High OD", ]

# Fit Poisson model
poisson_model <- brm(
  counts ~ s(x, bs = "cr"),
  data = high_od_data,
  family = poisson(),
  chains = 2,
  iter = 1000
)

# Fit Negative Binomial model
nb_model <- brm(
  counts ~ s(x, bs = "cr"),
  data = high_od_data,
  family = negbinomial(),
  chains = 2,
  iter = 1000
)

# Compare models
print("Poisson Model Summary:")
summary(poisson_model)
print("\nNegative Binomial Model Summary:")
summary(nb_model)

# Plot predictions
plot_predictions <- function(model, data, title) {
  predictions <- predict(model)

  ggplot() +
    geom_point(data = data, aes(x = x, y = counts), alpha = 0.5) +
    geom_line(aes(x = data$x, y = predictions[,"Estimate"]), color = "red") +
    geom_ribbon(aes(x = data$x,
                    ymin = predictions[,"Q2.5"],
                    ymax = predictions[,"Q97.5"]),
                alpha = 0.2, fill = "red") +
    theme_minimal() +
    labs(title = title)
}

# Plot both models
p1 <- plot_predictions(poisson_model, high_od_data, "Poisson Model Fit")
p2 <- plot_predictions(nb_model, high_od_data, "Negative Binomial Model Fit")

gridExtra::grid.arrange(p1, p2, ncol = 1)
