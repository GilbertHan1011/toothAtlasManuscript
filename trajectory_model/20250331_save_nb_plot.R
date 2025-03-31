#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)

# 1. Generate example data with different shape parameters
set.seed(123)

# Function to generate data with different shape parameters
generate_nb_data <- function(n = 1000, mean_value = 10) {
  # Different shape parameters
  shapes <- c(0.5, 2, 10, 50, 100)  # from high to low overdispersion

  data <- lapply(shapes, function(shape) {
    rnbinom(n, mu = mean_value, size = shape)
  })

  data.frame(
    value = unlist(data),
    shape = rep(paste("Shape =", shapes), each = n)
  ) %>%
    mutate(shape = factor(shape, levels = paste("Shape =", shapes)))
}

# Generate data
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

# Display the plot
print(p1)

# Save the plot to a file
# You can change the format by using different file extensions: .png, .pdf, .jpg, etc.
ggsave("../../results/trajectory/20241129_conserved_trajectory_model/negative_binomial_distributions.pdf", plot = p1, width = 10, height = 12, dpi = 300)
cat("Plot saved as 'negative_binomial_distributions.png'\n") 