# Script to create correlation plots from benchmark results
# --------------------------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Source the plotting functions
source("trajectory_model/plot_correlation_by_noise.R")

# Create output directory
dir.create("results/trajectory", recursive = TRUE, showWarnings = FALSE)

# Load benchmark results (assuming they're saved as RDS)
# If you have the results loaded as benchmark_results_static1, 
# comment this line out and use that variable directly
benchmark_results <- readRDS("results/trajectory/benchmark_results.rds")

# Create correlation by noise ratio plot
cat("Creating correlation by noise ratio plot...\n")
ratio_plot <- plot_correlation_by_noise(
  benchmark_results,
  "results/trajectory/correlation_by_noise_ratio.pdf"
)

# Create correlation by noise level plot
cat("Creating correlation by noise level plot...\n")
level_plot <- plot_correlation_by_noise_level(
  benchmark_results,
  "results/trajectory/correlation_by_noise_level.pdf"
)

# If you want to directly compare performance between methods
# You can create a customized direct comparison plot

# First, prepare the data
processed_data <- benchmark_results$summary %>%
  # For each unique simulation configuration
  group_by(simulation, signal_noise_sd, noise_sd) %>%
  # Only keep one entry per configuration (since bayes is the same for both)
  summarize(
    signal_noise_sd = first(signal_noise_sd),
    noise_sd = first(noise_sd),
    bayes_corr = first(bayes_correlation),
    pearson_corr = cor_correlation[method == "pearson"],
    dtw_corr = cor_correlation[method == "dtw"],
    .groups = "drop"
  ) %>%
  # Calculate the performance difference
  mutate(
    pearson_vs_bayes = pearson_corr - bayes_corr,
    dtw_vs_bayes = dtw_corr - bayes_corr,
    dtw_vs_pearson = dtw_corr - pearson_corr,
    noise_ratio = signal_noise_sd / noise_sd
  )

# Create a plot showing the advantage of each method
cat("Creating method comparison plot...\n")
comparison_data <- processed_data %>%
  pivot_longer(
    cols = c(pearson_vs_bayes, dtw_vs_bayes, dtw_vs_pearson),
    names_to = "comparison",
    values_to = "advantage"
  ) %>%
  mutate(
    comparison_label = case_when(
      comparison == "pearson_vs_bayes" ~ "Pearson vs Bayes",
      comparison == "dtw_vs_bayes" ~ "DTW vs Bayes",
      comparison == "dtw_vs_pearson" ~ "DTW vs Pearson"
    )
  )

comparison_plot <- ggplot(comparison_data, 
                         aes(x = noise_ratio, y = advantage, color = comparison_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_smooth(method = "loess", se = TRUE) +
  geom_point() +
  facet_wrap(~signal_noise_sd, 
             labeller = labeller(signal_noise_sd = function(x) paste("Signal Noise:", x))) +
  labs(
    title = "Method Advantage by Noise Ratio",
    subtitle = "Positive values indicate the first method outperforms the second",
    x = "Noise Ratio (signal_noise/noise)",
    y = "Correlation Advantage",
    color = "Comparison"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  scale_color_brewer(palette = "Set1")

# Save the comparison plot
ggsave("results/trajectory/method_comparison.pdf", comparison_plot, width = 10, height = 7)

cat("All plots created successfully! Check the results/trajectory directory.\n") 