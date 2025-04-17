sim_data <- generate_simulation_data(
  n_true_arrays = 25,
  n_noise_arrays = 15,
  n_positions = 100,
  base_expression = 1.5,
  signal_strength = 10,
  signal_noise_sd = 1,
  noise_sd = 7,
  seed = 42
)



# Calculate correlations - handle potential errors
cor_results <- 
  calculate_trajectory_correlations(
    trajectories = sim_data$data,
    method = "dtw",
    weights = NULL,  # No longer using array weights since the function now creates its own
    smooth_first =F,
    benchmark = TRUE,
    true_signal = sim_data$true_signal,
    true_arrays = sim_data$true_arrays
  )
dir.create("results/trajectory/benchmark_plot")
pdf("results/trajectory/benchmark_plot/20250415_pairwise_similarity.pdf")
Heatmap(cor_results$correlation_matrix,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,col=colorRamp2(c(0, 0.5, 1), c("DeepSkyBlue3", "white", "red")))
dev.off()

rowMeanCor <- rowMeans(cor_results$correlation_matrix)
barplot(rowMeanCor)

# Convert the means to a data frame for ggplot
plot_data <- data.frame(
  Variable = names(rowMeanCor),
  Mean_Correlation = rowMeanCor,
  id = 1:40,
  group = ifelse(1:40 <= 25, "Group 1", "Group 2")  # Add grouping variable
)

ggplot(plot_data, aes(x = reorder(Variable, id), 
                      y = Mean_Correlation,
                      fill = group)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Group 1" = "#FF4444", "Group 2" = "#4285f4")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x  = element_blank()  # Remove x-axis label
  ) +
  labs(
    title = "Mean Correlations by Variable",
    y = "Mean Correlation"  # Removed x = "Variables"
  )
ggsave("results/trajectory/benchmark_plot/20250415_benchmark_hist.pdf",width = 6,height = 4)


cor_results2 <- 
  calculate_trajectory_correlations(
    trajectories = sim_data$data,
    method = "pearson",
    weights = NULL,  # No longer using array weights since the function now creates its own
    smooth_first =F,
    benchmark = TRUE,
    true_signal = sim_data$true_signal,
    true_arrays = sim_data$true_arrays
  )


pdf("results/trajectory/benchmark_plot/20250415_pairwise_similarity_pearson.pdf")
Heatmap(cor_results2$correlation_matrix,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,col=colorRamp2(c(0, 0.5, 1), c("DeepSkyBlue3", "white", "red")))
dev.off()



pdf("results/trajectory/benchmark_plot/20250415_pairwise_similarity_pearson.pdf")
Heatmap(cor_results2$correlation_matrix,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,col=colorRamp2(c(0, 0.5, 1), c("DeepSkyBlue3", "white", "red")))
dev.off()

rowMeanCor <- rowMeans(cor_results2$correlation_matrix)
barplot(rowMeanCor)

# Convert the means to a data frame for ggplot
plot_data <- data.frame(
  Variable = names(rowMeanCor),
  Mean_Correlation = rowMeanCor,
  id = 1:40,
  group = ifelse(1:40 <= 25, "Group 1", "Group 2")  # Add grouping variable
)

ggplot(plot_data, aes(x = reorder(Variable, id), 
                      y = Mean_Correlation,
                      fill = group)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Group 1" = "#FF4444", "Group 2" = "#4285f4")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x  = element_blank()  # Remove x-axis label
  ) +
  labs(
    title = "Mean Correlations by Variable",
    y = "Mean Correlation"  # Removed x = "Variables"
  )
ggsave("results/trajectory/benchmark_plot/20250415_benchmark_hist2.pdf",width = 6,height = 4)






cor_end_time <- Sys.time()
cor_runtime <- as.numeric(difftime(cor_end_time, cor_start_time, units = "secs"))

# Store results
result_name <- sprintf("sim_%d_signal_noise_%.2f_noise_sd_%.2f_method_%s_smooth_%s",
                       sim, signal_noise, noise_sd, method, smooth)
