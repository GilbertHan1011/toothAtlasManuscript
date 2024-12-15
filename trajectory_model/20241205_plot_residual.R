library(ggplot2)
library(gridExtra)
library(mgcv)

# Assuming your GAM model is called 'gam_model'

# 1. Create diagnostic plots
plot_residuals <- function(model) {
  # Get model diagnostics
  fitted_vals <- fitted(model)
  residuals <- residuals(model, type = "deviance")
  se <- sqrt(model$sig2)  # Standard errors

  # Create data frame for plotting
  df_diagnostics <- data.frame(
    fitted = fitted_vals,
    residuals = residuals,
    observed = fitted_vals + residuals,
    se = se
  )

  # Residuals vs Fitted
  p1 <- ggplot(df_diagnostics, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    labs(title = "Residuals vs Fitted",
         x = "Fitted values",
         y = "Residuals")

  # Fitted vs Observed
  p2 <- ggplot(df_diagnostics, aes(x = observed, y = fitted)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = "Fitted vs Observed",
         x = "Observed values",
         y = "Fitted values")

  # Residuals QQ Plot
  p3 <- ggplot(df_diagnostics, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = "Normal Q-Q Plot",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles")

  # Scale-Location Plot
  p4 <- ggplot(df_diagnostics, aes(x = fitted, y = sqrt(abs(residuals)))) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    labs(title = "Scale-Location Plot",
         x = "Fitted values",
         y = "√|Residuals|")

  # Combine all plots
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

# Use the function
plot_residuals(gam_model)


# -------- Example 1: Simulated Data --------
set.seed(123)

# Generate sample data
n <- 200
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
# Generate response with non-linear relationships
f1 <- function(x) sin(2 * pi * x)
f2 <- function(x) exp(2 * x)
y <- f1(x1) + f2(x2) + rnorm(n, 0, 0.3)

# Fit GAM model
gam_model <- gam(y ~ s(x1) + s(x2), method = "REML")

# Basic summary
summary(gam_model)

# Plot diagnostics
par(mfrow = c(2, 2))
gam.check(gam_model)

# Plot smooth terms
plot(gam_model, pages = 1)



library(ggplot2)
library(gridExtra)
library(dplyr)

# Create sample data
set.seed(123)
n <- 100
x <- seq(1, n, 1)
true_mean <- 10
true_var <- 4

# Generate data with known mean and variance
y <- rnorm(n, mean = true_mean + 0.05*x, sd = sqrt(true_var))
df <- data.frame(x = x, y = y)

# Fit a simple linear model
fit <- lm(y ~ x, data = df)
df$fitted <- fitted(fit)
df$residuals <- residuals(fit)
df$mean_line <- mean(y)

# 1. Mean Plot
p1 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(y), color = "red", size = 1) +
  annotate("text", x = max(x), y = mean(y),
           label = "Mean", color = "red", hjust = 1) +
  labs(title = "Mean: Average of all y values",
       x = "X", y = "Y") +
  theme_minimal()

# 2. Variance Plot
p2 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(y), color = "red", size = 1) +
  geom_segment(aes(xend = x, yend = mean_line),
               color = "blue", alpha = 0.3) +
  labs(title = "Variance: Spread around the mean",
       x = "X", y = "Y") +
  theme_minimal()

# 3. Residuals Plot
p3 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fitted), color = "red", size = 1) +
  geom_segment(aes(xend = x, yend = fitted),
               color = "blue", alpha = 0.3) +
  labs(title = "Residuals: Distance from fitted line",
       x = "X", y = "Y") +
  theme_minimal()

# 4. Fitness Plot
p4 <- ggplot(df, aes(x = y, y = fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1,
              color = "red", linetype = "dashed") +
  labs(title = "Fitness: Observed vs Fitted values",
       x = "Observed Values", y = "Fitted Values") +
  theme_minimal()

# Combine all plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Alternative: More detailed single plot showing all components
p_combined <- ggplot(df, aes(x = x, y = y)) +
  # Original data points
  geom_point(alpha = 0.5) +
  # Fitted line
  geom_line(aes(y = fitted), color = "red", size = 1) +
  # Mean line
  geom_hline(yintercept = mean(y), color = "blue",
             linetype = "dashed") +
  # Residuals
  geom_segment(aes(xend = x, yend = fitted),
               color = "green", alpha = 0.3) +
  # Variance from mean
  geom_segment(aes(xend = x, yend = mean_line),
               color = "purple", alpha = 0.3) +
  # Add labels
  annotate("text", x = max(x), y = mean(y),
           label = "Mean", color = "blue", hjust = 1) +
  annotate("text", x = max(x), y = max(fitted),
           label = "Fitted Line", color = "red", hjust = 1) +
  labs(title = "All Components in One Plot",
       subtitle = "Green segments = Residuals, Purple segments = Variance",
       x = "X", y = "Y") +
  theme_minimal()

# Display the combined plot
print(p_combined)

# Add statistical summary
summary_stats <- data.frame(
  Metric = c("Mean", "Variance", "Residual SD", "R-squared"),
  Value = c(
    mean(df$y),
    var(df$y),
    sd(df$residuals),
    summary(fit)$r.squared
  )
)

# Print summary statistics
print(summary_stats)


library(ggplot2)
library(gridExtra)
library(dplyr)

# Create sample data
set.seed(123)
n <- 100
x <- seq(1, n, 1)
true_mean <- 10
true_sd <- 2

# Generate data
y <- rnorm(n, mean = true_mean + 0.05*x, sd = true_sd)
df <- data.frame(x = x, y = y)

# Fit model and calculate statistics
fit <- lm(y ~ x, data = df)
df$fitted <- fitted(fit)
df$residuals <- residuals(fit)
df$mean_line <- mean(y)

# Calculate standard deviation bands
df$sd_upper <- mean(y) + sd(y)
df$sd_lower <- mean(y) - sd(y)
df$sd_upper2 <- mean(y) + 2*sd(y)
df$sd_lower2 <- mean(y) - 2*sd(y)

# 1. Mean and Standard Deviation Plot
p1 <- ggplot(df, aes(x = x, y = y)) +
  # Add SD bands
  geom_ribbon(aes(ymin = sd_lower2, ymax = sd_upper2),
              fill = "lightblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = sd_lower, ymax = sd_upper),
              fill = "lightblue", alpha = 0.3) +
  # Add points and mean line
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(y), color = "red", size = 1) +
  geom_hline(yintercept = df$sd_upper, color = "blue",
             linetype = "dashed") +
  geom_hline(yintercept = df$sd_lower, color = "blue",
             linetype = "dashed") +
  # Add labels
  annotate("text", x = max(x), y = mean(y),
           label = "Mean", color = "red", hjust = 1) +
  annotate("text", x = max(x), y = df$sd_upper,
           label = "Mean + SD", color = "blue", hjust = 1) +
  annotate("text", x = max(x), y = df$sd_lower,
           label = "Mean - SD", color = "blue", hjust = 1) +
  labs(title = "Mean and Standard Deviation",
       subtitle = "Light blue bands show ±1 and ±2 SD",
       x = "X", y = "Y") +
  theme_minimal()

# 2. Variance Plot
p2 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(y), color = "red", size = 1) +
  geom_segment(aes(xend = x, yend = mean_line),
               color = "blue", alpha = 0.3) +
  labs(title = "Variance",
       subtitle = "Blue lines show distance from mean",
       x = "X", y = "Y") +
  theme_minimal()

# 3. Residuals Plot
p3 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fitted), color = "red", size = 1) +
  geom_segment(aes(xend = x, yend = fitted),
               color = "green", alpha = 0.3) +
  labs(title = "Residuals",
       subtitle = "Green lines show distance from fitted line",
       x = "X", y = "Y") +
  theme_minimal()

# 4. Model Fit Plot
p4 <- ggplot(df, aes(x = y, y = fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1,
              color = "red", linetype = "dashed") +
  labs(title = "Model Fit",
       subtitle = "Observed vs Fitted values",
       x = "Observed Values", y = "Fitted Values") +
  theme_minimal()

# Combine all plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Print summary statistics
summary_stats <- data.frame(
  Metric = c("Mean", "Standard Deviation", "Variance",
             "Residual SD", "R-squared"),
  Value = c(
    mean(df$y),
    sd(df$y),
    var(df$y),
    sd(df$residuals),
    summary(fit)$r.squared
  )
)

print(summary_stats)
