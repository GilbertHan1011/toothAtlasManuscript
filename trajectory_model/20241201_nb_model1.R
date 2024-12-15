#!/usr/bin/env Rscript

# Load libraries
library(brms)
library(tidyverse)
library(mgcv)
library(cmdstanr)
library(foreach)
library(doParallel)
library(progressr)

# Source utility functions
source("script/utils/trajectory_model_util.R")

# Load data
inputDataMat <- readRDS("processed_data/trajectory/20241130_data_varmat.Rds")
genes <- dimnames(inputDataMat)[[3]]

# Create output directory
dir.create("processed_data/trajectory/model_fits", recursive = TRUE, showWarnings = FALSE)

# Setup parallel processing
n_cores <- parallel::detectCores() - 1  # Leave one core free
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Create progress handler
handlers(global = TRUE)
handlers("progress")

# Parallel processing with progress bar
with_progress({
  p <- progressor(along = 1:length(genes))

  results <- foreach(
    gene_index = 1:length(genes),
    .packages = c("brms", "tidyverse", "mgcv", "cmdstanr"),
    .errorhandling = 'pass'
  ) %dopar% {
    tryCatch({
      # Prepare data
      preparedData <- prepare_data_for_gam(inputDataMat[,,gene_index])

      # Fit model
      fitModel <- bayesian_gam_regression_nb_shape(
        preparedData$x,
        preparedData$y,
        preparedData$array_idx,
        n_knots = 5
      )

      # Save results
      saveRDS(
        fitModel,
        file = sprintf(
          "processed_data/trajectory/model_fits/gene_%s_fit.Rds",
          genes[gene_index]
        )
      )

      # Update progress
      p(sprintf("Processed gene %s", genes[gene_index]))

      # Return success status
      list(
        status = "success",
        gene = genes[gene_index],
        message = "Successfully fitted model"
      )

    }, error = function(e) {
      # Return error status
      list(
        status = "error",
        gene = genes[gene_index],
        message = e$message
      )
    })
  }
})

# Stop cluster
stopCluster(cl)

# Process results
success_count <- sum(sapply(results, function(x) x$status == "success"))
error_count <- sum(sapply(results, function(x) x$status == "error"))

# Print summary
cat(sprintf("\nProcessing complete:\n"))
cat(sprintf("Successfully processed: %d genes\n", success_count))
cat(sprintf("Failed processing: %d genes\n", error_count))

# Log errors if any
errors <- results[sapply(results, function(x) x$status == "error")]
if(length(errors) > 0) {
  error_log <- data.frame(
    gene = sapply(errors, function(x) x$gene),
    error = sapply(errors, function(x) x$message)
  )

  write.csv(
    error_log,
    "processed_data/trajectory/model_fits/error_log.csv",
    row.names = FALSE
  )

  cat("\nError details have been saved to 'error_log.csv'\n")
}
