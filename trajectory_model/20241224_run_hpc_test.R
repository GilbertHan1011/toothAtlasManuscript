
# Load libraries
library(brms)
library(tidyverse)
library(mgcv)
library(cmdstanr)
library(foreach)
library(doParallel)
library(progressr)

# Source utility functions
source("~/hlt_data/traj_util.R")

# Load data
inputDataMat <- readRDS("~/hlt_data/data/20241223_input_osteogenicData_reshaped.Rds")
genes <- dimnames(inputDataMat)[[3]]

sample_name = dimnames(inputDataMat)[[1]]

Mes_index = grep("Mes_",sample_name)
mesMat = inputDataMat[Mes_index,,]
# Create output directory
# Setup parallel processing
n_cores <- 25  # Leave one core free
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
      preparedData <- prepare_data_for_gam(mesMat[,,gene_index])

      # Fit model
      fitModel <- bayesian_gam_regression_nb_shape(
        preparedData$x,
        preparedData$y,
        preparedData$array_idx,
        n_knots = 5
      )
      weight = fitModel$array_weights
      # Save results
      write.csv(
        weight,
        file = sprintf(
          "~/hlt_data/result/20241224_test/gene_%s_fit_weight.csv",
          genes[gene_index]
        )
      )
      saveRDS(
        fitModel,
        file = sprintf(
          "~/hlt_data/result/20241224_test/gene_%s_fit_weight.Rds",
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
