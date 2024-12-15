#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
gene_index <- as.numeric(args[1])

# Load libraries
library(brms)
library(tidyverse)
library(mgcv)
library(cmdstanr)

# Source utility functions
source("script/utils/trajectory_model_util.R")

# Load data
inputDataMat <- readRDS("processed_data/trajectory/20241130_data_varmat.Rds")
genes <- dimnames(inputDataMat)[[3]]

# Create output directory if it doesn't exist
dir.create("processed_data/trajectory/model_fits", recursive = TRUE, showWarnings = FALSE)

# Fit model for specific gene
tryCatch({
  # Prepare data
  preparedData <- prepare_data_for_gam(inputDataMat[,,gene_index])
  print(head(preparedData$x))
  # Fit model
  fitModel <- bayesian_gam_regression_nb_shape(preparedData$x,
                                      preparedData$y,
                                      preparedData$array_idx,
                                      n_knots = 5)

  # Save results
  saveRDS(fitModel,
          file = sprintf("processed_data/trajectory/model_fits/gene_%s_fit.Rds",
                         genes[gene_index]))

  # Print success message
  cat(sprintf("Successfully fitted model for gene %s\n", genes[gene_index]))

}, error = function(e) {
  # Log error
  cat(sprintf("Error fitting model for gene %s: %s\n",
              genes[gene_index],
              e$message))
})
