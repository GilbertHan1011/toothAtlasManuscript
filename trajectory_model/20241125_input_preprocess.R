mes <- readRDS("processed_data/integrated_data/20241106_mesenchyme.Rds")
pseudo <- read.csv("processed_data/trajectory/20241124_pseudotime_predicted.csv",row.names = 1)
varGene <- read.csv("processed_data/framework/geneMeta/20241125_mes_vargene_2000.csv",row.names = 1) %>% unlist()
mes <- mes[varGene,rownames(pseudo)]
saveRDS(mes,"processed_data/integrated_data/20241125_mes_var.Rds")
mesDf <- mes@assays$originalexp@data %>% as.matrix()
pseudotime <- pseudo$lightGBM
batch <- mes@meta.data$Project
# Function to bin values into 100 regions
bin_to_100 <- function(x) {
  # Create 100 bins
  bins <- cut(x,
              breaks = seq(min(x), max(x), length.out = 101), # 101 breaks make 100 bins
              labels = FALSE,  # Return numeric labels
              include.lowest = TRUE  # Include the lowest value in first bin
  )
  return(bins)
}

# Apply to your pseudotime
pseudotime_binned <- bin_to_100(pseudotime)
plot(pseudotime,pseudotime_binned)
metaDf <- data.frame(batch, pseudotime_binned)
metaDf$bin <- paste0(metaDf$batch,"_",metaDf$pseudotime_binned)


#== calculate means of bins------------
calculate_bin_means_fast1 <- function(expression_matrix, bin_labels) {
  # Convert to factors once
  bin_factors <- factor(bin_labels, levels = sort(unique(bin_labels)))

  # Use apply + tapply
  result_matrix <- t(apply(expression_matrix, 1, function(x) {
    tapply(x, bin_factors, mean)
  }))

  return(result_matrix)
}
binned_means <- calculate_bin_means_fast1(mesDf, metaDf$bin)

parts <- strsplit(colnames(binned_means), "_")
prefixes <- sapply(parts, function(x) paste(x[1:2], collapse="_"))  # Gets "Atlas_Jan"
numbers <- as.numeric(sapply(parts, function(x) x[3]))  # Gets 1, 10, 11, etc.
plot(numbers,binned_means["Bglap",])
plot(numbers,binned_means["Sp6",])

thred = 0.1
geneNum <- round(thred * ncol(binned_means))

filteredGene <- rownames(binned_means)[rowSums(binned_means>0)>geneNum]

batch_thred <- 0.3
n_bin = 100
bath_thred_real <- batch_thred*n_bin
batchName <- names(table(prefixes) >bath_thred_real)[table(prefixes) >bath_thred_real]
binned_means_filter <- binned_means[filteredGene,prefixes%in%batchName]

reshape_to_3d <- function(matrix_data, prefixes, numbers, n_bins=100) {
  # Get unique prefixes
  unique_prefixes <- unique(prefixes)

  # Create empty array with NAs
  result <- array(NA,
                  dim = c(length(unique_prefixes), n_bins, nrow(matrix_data)),
                  dimnames = list(unique_prefixes,        # prefixes
                                  1:n_bins,                 # numbers
                                  rownames(matrix_data)))   # genes

  # Fill the array
  for(i in seq_along(prefixes)) {
    prefix <- prefixes[i]
    number <- numbers[i]
    if(number <= n_bins) {  # Check if number is within bounds
      result[prefix, number, ] <- matrix_data[, i]
    }
  }

  return(result)
}

# Reshape the data


parts <- strsplit(colnames(binned_means_filter), "_")
prefixes <- sapply(parts, function(x) paste(x[1:2], collapse="_"))  # Gets "Atlas_Jan"
numbers <- as.numeric(sapply(parts, function(x) x[3]))  # Gets 1, 10, 11, etc.

# Check the data before reshaping
print("Range of numbers:")
print(range(numbers))
print("Unique prefixes:")
print(unique(prefixes))
print("Matrix dimensions:")
print(dim(binned_means_filter))

# Reshape the data
#reshaped_data <- reshape_to_3d(binned_means_filter, prefixes, numbers)

reshaped_data <- reshape_to_3d(binned_means_filter, prefixes, numbers)

mes$pseudo <- pseudo$lightGBM
source("script/utils/trajectory_model_util.R")
results <- scRNA_2_mat(mes,assay = "originalexp",slot = "count",pseudo_col = "pseudo",project_col = "Project")

reshapedData <- results$reshaped_data


