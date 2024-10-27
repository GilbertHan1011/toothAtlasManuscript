#== annotation utils----------------
#> Author: Gilbert Han
#> 2023.4.9
#>

#=== plot percentage bar plot---------------------------------------

percentPlotFun <- function(object, cluster_id, group_id, group_levels=NULL,color=brewer.pal(8, "Spectral")){
  require(RColorBrewer)
  require(scales)
  require(tidyverse)
  df <- table(object@meta.data[[cluster_id]], object@meta.data[[group_id]]) %>% t() %>% as.matrix()
  df <- t(apply(df, 1, function(x) x / sum(x)))
  df <- as.data.frame(df)

  df_long <- df %>%
    rownames_to_column(group_id) %>%
    pivot_longer(cols = -group_id, names_to = cluster_id)
  num <- length(unique(object@meta.data[[group_id]]))
  mycolor <- colorRampPalette(color)(num)

  if (!is.null(group_levels)){
    df_long[[group_id]] <- factor(df_long[[group_id]],levels = group_levels)
  }

  # Create bar plot
  ggplot(df_long, aes(x = !!sym(cluster_id), y = value, fill = !!sym(group_id))) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = percent_format()) +
    labs(x = "Cluster", y = "Age") +
    scale_fill_manual(values = mycolor) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 15, face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1, size = 15, face = "bold"),
          title = NULL, panel.background = element_blank())
}

#== make annoation tree------------------------

anno_row <- function(object,cluster_id, group_id, threshold_add=0.2) {
  # Create the table
  df <- table(object@meta.data[[cluster_id]], object@meta.data[[group_id]]) %>% t() %>% as.matrix()

  # Normalize the table by row
  df <- t(apply(df, 1, function(x) x / sum(x)))
  df <- as.data.frame(df)
  if(nrow(df)==1){
    max_row <- rep(rownames(df),ncol(df))
  }else{
  # Divide the dataframe by the column sums
  df <- apply(df, 2, function(x) x / sum(x))
  df <- as.data.frame(df)
  # Get row names of maximum values in each column
  max_row <- apply(df, 2, function(x) rownames(df)[which.max(x)])

  # Calculate threshold
  threshold <- 1 / nrow(df) + threshold_add

  # Get max_boolean and replace row names
  max_boolean <- apply(df, 2, function(x) max(x) > threshold)
  max_row[!max_boolean] <- "."

  return(max_row)
  }
}

#== this function return max group in given meta---------------------------

anno_dict <- function(object,cluster_id, group_select,threshold_add=0.2){
  row_name <- group_select
  col_name <- table(object@meta.data[[cluster_id]])%>%as.matrix()%>%rownames()
  annotree <-  data.frame(matrix("", nrow=length(row_name), ncol=length(col_name)))

  rownames(annotree) <- row_name
  colnames(annotree) <- col_name
  for (i in row_name){
    annotree[i,] <- anno_row(object,cluster_id,i,threshold_add)
  }
  return(annotree)
}

#== this function combine every given cluster levels
anno_bind <- function(object,cluster_select, group_select,threshold_add=0.2){
  annolist <- lapply(cluster_select,function(x) anno_dict(object,cluster_id = x, group_select=group_select,threshold_add = threshold_add)%>%t)
  annoBind <- do.call(rbind,annolist)
  annoBind%>%as.data.frame()
}

#== load json parameters---------
load_json <- function(jsonPath) {
  parameter_list <- jsonlite::read_json(jsonPath)
  parameter_list <- lapply(parameter_list, function(x) {
    if (is.list(x)) {
      return(unlist(x))
    } else {
      return(x)
    }
  })

  # Assign parameter_list to the global environment
  assign("parameter_list", parameter_list, envir = .GlobalEnv)

  # List of specific parameters to load into global environment
  params_to_load <- c("n_cores_markers", "assay_markers", "assay_slot", "logfc.threshold",
                      "min.pct", "min.diff.pct", "max.cells.per.ident", "min.cells.feature",
                      "min.cells.group", "base", "only.pos", "batch_var", "test.use",
                      "specificity_base", "add_batch_as_latent", "additional_suffix")

  # Assign specified parameters to the global environment
  for (param_name in params_to_load) {
    if (param_name %in% names(parameter_list)) {
      assign(param_name, parameter_list[[param_name]], envir = .GlobalEnv)
    }
  }

  # Print a message about test.use
  message("test.use ", get("test.use", envir = .GlobalEnv), " ... ", parameter_list$test.use)
}

