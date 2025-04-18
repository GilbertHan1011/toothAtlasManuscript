# Load required packages
library(igraph)
library(RColorBrewer)
library(dplyr)

network <- read.csv("~/Desktop/disk2/202409_tooth/process/trajectory/20250417_mineralization_downstream/output_cytoscape.csv")
rownames(mineFitHM)[order(apply(mineFitHM, 1, which.max))]
df = data.frame(gene = rownames(mineFitHM_ordered),geneorder = 1:length(rownames(mineFitHM_ordered)))
rownames(df) <- rownames(mineFitHM_ordered)
df= df[rownames(mineFitHM),]
#geneOrder = data.frame(gene = rownames(mineFitHM),order = match(rownames(mineFitHM_ordered), rownames(mineFitHM)))
geneOrder = data.frame(
  gene = rownames(mineFitHM),
  order = match(rownames(mineFitHM), rownames(mineFitHM_ordered))
)
colnames(geneOrder) <- c("Target","or")
network <- left_join(network,geneOrder)
write.csv(network,"~/Desktop/disk2/202409_tooth/process/trajectory/20250417_mineralization_downstream/output_cytoscape_mod.csv",quote = F)
write.csv(geneOrder, "~/Desktop/disk2/202409_tooth/process/trajectory/20250417_mineralization_downstream/20250418_hm_geneorder.csv")

# Function to plot network from file
plot_tf_gene_network <- function(network, output_file = "network_plot.pdf",
                                 width = 10, height = 8,num_gene = 3) {

  # Create graph object
  g <- graph_from_data_frame(d = network[, c("Source", "Target", "Fold_Enrichment", "P_value")],
                             directed = TRUE)

  # Get unique nodes and their types (TF or gene)
  nodes <- V(g)
  node_type <- ifelse(names(nodes) %in% unique(network$Source), "TF", "Gene")

  # Get node order values (for genes only, TFs will have NA)
  node_order <- rep(NA, length(nodes))

  # For each gene, find its order from the network data
  for (i in 1:length(nodes)) {
    if (node_type[i] == "Gene") {
      gene_name <- names(nodes)[i]
      # Get the order value from the first occurrence of this gene
      node_order[i] <- network$order[network$Target == gene_name][1]
    }
  }

  # Find genes co-regulated by at least three TFs
  gene_counts <- table(network$Target)
  coregulated_genes <- names(gene_counts)[gene_counts >= num_gene]

  # Set node colors based on order (use a color gradient)
  node_color <- rep("red", length(nodes))  # Default color for TFs

  # For genes, use a color palette based on order
  gene_indices <- which(node_type == "Gene")
  if (length(gene_indices) > 0) {
    # Create color palette based on range of order values
    order_values <- node_order[gene_indices]
    order_values <- order_values[!is.na(order_values)]

    if (length(order_values) > 0) {
      # Get a color palette
      n_colors <- 100
      color_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n_colors)

      # Map order values to color indices
      min_order <- min(order_values, na.rm = TRUE)
      max_order <- max(order_values, na.rm = TRUE)

      for (i in gene_indices) {
        if (!is.na(node_order[i])) {
          # Scale the order value to an index in the color palette
          color_idx <- floor((node_order[i] - min_order) / (max_order - min_order) * (n_colors - 1)) + 1
          node_color[i] <- color_palette[color_idx]
        }
      }
    }
  }

  # Set node sizes (larger for TFs and co-regulated genes)
  node_size <- rep(8, length(nodes))  # Default size for genes
  node_size[node_type == "TF"] <- 15  # Larger size for TFs

  # Increase size for co-regulated genes
  for (i in 1:length(nodes)) {
    if (node_type[i] == "Gene" && names(nodes)[i] %in% coregulated_genes) {
      node_size[i] <- 12  # Intermediate size for co-regulated genes
    }
  }

  # Create node labels (TF names and co-regulated gene names)
  node_labels <- rep(NA, length(nodes))
  node_labels[node_type == "TF"] <- names(nodes)[node_type == "TF"]

  # Add labels for co-regulated genes
  for (i in 1:length(nodes)) {
    if (node_type[i] == "Gene" && names(nodes)[i] %in% coregulated_genes) {
      node_labels[i] <- names(nodes)[i]
    }
  }

  # Set edge widths based on Fold_Enrichment
  edge_width <- E(g)$Fold_Enrichment / max(E(g)$Fold_Enrichment) * 5

  # Set edge colors based on p-value
  edge_color <- "gray40"

  # Set layout (try different layouts to see which works best)
  layout <- layout_with_fr(g)

  # Open PDF file for output
  pdf(output_file, width = width, height = height)

  # Plot the network
  plot(g,
       vertex.color = node_color,
       vertex.size = node_size,
       vertex.label = node_labels,
       vertex.label.color = "black",
       vertex.label.dist = 0,
       vertex.label.family = "sans",
       vertex.label.cex = 0.8,
       edge.width = edge_width,
       edge.color = edge_color,
       edge.arrow.size = 0.5,
       layout = layout,
       main = "TF-Gene Network")

  # Add legend for node colors (order)
  if (length(unique(node_order[!is.na(node_order)])) > 1) {
    legend_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(5)
    legend_labels <- round(seq(min_order, max_order, length.out = 5))

    legend("bottomright",
           legend = c(legend_labels, "Co-regulated genes (≥3 TFs)"),
           col = c(legend_colors, "black"),
           pch = c(rep(19, 5), 1),
           pt.cex = c(rep(1, 5), 1.5),
           title = "Legend",
           cex = 0.8)
  }

  # Close the PDF device
  dev.off()

  # Print summary of co-regulated genes
  if (length(coregulated_genes) > 0) {
    cat("Genes co-regulated by at least 3 TFs:", length(coregulated_genes), "\n")
    cat("Co-regulated genes:", paste(coregulated_genes, collapse=", "), "\n")
  } else {
    cat("No genes are co-regulated by at least 3 TFs\n")
  }

  cat("Network plot saved to:", output_file, "\n")

  # Return information about co-regulated genes for further analysis
  return(list(
    coregulated_genes = coregulated_genes,
    regulation_counts = gene_counts
  ))
}

plot_tf_gene_network(network = network,output_file = "results/trajectory/20250415_trajdtw_fit/20250418_network.pdf")
plot_tf_gene_network(network = network,output_file = "results/trajectory/20250415_trajdtw_fit/20250418_network_2.pdf",num_gene = 2)

plot_tf_gene_network(network = network,output_file = "results/trajectory/20250415_trajdtw_fit/20250418_network_3.pdf",num_gene = 1)
