# Trajectory Analysis Downstream Enrichment

This document demonstrates the downstream gene ontology enrichment analysis of gene clusters identified from trajectory analysis. The purpose of this analysis is to identify biological processes, molecular functions, and cellular components associated with each gene cluster along the developmental trajectory.

## Workflow Overview

1. **Data Input**: The analysis starts by loading gene clusters obtained from tradeSeq trajectory analysis.
2. **GO Enrichment**: Using clusterProfiler, gene ontology enrichment analysis is performed for each gene cluster.
3. **Representative GO Terms Selection**: A custom function `make_cluster_comparison` selects the most representative GO terms for each cluster by:
   - Converting p-values to -log10 scale
   - Comparing each cluster's enrichment against the average of other clusters
   - Selecting the top terms per cluster
4. **Visualization**: Results are displayed in a dotplot with a custom color gradient, where:
   - Dot size represents gene count
   - Color intensity represents statistical significance (-log10 p-value)
   - Red indicates high significance
   - Grey indicates low significance
   - Values exceeding the scale limit are displayed in red

The dotplot visualization provides an intuitive way to understand the biological functions associated with different gene clusters in the developmental trajectory.

```R
library(clusterProfiler)
genelist <- read.csv("process/trajectory/tradeseq/20250403_gene_cluster.csv",row.names = 1)
trajList <- split(genelist$gene,genelist$sigGene_ordered)
ck <- compareCluster(geneCluster = trajList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")


make_cluster_comparison <- function(goCompare,goPerGroup = 5) {
  ckRes <- goCompare@compareClusterResult
  
  ckResSub <- ckRes[c("Description", "Cluster", "p.adjust")]
  
  ckWide <- pivot_wider(ckResSub, names_from = Cluster, values_from = p.adjust)
  
  ckWide <- ckWide %>% column_to_rownames("Description")
  
  ckWideMod <- ckWide %>%
    log10() * -1
  
  ckWideMod[is.na(ckWideMod)] <- 0
  
  ckRowList <- list()
  
  for (i in 1:dim(ckWideMod)[2]) {
    ckRow <- ckWideMod[i] - rowMeans(ckWideMod[-i])
    ckRowList[[i]] <- ckRow
  }
  
  ckRowDf <- data.frame(do.call(cbind, ckRowList))
  
  ckSelect <- c()
  
  for (i in 1:dim(ckWideMod)[2]) {
    ckSelect <- c(ckSelect, ckRowDf %>% dplyr::arrange(desc(across(i))) %>% rownames %>% .[1:goPerGroup])
  }
  ckSelect <- ckSelect %>% unique
  ckResSelect <- ckRes %>% dplyr::filter(Description %in% ckSelect)
  
  ckResSelect$p.adjust <- log10(ckResSelect$p.adjust)
  
  ckResSelect$Description <- factor(ckResSelect$Description, levels = ckSelect)
  
  ckResSelect <- dplyr::arrange(ckResSelect, Description)
  
  goCompare@compareClusterResult <- ckResSelect
  
  return(goCompare)
}
res_select <- make_cluster_comparison(ck,goPerGroup = 5)

dotplot(res_select)+ scale_fill_gradient(low = "red",high = "grey90",  limits = c(-50, 0)) 
dotplot(res_select) +
  scale_fill_gradientn(colors = c("red", "gray90"), 
                       limits = c(-20, 0), 
                       na.value = "red") +  # Set NA values to max color (grey90)
  labs(fill = "-log10(p.adjust)") +
  theme_minimal()
ggsave("results/final_plot/20250403_traj_dotplot.pdf",width = 8,height = 6)
```