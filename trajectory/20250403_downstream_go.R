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
