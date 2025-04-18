odonto = readRDS("/home/gilberthan/Desktop/disk2/202409_tooth/processed_data/integrated_data/20250326_mesenchyme.Rds")
odonto$odonto <- "other"
odonto$odonto[odonto$C9_named%in%c("Odontoblast")] = "odonto"
odonto$odonto[odonto$C9_named%in%c("Apical papilla","Dental mes")] = "mes"
Idents(odonto) <- odonto$odonto
odontoMarker <- FindMarkers(odonto,ident.1 = "odonto",ident.2 = "mes")
write.csv(odontoMarker,"process/trajectory/20250415_GO_downstream/20250416_odonto_arker_seurat.csv")
odontoMarker <- odontoMarker%>%rownames_to_column("gene")%>%left_join(conserve_odonto)



line_thred = 21
left_thred = 3
right_thred = 8
# First, let's calculate the counts of points in each region
# Filter points inside the box (x from 3 to 6, y from 0 to 100)
in_box <- odontoMarker %>%
  filter(avg_log2FC >= left_thred & avg_log2FC <= right_thred &
           -raw_score >= 0 & -raw_score <= 150)


# Count points above and below y = 20 within the box
above_line <- sum(-in_box$raw_score > line_thred)
below_line <- sum(-in_box$raw_score <= line_thred)

# Create a data frame for the count labels
count_labels <- data.frame(
  x = c(4.5, 4.5),
  y = c(60, 10),
  label = c(paste0("Above: ", above_line), paste0("Below: ", below_line))
)

# Now create the plot with all elements
ggplot(odontoMarker, aes(x = avg_log2FC, y = -raw_score)) +
  # Original points
  geom_point(alpha = 0.7, size = 2, color = "DeepSkyBlue3") +

  # Add red box from x=(3,6) and y=(0,100)
  geom_rect(aes(xmin = left_thred, xmax = right_thred, ymin = 0, ymax = 150),
            fill = NA, color = "red", linewidth = 1) +

  # Add red line at y = 20
  geom_hline(yintercept = 20, color = "red", linewidth = 1) +

  # Original reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +

  # Add count labels
  geom_text(data = count_labels, aes(x = x, y = y, label = label),
            color = "red", fontface = "bold", size = 5) +

  # Rest of the original styling
  labs(
    title = "Relationship between DE and Conservation",
    subtitle = paste("Points in box:", nrow(in_box),
                     "| Above y=20:", above_line,
                     "| Below y=20:", below_line),
    x = "Average log2 Fold Change",
    y = "Conservation Score (-raw_score)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "darkred"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  # Set axis limits to ensure the entire box is visible
  coord_cartesian(xlim = c(min(de_gene_join_run4$avg_log2FC, 0),
                           max(de_gene_join_run4$avg_log2FC, 7)),
                  ylim = c(min(-de_gene_join_run4$raw_score, 0),
                           max(-de_gene_join_run4$raw_score, 105)))
ggsave("results/trajectory/20250415_trajdtw_fit/20250416_odonto_benchmarking.pdf",width = 8,height = 6)

gene_above = in_box$gene[-in_box$raw_score > line_thred]
gene_below = in_box$gene[-in_box$raw_score <= line_thred]

bm_list <- list(conservate = gene_below, unconservate = gene_above)
bm_ck_res3 <- clusterProfiler::compareCluster(bm_list, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(bm_ck_res3)
ggsave("results/trajectory/20250415_trajdtw_fit/20250416_odonto_benchmarking_GO_enrichment.pdf",width = 6,height = 8)


tradeseq_odonto <- read.csv("process/trajectory/tradeseq/20250403_gene_cluster.csv",row.names = 1)
tradeseqgene_1 <- tradeseq_odonto$gene[tradeseq_odonto$sigGene_ordered%in%c("GC4","GC5","GC8")]
conservate_odonto_trade <- conserve_odonto[conserve_odonto$gene%in%tradeseqgene_1,]
conservate_odonto_trade <- na.omit(conservate_odonto_trade)
odonto_trade_conserve <- conservate_odonto_trade$gene%>%head(50)
odonto_trade_unconserve <- conservate_odonto_trade$gene%>%tail(50)

tradeseq_list <- list(conservate = odonto_trade_conserve, unconservate = odonto_trade_unconserve)
bm_ck_res4 <- clusterProfiler::compareCluster(tradeseq_list, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(bm_ck_res4,showCategory=10)
ggsave("results/trajectory/20250415_trajdtw_fit/20250416_odonto_tradeseq_GO_enrichment.pdf",width = 6,height = 10)
