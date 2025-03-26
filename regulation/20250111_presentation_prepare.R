## Presentation prepare

mes_origin <- readRDS("processed_data/integrated_data/20241024_mesenchyme.Rds")
mes_origin <- NormalizeData(mes_origin)
FeaturePlot(mes_origin,c("Sp7","Ambn"))


#== enrichment-----------
library(clusterProfiler)
library(org.Mm.eg.db)
library(simplifyEnrichment)
geneList1 = candidateList$gene[candidateList$prior == 1]
goRes <- enrichGO(gene         = geneList1,
         OrgDb         = org.Mm.eg.db,
         keyType       = 'SYMBOL',
         ont           = "BP",
         pAdjustMethod = "BH",
         pvalueCutoff  = 0.01,
         qvalueCutoff  = 0.05)
mat = GO_similarity(goRes@result$ID)
simplifyGO(mat)


#== Ovol2 centric-----------
ovolDf = readxl::read_xls("process/regulon/Ovol2_centric/Ovol2-TFs network.xls")

genes <- ovolDf$To
intersect(genes,geneList1)
intersect(genes,candidateList$gene)


goRes2 <- enrichGO(gene         = genes,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
dotplot(goRes2)

ggsave("results/regulation/Ovol2_centric/20250111_ovol2_regulon_enrichment.pdf",width = 6,height = 6)

View(goRes2@result)
