devtools::install_github("aet21/SCENT")
#library(clusterProfiler)
source("script/utils/seurat_utils.R")
#== try my own data------------------------------
mes.m <- as.matrix(mes@assays$originalexp@data)
load("data/geneinfo_2022.rda")
convertGene <- function(symbols) {
  mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])

  humansymbol2mousesymbol = mapper(geneinfo_2022, "entrez", "symbol_mouse")
  converted_symbols = humansymbol2mousesymbol[symbols %>% as.character()]

  return(converted_symbols)
}
#entriz = convertGene(rownames(mes.m))
mes_pheno <- mes$C9_named

geneTb2 <- convertGene(rownames(mes.m))%>%na.omit()

mes.m <- mes.m[names(geneTb2),]
rownames(mes.m) <- geneTb2[rownames(mes.m)]


#rownames(mes.m) <- geneTb$ENTREZID
ccat.mes <- CompCCAT(exp = mes.m, ppiA = net17Jan16.m);
boxplot(ccat.mes ~ mes_pheno, main = "SR potency estimates",
        xlab = "Cell Type", ylab = "SR")
mes$scent <- ccat.mes
FeaturePlot(mes,"scent")+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),
                                                values = c(0,0.4,0.55,0.65,1.0))
ggsave("results/trajectory/20241124_multipotent/20241124_scent_score_preknn.pdf",width = 6,height = 6)
scentMes <- as.data.frame(mes$scent)
write.csv(scentMes,"process/trajectory/20241124_multipotent/20241124_scent_result.csv")
