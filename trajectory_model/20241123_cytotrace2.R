library(tidyverse)
library(CytoTRACE2)
library(RColorBrewer)
devtools::install("~/soft/cytotrace2/cytotrace2_r/")

devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
mes <- readRDS("processed_data/integrated_data/20241106_mesenchyme.Rds")
cellName <- read.csv("processed_data/framework/attributeName/cellName_mes_keep20241119.csv",row.names = 1) %>% unlist()
mes <- mes[,cellName]
loadData_fromSeurat <- function (object, slot_type,assay)
{
  data <- as.data.frame(Seurat::GetAssayData(object = object,
                                             assay = "originalexp", slot = slot_type))
  return(data)
}

cytotrace2_result <- cytotrace2(mes, is_seurat = TRUE, slot_type = "data", species = 'mouse',assay = "originalexp")

test <- loadData_fromSeurat(mes,slot_type = "data",assay = "originalexp")

cytotrace2_result <- cytotrace2(test)
mes$cytotrace2 <- cytotrace2_result$CytoTRACE2_Score

p1 <- FeaturePlot(mes,"cytotrace2")+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),
                                                     values = c(0,0.5,0.7,0.85,1.0))
ggsave("results/trajectory/20241124_multipotent/20241124_cytotrace2_score.pdf",width = 6,height = 6)

mes$cytotrace2_preKNN <- cytotrace2_result$preKNN_CytoTRACE2_Score
p2 <- FeaturePlot(mes,"cytotrace2_preKNN")+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),
                                                     values = c(0,0.5,0.7,0.85,1.0))
ggsave("results/trajectory/20241124_multipotent/20241124_cytotrace2_score_preknn.pdf",width = 6,height = 6)
write.csv(cytotrace2_result,"process/trajectory/20241124_multipotent/20241124_cytotrace2_result.csv")

p1|p2
ggsave("results/trajectory/20241124_multipotent/20241124_cytotrace2_align.pdf",width = 8,height = 4)
