rm(list = ls())
library(tidyverse)
setwd("script/utils/scib/")
source('plotBestMethodsRNA.R') #== code from scib-reproducibility
source('plotSingleTaskRNA.R') #== code from scib-reproducibility
metricTest <- read.csv("../202409_tooth_raw/scib/metrics.csv")%>%column_to_rownames("X")
#setwd("script/scib_script/visualization/")

dir.create("../../../process/pre-intergration/scib_metrics/")
plotSingleTaskRNA(
  csv_metrics_path = "../../../../202409_tooth_raw/scib/metrics.csv",
  outdir = "../../../process/pre-intergration/scib_metrics/",
  weight_batch = 0.4
)

