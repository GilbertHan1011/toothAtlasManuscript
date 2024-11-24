rm(list = ls())

source("script/utils/seurat_utils.R")
library(destiny)
library(tidyverse)

cellName = read.csv("processed_data/framework/attributeName/cellName_mes_keep20241119.csv",row.names = 1) %>% unlist()
mes <- readRDS("processed_data/integrated_data/20241106_mesenchyme.Rds")
mes <- mes[,cellName]
find_sigmas(mes@reductions$X_SCANVI@cell.embeddings)

downstream_dm <- run_diffMap(t(mes@reductions$X_SCANVI@cell.embeddings),
                             mes$C9_named, sigma = 2)


plot_eigenVal(dm = downstream_dm)
my_color<-brewer.pal(6,'Spectral')
splom(~downstream_dm@eigenvectors[, 1:6], groups = mes$C9_named, col = my_color,
      key = list(space="right", points = list(pch = 19, col = my_color),
                 text = list(c(unique(mes$C9_named)))))

dm_vector <- downstream_dm@eigenvectors
write.csv(dm_vector,"processed_data/framework/embedding/20241119_destiny_diffmap.csv")
