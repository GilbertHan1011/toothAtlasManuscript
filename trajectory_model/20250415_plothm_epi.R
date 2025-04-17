epiHM <- read.csv("process/trajectory/20250414_epi_run_2/fitted_trajectories_optimized.csv",row.names = 1)
epiLabel <- read.csv("process/trajectory/20250414_epi_run_2/epi_raw_gene.csv",row.names = 1)
epiName <- rownames(epiLabel)[rownames(epiLabel)%in%rownames(epiHM)]
epiLabel <- epiLabel[epiName,]
epiHM <- t(epiHM)
epiHM <- epiHM[epiName,]
pdf("results/trajectory/20250415_trajdtw_fit/20250415_epi_fit_hm.pdf")
Heatmap(epiHM,cluster_rows = F,cluster_columns = F,row_split = epiLabel,show_row_names = F,show_column_names = F,col=colorRamp2(c(-1, 0, 1), c("DeepSkyBlue3", "white", "red")) )
dev.off()
