# devtools::install_github('theislab/kBET')
# devtools::install_local("../../../soft/kBET-master.zip")
library(kBET)
# environment-----------------

rna = zellkonverter::readH5AD("process/pre-intergration/big_data/20241019_lop10_hvg.h5ad")

outdir <- "results/preIntergration/20241020_kbet/"
cellName <- read.csv("processed_data/attributeName/cellName.csv",row.names = 1) %>% unlist
varName <- read.csv("processed_data/attributeName/varName.csv",row.names = 1) %>% unlist
#dir.create(outdir)

rownames(rna) <- varName
colnames(rna) <- cellName
rna_seurat <- as.Seurat(rna,counts = "counts", data = "logcounts")
saveRDS(rna_seurat,"process/pre-intergration/big_data/20241020_seurat.Rds")

for (i in unique(rna_seurat$Project)){
  study=i
  target=rna_seurat[,rna_seurat$Project==study]
  data = GetAssayData(object = target, slot = "counts",assay = "originalexp")
  data <- t(data)
  batch <- target$Sample
  subset_size <- 1000/length(batch) #subsample to 10% of the data

  if (length(batch)>2000){
    subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
    batch.estimate <- kBET(data[subset_id,], batch[subset_id], plot=TRUE, do.pca = TRUE, dim.pca = 2)
  }
  else{
    batch.estimate <- kBET(data, batch, plot=TRUE, do.pca = TRUE, dim.pca = 15)
  }
  saveRDS(batch.estimate,paste0(outdir,study,"_estimateResult.Rds"))
  plot.data <- data.frame(class=rep(c('observed', 'expected'),
                                    each=length(batch.estimate$stats$kBET.observed)),
                          data =  c(batch.estimate$stats$kBET.observed,
                                    batch.estimate$stats$kBET.expected))
  g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() +
    labs(x='Test', y='Rejection rate',title='kBET test results') +
    theme_bw() +
    scale_y_continuous(limits=c(0,1))
  ggsave(paste0(outdir,study,"_reject_box.pdf"))
}


kbetFile <- list.files(path = outdir,pattern = "*.Rds",full.names = T)
kbetName <- gsub(outdir,"",kbetFile) %>%  gsub("\\/","",.)
kbetName <- gsub("_estimateResult.Rds","",kbetName)
kbetList <- lapply(kbetFile,readRDS)
names(kbetList) <- kbetName
kbet <- lapply(kbetList,function(x) x$summary$kBET.observed)
kbetDf <- do.call(cbind,kbet)
boxplot(kbetDf)
ggplot2::geom_boxplot(kbetDf)
library(reshape2)

# Specify id.vars: the variables to keep but not split apart on
kbetlong=melt(kbetDf)[,-1]
ggplot(kbetlong,aes(x=Var2,y=value))+geom_boxplot()+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = "bold"))+
  xlab("Project")+ylab("rejection rate")
ggsave(paste0(outdir,"_overall.pdf"),width = 10,height = 8)

# # create batch metadata---------------------------------------------------
createSampleMeta <- function(seurat,meta=osteoMeta,name){
  Idents(seurat) <- seurat$Sample
  new.id <- meta[name]%>%unlist()
  names(new.id) <- meta["Sample"]%>%unlist()
  seurat <- RenameIdents(seurat,new.id)
  seurat@meta.data[name] <- Idents(seurat)
  Idents(seurat) <- seurat$Sample
  return(seurat)
}
coreData <- read.csv("data/metadata//metadata_latest.csv")
rna_seurat <- createSampleMeta(rna_seurat,coreData,"batch")
rna_seurat@meta.data[c("Sample","batch")] %>% unique %>% as.data.frame() %>% View
rna_seurat@meta.data[c("Sample","batch")] %>% View
meta <- rna_seurat@meta.data
write.csv(meta,"processed_data/metadata/20241020_metadata.csv")
write.csv(meta$batch,"processed_data/attributeName/batch.csv")
# osteoMerge$short_id <- as.character(osteoMerge$short_id)
# osteoMerge$short_id[osteoMerge$Project=="CalvariaP4_Ayturk"] <- "CAy"
# osteoMerge$short_id[osteoMerge$Project=="CranioSoxc_Angelozzi"] <- "CAn"
# bm <- RenameCells(osteoMerge[,osteoMerge$Project=="BMSC-Specification_Kishor"],add.cell.id="BMSCSpecification_Kishor_")
# View(colnames(osteoMerge)%>%as.data.frame())
# colName <- colnames(osteoMerge[,osteoMerge$Project!="BMSC-Specification_Kishor"])
# colName <- c(colName,colnames(bm))
# test <- RenameCells(osteoMerge,new.names=colName)
# osteoMerge <- test
# rm(test)
# write.csv(osteoMerge@meta.data,"../data/annodata/1.18_osteoMeta.csv")
# SaveH5Seurat(osteoMerge,"../data/annodata/1.18_osteoMerge.h5Seurat")
# Convert("../data/annodata/1.18_osteoMerge.h5Seurat",dest="h5ad")
#
# saveRDS(osteoMerge,"../data/annodata/1.18_osteo_merge.Rds")

