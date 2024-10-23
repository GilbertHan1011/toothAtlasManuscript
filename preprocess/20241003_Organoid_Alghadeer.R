fileName <- list.files("../202409_tooth_raw/Organoid_Alghadeer/",pattern = "GSM5*")
extracted <- gsub(".*?_(\\w+)_(\\w+)_\\w+_\\w+.*", "\\1_\\2", fileName) %>% unique
extracted <- extracted[1:15]

fileList <- lapply(extracted,function(x) list.files(path = "../202409_tooth_raw/Organoid_Alghadeer/",pattern = paste0(x)))
matrix <- lapply(fileList,`[`,3) %>% unlist
cellAnno <- lapply(fileList,`[`,1)%>% unlist
geneAnno <- lapply(fileList,`[`,2)%>% unlist
geneFile <- read.table(paste0("../202409_tooth_raw/Organoid_Alghadeer/",geneAnno[[1]]),header = F)
cellFile <- read.table(paste0("../202409_tooth_raw/Organoid_Alghadeer/",cellAnno[[1]]),header = F)
matFile <- read.table(paste0("../202409_tooth_raw/Organoid_Alghadeer/",matrix[[1]]))
count_matrix <- Matrix::sparseMatrix(

  i = matFile$V1,

  j = matFile$V2,

  x = matFile$V2,

  dims = c(nrow(geneFile), nrow(cellFile)),

  dimnames = list(geneFile$V1, cellFile$V1)

)
seurat <- CreateSeuratObject(count_matrix,project = extracted[[1]],min.cells = 3,min.features = 300)
loadFun <- function(i){
  geneFile <- read.table(paste0("../202409_tooth_raw/Organoid_Alghadeer/",geneAnno[[i]]),header = F)
  cellFile <- read.table(paste0("../202409_tooth_raw/Organoid_Alghadeer/",cellAnno[[i]]),header = F)
  matFile <- read.table(paste0("../202409_tooth_raw/Organoid_Alghadeer/",matrix[[i]]))
  count_matrix <- Matrix::sparseMatrix(

    i = matFile$V1,

    j = matFile$V2,

    x = matFile$V2,

    dims = c(nrow(geneFile), nrow(cellFile)),

    dimnames = list(geneFile$V1, cellFile$V1)
  )
  seurat <- CreateSeuratObject(count_matrix,project = extracted[[i]],min.cells = 3,min.features = 300)
  seurat$orig.ident <- extracted[[i]]
  seurat
}
seuratList <- lapply(1:15, loadFun)
seuratList <- lapply(seuratList, qcFun,Species = "Human")
seuratMerge <- merge(seuratList[[1]],seuratList[2:15])
seuratMerge <- runharmony(seuratMerge)


ens <- geneFile[[1]]
symbol <- geneFile[[2]]
names(symbol) <- ens

symbol <- symbol[rownames(seuratMerge)]
rownames(seuratMerge) <- symbol
DimPlot(seuratMerge)
FeaturePlot(seuratMerge,"SMPD3")
FeaturePlot(seuratMerge,"AMEL")
FeaturePlot(seuratMerge,"DSP")
FeaturePlot(seuratMerge,"VCAN")
seuratMerge$nFeature_RNA %>% hist(breaks = 100)
seuratMerge@assays$RNA@layers$scale.data <- NULL

saveRDS(seuratMerge,"preprocess_data/Organoid_Alghadeer.Rds")
