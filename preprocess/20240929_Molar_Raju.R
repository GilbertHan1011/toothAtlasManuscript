rm(list=ls())
library(Seurat)
library(dplyr)
library(purrr)
source("script/utils/seurat_utils.R")
fileNames <- list.files("../202409_tooth_raw/Molar_Raju/",full.names = T)
fileShort <- list.files("../202409_tooth_raw/Molar_Raju/",full.names = F)

countList <- lapply(fileNames,Read10X)
seuratList <- map2(countList,fileShort,function(x,y) CreateSeuratObject(x,min.cells = 3, min.features = 500,project = y))
seuratList <- lapply(seuratList,qcFun)
listLength <- length(seuratList)
seuratMerge <- merge(seuratList[[1]],seuratList[2:listLength])
seuratMerge <- runharmony(seuratMerge)
# Assuming `qcFun` and `runharmony` are defined in your sourced utils file
source("script/utils/seurat_utils.R")

# to make it more easy, I make a function
process_seurat_data <- function(data_directory, min_cells = 3, min_features = 500, qc_mito = 20,qc_rna = 300,species = "Mouse") {
  require(Seurat)
  require(dplyr)
  require(purrr)
  # List all files in the specified directory
  fileNames <- list.files(data_directory, full.names = TRUE)
  fileShort <- list.files(data_directory, full.names = FALSE)

  # Read the 10X data
  countList <- lapply(fileNames, Read10X)

  # Create Seurat objects with quality control
  seuratList <- map2(countList, fileShort, function(x, y) {
    CreateSeuratObject(x, min.cells = min_cells, min.features = min_features, project = y)
  })

  # Apply quality control function
  seuratList <- lapply(seuratList, qcFun,percent_mito=qc_mito,nFeature_RNA=qc_rna,Species=species)

  # Merge all Seurat objects
  seuratMerge <- merge(seuratList[[1]], seuratList[2:length(seuratList)])

  # Run Harmony
  seuratMerge <- runharmony(seuratMerge)

  return(seuratMerge)
}

markers <-
  c("Sox9",
    "Msx2",
    "Vcan",
    "Kit",
    "Pclo",
    "Igfbp5",
    "Postn",
    "Fst",
    "Igfbp2",
    "Smoc2",
    "Sfrp2",
    "Igfbp3",
    "Smpd3",
    "Gsc",
    "Wnt10a",
    "Cdk1",
    "Pitx2",
    "Krt14",
    "Dsp",
    "C1qa",
    "Aif1",
    "Napsa",
    "Cdh5",
    "Rgs5",
    "Acp5"
  )

FeaturePlot(seuratMerge,markers[1:9])
FeaturePlot(seuratMerge,markers[10:18])
FeaturePlot(seuratMerge,markers[1:9])
FeaturePlot(seuratMerge,markers[10:18])
saveRDS(seuratMerge,"preprocess_data/Molar_Raju.Rds")
