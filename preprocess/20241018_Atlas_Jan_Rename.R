#setwd("../../../")
seurat <- readRDS("processed_data/preprocess_data/Atlas_Jan_Mouse.Rds")
source("script/utils/seurat_utils.R")
geneHuman <- rownames(seurat)
geneMouse <- convertGene("data/geneinfo_2022.rda",geneHuman)
seurat <- seurat[!is.na(geneMouse)]
geneMouse <- geneMouse[!is.na(geneMouse)]
rownames(seurat) <- geneMouse

RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]              <- newnames
    if (length(RNA@scale.data)) rownames(RNA@scale.data)    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

seurat <- RenameGenesSeurat(seurat,geneMouse)
saveRDS(seurat,"processed_data/preprocess_data/Atlas_Jan_Mouse_renamed.Rds")
