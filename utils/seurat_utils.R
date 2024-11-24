# change orig ident-----------------------
renameFun <- function(x,id){
  Idents(x) <- x$orig.ident
  x <- RenameIdents(x,id)
  x$orig.ident <- Idents(x)
  return(x)
}

# Annotate Label
# first parameter : seurat; Second parameter : new ID, which order same with level; third para : the slot you want to store in metadata
renameLabel <- function(x,id,label){
  names(id) <- levels(x)
  x <- RenameIdents(x,id)
  x@meta.data[[label]] <- Idents(x)
  return(x)
}


#== create sample meta---------------------------

createSampleMeta <- function(seurat,meta=osteoMeta,name){
  Idents(seurat) <- seurat$Sample
  new.id <- meta[name]%>%unlist()
  names(new.id) <- meta["Sample"]%>%unlist()
  seurat <- RenameIdents(seurat,new.id)
  seurat@meta.data[name] <- Idents(seurat)
  Idents(seurat) <- seurat$Sample
  return(seurat)
}
#== seurat quality control----------------
qcFun <-  function(x,nCount_RNA=800,percent_mito=20,nFeature_RNA=300,Species="Mouse"){
  if(Species == "Mouse"){
    x <- PercentageFeatureSet(x, "^mt-", col.name = "percent_mito")
  } else if(Species == "Human"){
    x <- PercentageFeatureSet(x, "^MT-", col.name = "percent_mito")
  }
  x <- x[,(x$nCount_RNA > nCount_RNA & x$percent_mito < 20 & x$nFeature_RNA > 300)]
  #x <- subset(x, cells = selected_count)
  return(x)
}

#== dbl funciton-------------
dbFinderFun <- function(x){
  require(scDblFinder)
  require(BiocParallel)
  print(unique(x$orig.ident))
  DefaultAssay(x) <- "RNA"
  sce <- as.SingleCellExperiment(x)
  sce <- scDblFinder(sce, samples="orig.ident", BPPARAM=MulticoreParam(10))
  droplet_class = sce$scDblFinder.class
  x$scDblFinder_class <- droplet_class
  return(x)
}

#== run seurat pipeline-------------------
runSeurat <- function(x,dim=30,seed = 0){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  x <- FindNeighbors(x, dims = 1:dim)
  x <- FindClusters(x, resolution = 0.5,  random.seed= seed)
  x <- RunUMAP(x, dims = 1:dim,random.seed = seed)
}

# seurat cca merge function-----------------------------------
ccaMergeFun <- function(seuratList){
  seuratList <- lapply(X = seuratList, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  features <- SelectIntegrationFeatures(object.list = seuratList)
  seurat.anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features)
  seurat.combined <- IntegrateData(anchorset = seurat.anchors)
  DefaultAssay(seurat.combined) <- "integrated"
  seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
  seurat.combined <- RunPCA(seurat.combined, npcs = 50, verbose = FALSE)
  seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:30)
  seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:30)
  seurat.combined <- FindClusters(seurat.combined, resolution = 0.5)
  return(seurat.combined)
}

#== harmony function---------------------------------------
runharmony <- function(x,dim=30){
  require(harmony)
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  x <- RunHarmony(x, "orig.ident",reduction.use = "pca")
  x <- FindNeighbors(x, dims = 1:dim,reduction = "harmony")
  x <- FindClusters(x, resolution = 0.5)
  x <- RunUMAP(x, dims = 1:dim,reduction = "harmony")

}


#== use for large dataset to matrix--------------------
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])

  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x

  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)

}


#== trajectory function-------------------
run_diffMap <- function(data=data, condition=condition, sigma="local", k = 20){
  destinyObj <- as.ExpressionSet(as.data.frame(t(data)))
  destinyObj$condition <- factor(condition)
  dm <- DiffusionMap(destinyObj, sigma, k)
  return(dm)
}

plot_eigenVal <- function(dm=dm){
  linepad <- .5
  plot(
    eigenvalues(dm),
    ylim = 0:1,
    pch = 20,
    xlab ='Diffusion component (DC)',
    ylab ='Eigenvalue'
  )
}

plot_dm_3D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours, size=0.5, legend = NULL){
  cond <- factor(condition)
  col <- factor(condition)
  levels(col) <- colours
  col <- as.vector(col)
  DCs <- paste("DC",dc, sep="")

  data <- data.frame(
    dm@eigenvectors[,DCs[1]],
    dm@eigenvectors[,DCs[2]],
    dm@eigenvectors[,DCs[3]]
  )
  colnames(data) <- DCs

  plot3d(
    data,
    bg=col,
    col=col,
    size=size,
    box = FALSE
  )

  if (is.null(legend)==FALSE){
    legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")
  }
}

plot_dm_3D_in_2D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours, outline.color = "grey20", size=1, pch = 21, theta = 40, phi = 40, bty ="b"){
  #cond <- factor(condition)
  cols <- factor(condition)
  levels(cols) <- colours
  cols <- as.vector(cols)
  DCs <- paste("DC",dc, sep="")

  data <- data.frame(
    dm@eigenvectors[,DCs[1]],
    dm@eigenvectors[,DCs[2]],
    dm@eigenvectors[,DCs[3]]
  )
  colnames(data) <- DCs

  if(pch == 21){
    scatter3D(x=as.matrix(data[,DCs[1]]), y=as.matrix(data[,DCs[2]]), z=as.matrix(data[,DCs[3]]), bg = cols, col = outline.color, pch=pch, cex = size, xlab = DCs[1], ylab = DCs[2], zlab = DCs[3], theta =theta, phi = phi, bty = bty)
  } else {
    scatter3D(x=as.matrix(data[,DCs[1]]), y=as.matrix(data[,DCs[2]]), z=as.matrix(data[,DCs[3]]),
              colvar =as.numeric(factor(condition)), col = colours, pch=pch, cex = size,
              xlab = DCs[1], ylab = DCs[2], zlab = DCs[3], theta =theta, phi = phi, bty = bty)
  }

}

#== dietseurat while remain reduction information-------------------

dietFun <- function(seurat){
  reduction <- seurat@reductions
  seurat <- DietSeurat(seurat)
  seurat@reductions <- reduction
  return(seurat)
}

#== runTricycle to inference cellcycle information-----------------
runCellcycle <- function(seurat, species = "mouse", gname.type ="SYMBOL"){
  require(tricycle)
  sce <- as.SingleCellExperiment(seurat)
  sce <- project_cycle_space(sce,species = species,gname.type = gname.type)
  sce <- estimate_cycle_position(sce)
  sce <- estimate_Schwabe_stage(sce,
                                gname.type = gname.type,
                                species = species)
  seurat$CCStage <- sce$CCStage
  seurat$tricyclePosition <- sce$tricyclePosition
  return(seurat)
}

# Wrapped function for easy processing scRNA datasets
process_seurat_data <- function(data_directory, min_cells = 3, min_features = 500, qc_mito = 20,qc_rna = 300,species = "Mouse") {
  require(Seurat)
  require(dplyr)
  require(purrr)
  # List all files in the specified directory
  fileNames <- list.files(data_directory, full.names = TRUE)
  fileShort <- list.files(data_directory, full.names = FALSE)

  # Read the 10X data
  countList <- lapply(fileNames, Read10X)
  print("creating seurat object")
  # Create Seurat objects with quality control
  seuratList <- map2(countList, fileShort, function(x, y) {
    CreateSeuratObject(x, min.cells = min_cells, min.features = min_features, project = y)
  })

  # Apply quality control function
  seuratList <- lapply(seuratList, qcFun,percent_mito=qc_mito,nFeature_RNA=qc_rna,Species=species)

  # Merge all Seurat objects
  seuratMerge <- merge(seuratList[[1]], seuratList[2:length(seuratList)])
  print("Running harmony")
  # Run Harmony
  seuratMerge <- runharmony(seuratMerge)

  return(seuratMerge)
}

#== rename gene name-------------------
convertGene <- function(path, symbols) {
  load(path)
  mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])

  humansymbol2mousesymbol = mapper(geneinfo_2022, "symbol_mouse", "symbol")
  converted_symbols = humansymbol2mousesymbol[symbols %>% as.character()]

  return(converted_symbols)
}

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
