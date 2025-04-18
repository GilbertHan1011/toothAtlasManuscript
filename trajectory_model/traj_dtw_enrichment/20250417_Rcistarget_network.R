#############################################################
# Building a Gene Regulatory Network (GRN) using RcisTarget
#############################################################

devtools::install_url("https://scenic.aertslab.org/scenic_paper/downloads/Rpackages/RcisTarget_0.99.0.tar.gz")
# Install required packages (uncomment if needed)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("RcisTarget", "DT", "visNetwork"))
#
# # Install the appropriate motif database for your organism
# # For human (hg19):
# BiocManager::install("RcisTarget.hg19.motifDatabases.20k")
# # For mouse (mm9):
BiocManager::install("RcisTarget.mm9.motifDatabases.20k")
install.packages("~/Desktop/disk1/database/RcisTarget/RcisTarget.mm9.motifDatabases.20k_0.1.1.tar.gz")
library(RcisTarget.mm9.motifDatabases.20k)
# Load required libraries
library(RcisTarget)
library(data.table)
library(reshape2)
library(visNetwork)
library(DT)

##########################
# 1. Prepare your gene list
##########################

# Option 1: Load your gene list from a file
# geneListFile <- "path/to/your/gene_list.txt"
# genes <- read.table(geneListFile, header=FALSE, stringsAsFactors=FALSE)[,1]



# Create a named list (required format for RcisTarget)
geneLists <-mine_conservation$gene
geneLists <- list(geneListName=geneLists)
cat("Gene list prepared with", length(genes), "genes\n")

##########################
# 2. Load motif databases
##########################

# Choose the appropriate database for your organism
# For human (hg19):
library(RcisTarget.hg19.motifDatabases.20k)
data("mm9_500bpUpstream_motifRanking")
data("mm9_10kbpAroundTss_motifRanking")
#data(hg19_500bpUpstream_motifRanking)  # Promoter regions
# Alternative: data(hg19_10kbpAroundTss_motifRanking)  # Wider regions
data("motifAnnotations_mgi_v9")
data("mm9_direct_motifAnnotation")
# Load motif annotation
#data(motifAnnotations_hgnc)
ranking <- importRankings("~/Desktop/disk1/database/RcisTarget/mm9-tss-centered-10kb-7species.mc9nr.feather")
motifRankings <- ranking
motifAnnotation <- motifAnnotations_mgi_v9

cat("Loaded motif database with", nrow(motifAnnotation), "motifs\n")

##########################
# 3. Run motif enrichment
##########################

# Full RcisTarget analysis (one-step approach)
motifEnrichmentTable_wGenes <- cisTarget(
  geneLists,
  mm9_10kbpAroundTss_motifRanking,
  motifAnnot_direct = motifAnnotation
)

# Alternative: Step by step approach (gives more control)
# Step 1: Calculate enrichment scores (AUC)
# motifs_AUC <- calcAUC(geneLists, motifRankings)
#
# # Step 2: Select significant motifs and add TF annotation
# motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
#                                           motifAnnot=motifAnnotation,
#                                           nesThreshold=3.0)
#
# # Step 3: Identify significant genes for each motif
# motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
#                                                   geneSets=geneLists,
#                                                   rankings=motifRankings,
#                                                   method="aprox")

cat("Found", nrow(motifEnrichmentTable_wGenes), "enriched motifs\n")


motifs_AUC <- calcAUC(geneLists, mm9_10kbpAroundTss_motifRanking)
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, motifAnnot_direct=mm9_direct_motifAnnotation)
motifEnrichmentTable_wGenes$motif


signifMotifNames <- motifEnrichmentTable_wGenes$motif[1:10]

incidenceMatrix <- getSignificantGenes(geneLists$geneListName,
                                       mm9_10kbpAroundTss_motifRanking,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000,
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")

library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),
                    label=c(motifs, genes),
                    title=c(motifs, genes), # tooltip
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE,
                                        nodesIdSelection = TRUE)
write.csv(edges,"results/trajectory/20250415_trajdtw_fit/20250417_minerlization_network.csv")

mm9_10kbpAroundTss_motifRanking

motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                           motifAnnot_direct=mm9_direct_motifAnnotation,
                                           motifAnnot_inferred=mm9_inferred_motifAnnotation)


anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_direct, motifEnrichmentTable$geneSet),
                      function(x) unique(unlist(strsplit(x, "; "))))


motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings,
                                                   geneSets=geneLists)
dim(motifEnrichmentTable_wGenes)

anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_direct, motifEnrichmentTable$geneSet),
                      function(x) unique(unlist(strsplit(x, "; "))))

View(motifEnrichmentTable)

wri
