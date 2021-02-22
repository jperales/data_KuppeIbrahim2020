## Setting the environment
### Internal variables
set.seed(1234)
OUTDIR <- "./data/CD10negative/"

### Load libraries
library("Matrix")
library("SingleCellExperiment")
library("scran")

## Load data
dat <- Matrix::readMM("./data/CD10negative/kidneyMap_UMI_counts.mtx")
rowDat <- read.table("./data/CD10negative/kidneyMap_UMI_counts_rowData.txt",
		     sep=",", header=TRUE, stringsAsFactors = FALSE)
colDat <- read.table("./data/CD10negative/kidneyMap_UMI_counts_colData.txt", 
		     sep=",",header=TRUE, stringsAsFactors = FALSE)
umapCoords <- read.table("data/CD10negative/kidneyMap_UMI_umapCoords.csv",
			 sep=",")
umapCoords <- as.matrix(umapCoords)

# Genes
rownames(dat) <- rowDat$ENSEMBL.ID
rownames(rowDat) <- rowDat$ENSEMBL.ID

# Cells
colnames(dat) <- paste0("cell",1:ncol(dat))
rownames(colDat) <- paste0("cell",1:ncol(dat))
rownames(umapCoords) <- paste0("cell",1:ncol(dat))

# Metafeatures
colnames(umapCoords) <- c("UMAP_1","UMAP_2")

# Summary of cell metadata
## Create a Single-Cell Experiment
sce <- SingleCellExperiment(assays=list("counts"=dat),
			    colData=colDat,
			    rowData=rowDat)

## Normalize data
#NOTE: Params defined by M.Ibrahim
sce = scran::computeSumFactors(sce, 
			       sizes = seq(10, 200, 20),  
			       clusters = sce$Annotation.Level.3, 
			       positive = TRUE)
sce <- logNormCounts(sce)

## Add original UMAP coords
reduceDim(sce, "UMAP") <- umapCoords

# Save data
saveRDS(sce, file=paste0(OUTDIR,"/sce.rds"))
