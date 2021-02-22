## Setting the environment
### Internal variables
set.seed(1234)
OUTDIR <- "./data/CD10negative/"

### Load libraries
library("Matrix")
library("Seurat")

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

## Create a Single-Cell Experiment
S <- Seurat::CreateSeuratObject(counts = dat,
				project = "CD10negative")

## Add metadata
S <- AddMetaData(S, colDat)

## Normalize data
S <- Seurat::NormalizeData(S)

## Add original UMAP coords
S[["umap"]] <- CreateDimReducObject(embeddings = umapCoords, 
				       key = "UMAP_", 
				       assay = DefaultAssay(S))

# Save data
saveRDS(S, file=paste0(OUTDIR,"/SeuratObject.rds"))
