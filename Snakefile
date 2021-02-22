
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
# https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html

HTTP = HTTPRemoteProvider()

rule all:
	input:
		"data/CD10negative/kidneyMap_UMI_counts_colData.txt",
		"data/CD10negative/kidneyMap_UMI_counts_rowData.txt",
		"data/CD10negative/kidneyMap_UMI_counts.mtx",
		"data/CD10negative/kidneyMap_UMI_umapCoords.csv",
		"data/CD10negative/sce.rds",
		"data/CD10negative/SeuratObject.rds"

### 1 download archived data from publication for CD10negative kidney cells
rule dn_tar:
	input:
		HTTP.remote("zenodo.org/record/4059315/files/Human_CD10negative.tar.gz", keep_local=True)
	output:
		"data/Human_CD10negative.tar.gz"
	message:
		"Downloading raw"
	shell:
		'(test -d data || mkdir data)'
		' && mv {input} {output[0]}'

rule tar_xf:
	input:
		"data/Human_CD10negative.tar.gz"
	output:
		"data/CD10negative/kidneyMap_UMI_counts_colData.txt",
		"data/CD10negative/kidneyMap_UMI_counts_rowData.txt",
		"data/CD10negative/kidneyMap_UMI_counts.mtx"
	message:
		"Extracting achived files from TAR"
	shell:
		'tar -xf {input[0]} -C data'
		' && rm {input[0]}'

rule dn_umapCoords:
	input:
		HTTP.remote("raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_CD10negative_umapCoords.csv", keep_local=True)
	output:
		"data/CD10negative/kidneyMap_UMI_umapCoords.csv"
	message:
		"Downloading original UMAP coords"
	shell:
		'(test -d data/CD10negative || mkdir -p data/CD10negative)'
		' && mv {input} {output[0]}'

### 2 Transform processed data into R objects for scran and seurat
### 2.1 Build SingleCellExperiment
rule build_sce:
	input:
		"data/CD10negative/kidneyMap_UMI_counts_colData.txt",
		"data/CD10negative/kidneyMap_UMI_counts_rowData.txt",
		"data/CD10negative/kidneyMap_UMI_counts.mtx"
	output:
		"data/CD10negative/sce.rds"
	message:
		"Building SingleCellExperiment"
	conda:
		"envs/osca.yaml"
	script:
		"scripts/build_sce.R"

### 2.2 Build Seurat Object
rule build_SeuratObject:
	input:
		"data/CD10negative/kidneyMap_UMI_counts_colData.txt",
		"data/CD10negative/kidneyMap_UMI_counts_rowData.txt",
		"data/CD10negative/kidneyMap_UMI_counts.mtx"
	output:
		"data/CD10negative/SeuratObject.rds"
	message:
		"Building SeuratObject"
	conda:
		"envs/seurat.yaml"
	script:
		"scripts/build_SeuratObject.R"

