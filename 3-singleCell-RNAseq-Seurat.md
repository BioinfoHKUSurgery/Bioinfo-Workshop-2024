# Single cell/nuclei RNA-sequencing (scRNA-seq/snRNA-seq)
This hands-on session will cover the basic workflow of single cell/nuclei RNA-sequencing after alignment, from quality control processing, dimension reduction visualization, clustering to integration, following the [Vignette in Seurat v5](https://satijalab.org/seurat/articles/pbmc3k_tutorial).

<img src="https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/830cc3d1-63f6-46b8-a662-ea64faea95bc" width=600 >

Figure modified from: https://hbctraining.github.io/scRNA-seq_online/lessons/postQC_workflow.html


## Setup a new working directory
First, we will start by setting a working directory to help organize the data, codes and output used/generated in this hands-on session. You can create a folder using *File Explorer* in Windows or *Finder* in mac:
- Create a folder named `scRNAseq_test` in your preferred directory
- Create a `Data` folder under `scRNAseq_test`

From RStudio, we can then set the working directory via `Session > Set Working Directory > Choose Directory`.

## Installation
### [Seurat v5](https://satijalab.org/seurat/)
Seurat can be installed like other packages in R using: 
```r
install.packages('Seurat')
library(Seurat)
```

## Dataset
In this tutorial, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. The raw data can be downloaded from [here](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz). There are 2,700 single cells that were sequenced on the Illumina NextSeq 500.

The count matrix were obtained after alignment to transcriptome using Cell Ranger.

### Data format
- Count matrix
  - barcodes.tsv
  - genes.tsv
  - matrix.mtx

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/0ed7d89c-cda4-49d0-a91d-bae0f06f4376)

- .h5 Seurat object

## 1. Create Seurat object
```r
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~\Downloads\pbmc3k_filtered_gene_bc_matrices\filtered_gene_bc_matrices\hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

## 2. Data processing - Quality control 
## 3. Normalization and scaling
## 4. Dimension reduction visualization
## 5. Clustering
## 6. Marker genes identification
## 7. Integration
