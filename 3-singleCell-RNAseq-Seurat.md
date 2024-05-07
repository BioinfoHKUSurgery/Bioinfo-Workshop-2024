# Single cell/nuclei RNA-sequencing (scRNA-seq/snRNA-seq)
This hands-on session will cover the basic workflow of single cell/nuclei RNA-sequencing after alignment, from quality control processing, dimension reduction visualization, clustering to integration, following the [Vignette in Seurat v5](https://satijalab.org/seurat/articles/pbmc3k_tutorial).
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/830cc3d1-63f6-46b8-a662-ea64faea95bc)
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




```r
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~\Downloads\pbmc3k_filtered_gene_bc_matrices\filtered_gene_bc_matrices\hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

## Dataset
- Software and data to be downloaded
- Seurat
`R` is a free, open-source programming language used for statistical computing and graphical presentation while `RStudio` is an integrated development environment (IDE) featureing tools for plotting, viewing history, debugging and managing your workspace for R and Python. 

In the era of data science, there are several pros and cons of using R programming in analyzing clinical or biological data:

<img src="https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/3d928712-a632-43e1-b7da-0bcdfc445d96" width=600 >


## 1. Installation
For this and the coming tutorials, you will need to download the data in tab-delimited text file format
- [Pheno.txt](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/blob/be74dc2428d340c467562eab9ce580caea905a5d/Data/1-Introduction-to-R/Pheno.txt) 
- [Mutations.txt](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/blob/be74dc2428d340c467562eab9ce580caea905a5d/Data/1-Introduction-to-R/Mutations.txt)

and install `R` and `RStudio`
- R: [https://cloud.r-project.org/](https://cloud.r-project.org)
- RStudio: [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)
