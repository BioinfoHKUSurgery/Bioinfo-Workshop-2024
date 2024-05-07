# Single cell/nuclei RNA-sequencing (scRNA-seq/snRNA-seq)

## Installation
### Setup a new working directory

### Seurat v5 ()
  
```{R}
install.packages('Seurat')
library(Seurat)
```

https://satijalab.org/seurat/articles/pbmc3k_tutorial

```{R}
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
