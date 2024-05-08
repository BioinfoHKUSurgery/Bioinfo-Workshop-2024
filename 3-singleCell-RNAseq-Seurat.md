# Single cell/nuclei RNA-sequencing (scRNA-seq/snRNA-seq)
This hands-on session will cover the basic workflow of single cell/nuclei RNA-sequencing after alignment, from quality control processing, dimension reduction visualization, and clustering.

<img src="https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/830cc3d1-63f6-46b8-a662-ea64faea95bc" width=600 >

Figure modified from: https://hbctraining.github.io/scRNA-seq_online/lessons/postQC_workflow.html


## Setup a new working directory
First, we will start by setting a working directory to help organize the data, codes and output used/generated in this hands-on session. You can create a folder using *File Explorer* in Windows or *Finder* in mac:
- Create a folder named `scRNAseq_test` in your preferred directory
- Create a `Data` folder under `scRNAseq_test`

From RStudio, we can then set the working directory as `scRNAseq_test` via `Session > Set Working Directory > Choose Directory`.

## Installation
### [Seurat v5](https://satijalab.org/seurat/)
Seurat can be installed like other packages in R using: 
```r
install.packages('Seurat')
library(Seurat)
library(ggplot2)

setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("presto","glmGamPoi"))
```
    
## Dataset
In this tutorial, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. The data can be downloaded from [here](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz). There are 2,700 single cells that were sequenced on the Illumina NextSeq 500.

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/0ed7d89c-cda4-49d0-a91d-bae0f06f4376)

The count matrix were obtained after alignment to transcriptome using Cell Ranger.

## 0. Cell Ranger
Cell Ranger is an analysis pipeline for processing the Chromium 10X single-cell data. It includes functions to align reads, generate feature-barcode matrices, perform clustering, integration and other secondary analysis, and more. Usually, it is used for alignment to create a feature-barcode matrix that can be further processed by `Seurat`.

### *0.1 Cell Ranger workflow* 
The detailed workflow of Cell Ranger can be found [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-what-is-cell-ranger#workflows).

#### One sample, one GEM well, one flow cell
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/30e8b935-ccb7-41a7-b545-839e065c2377)
#### One sample, one GEM well, multiple flow cells
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/fd9a333c-d5de-4308-8a13-9674e1e0576a)
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/9101a050-ae0f-4169-a246-64f4b1eea733)

The Cell Ranger workflow starts with demultiplexing the raw base call (BCL) files for each flow cell directory. You may use `cellranger mkfastq` or one of Illumina's demultiplexing software `bcl2fastq`. 

If you are beginning with FASTQ files that have already been demultiplexed, you can directly run `cellranger count` for alignment, filtering, barcode counting, and UMI counting to generate feature-barcode matrices. 

### *0.2 Cell Ranger count output*
The `outs` folder contains the pipeline output files that can be used for downstream analysis in `Seurat`.
<img src="https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/16cb5ed9-6341-4224-97e1-f7a777be0f3e" width=600 >

#### **0.2.1 Feature Barcode Matrices (MEX Format)**:
It contains gzipped TSV files with feature and barcode sequences corresponding to row and column indices respectively.
<img src="https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/a8698755-171c-45f2-86be-e34641f3ce50" width=300 >

```
filtered_feature_bc_matrix
├── hg19
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```
`filtered_feature_bc_matrix`: Includes only cell associated barcodes that cellranger determined as coming from real cells instead of background. Background barcodes are only included in the `raw_feature_bc_matrix`. For more information: https://kb.10xgenomics.com/hc/en-us/articles/115003480523-How-are-barcodes-classified-as-cell-associated-.
- matrix.mtx.gz: reads count as sparse matrices where each row indicates one feature and each column indicates one cell. The row and column indices correspond to `features.tsv.gz` and `barcodes.tsv.gz` files respectively
- `features.tsv.gz`: feature ID (i.e. Ensembl gene ID) and gene name
- `barcodes.tsv.gz`: barcode sequences

#### **0.2.2 Molecule Info H5**:
More recent versions of cellranger now also output using the .h5 file format, which can be read in using the `Read10X_h5()` function in `Seurat`.

#### **0.2.3 BAM**:
It also outputs an indexed BAM file containing position-sorted reads aligned to the genome and transcriptome, as well as unaligned reads. It is required for some of the secondary analyses, like RNA velocity.

#### **0.2.4 Web Summary**:
A summary HTML file contains `summary` metrics and automated secondary `Analysis` results. 

##### 0.2.4.1 Key metrics
`Estimated Number of Cells`: number of barcodes associated with cells

`Mean Reads per Cell`: total number of sequenced reads divided by the number of cells. A minimum of 20,000 read pairs per cell is recommended.

`Median Genes per Cell`: median number of genes detected per cell-associated barcode, which is dependent on cell type and sequencing depth.   

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/6056d707-6a81-4023-8354-a3cf1a748414)

##### 0.2.4.2 Barcode rank plot
The most important and informative plot in the Gene Expression Web Summary is the `Barcode Rank Plot` under the `Cells` dashboard, which shows the distribution of UMI counts in barcodes and barcodes inferred to be associated with cells. The y-axis is the number of UMI counts mapped to each barcode and the x-axis is the number of barcodes below that value. 

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/eb7edced-69c1-4c04-9199-3ec1dde46200)

**Typical sample (left)**: distinctive shape, which is referred to as a "cliff and knee". The blue-to-gray transition (green arrow) is referred to as the cliff; the solid gray (blue arrow) is the knee. The steep cliff is indicative of good separation between the cell-associated barcodes and the barcodes associated with empty GEMs. 

**Heterogeneous sample (right)**: bimodal plot with two cliff and knee distributions. However, there should still be clear separation between the barcodes called as cells (blue) and barcodes called as background (gray).

More examples on the barcode rank plot can be found under [Technical Note: Interpreting Cell Ranger Web Summary Files for Single Cell Gene Expression Assays](https://www.10xgenomics.com/support/single-cell-gene-expression/documentation/steps/sequencing/interpreting-cell-ranger-web-summary-files-for-single-cell-gene-expression-assays). 

##### 0.2.4.3 Sequencing metrics

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/afab78f2-848e-48e1-a361-947ec425b913)

##### 0.2.4.4 Mapping metrics
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/d00256ae-0603-4e12-8c4f-24bb90d425e1)

> [!WARNING]  
> Warning messages are reported if there are potential issues in quality of the data. 

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/f101c4f1-e4d8-4238-9f6f-73f659411eae)

Further details on quality assessment based on web summary.html can be found [here](https://www.10xgenomics.com/analysis-guides/quality-assessment-using-the-cell-ranger-web-summary)

## 1. Analysis using Seurat
This part of the tutorial is adapted from [Vignette in Seurat v5](https://satijalab.org/seurat/articles/pbmc3k_tutorial) with some modifications.

### 1.1 Create Seurat object
We next read in the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. For example, in Seurat v5, the count matrix is stored in `pbmc[["RNA"]]$counts`.

```r
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "Data/filtered_gene_bc_matrices/hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
## 1 layer present: counts
```
```r
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
## 3 x 30 sparse Matrix of class "dgCMatrix"
##  [[ suppressing 30 column names ‘AAACATACAACCAC-1’, ‘AAACATTGAGCTAC-1’, ‘AAACATTGATCAGC-1’ ... ]]
##                                                                   
## CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
## TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
## MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .
```
The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.

### 1.2 Data processing - Quality control 
The steps below represent the standard pre-processing workflow for the selection and filtration of cells based on QC metrics. A few QC metrics commonly used by the community include:

- `nFeature_RNA`: Number of unique genes detected in each cell
    - Low-quality cells or empty droplets will often have very few genes
    - Cell doublets or multiplets may exhibit an aberrantly high gene count
- `nCount_RNA`: Total number of molecules detected within a cell (correlates strongly with unique genes)
- `percent.mt`: Percentage of reads that map to the mitochondrial genome
    - Low-quality / dying cells often exhibit extensive mitochondrial contamination

`nFeature_RNA` and `nCount_RNA` are automatically calculated during CreateSeuratObject(). We can calculate mitochondrial QC metrics `percent.mt` with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features (e.g. starting with MT- as a set of mitochondrial genes).
```r
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

You can find the QC metrics stored in the object meta data
```r
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
##                  orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACATACAACCAC-1     pbmc3k       2419          779  3.0177759
## AAACATTGAGCTAC-1     pbmc3k       4903         1352  3.7935958
## AAACATTGATCAGC-1     pbmc3k       3147         1129  0.8897363
## AAACCGTGCTTCCG-1     pbmc3k       2639          960  1.7430845
## AAACCGTGTATGCG-1     pbmc3k        980          521  1.2244898
```

In the example below, we visualize QC metrics, and use these to filter cells.
- We filter cells that have unique feature counts over 2,500 or less than 200
- We filter cells that have >5% mitochondrial counts
```r
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/316ecc9d-97e7-4405-b852-b783a78afb30)

```r
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/4e6fc677-aab6-4ad4-897f-df57c5810cf6)

```r
# Retain cells with nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
## An object of class Seurat 
## 13714 features across 2638 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
##  1 layer present: counts
```

### 1.3 Normalization and scaling (Replace with SCTransform)
After removing unwanted cells from the dataset, the next step is to normalize the data. 

Standard workflow runs `NormalizeData`, `ScaleData`, and `FindVariableFeatures`:
- `NormalizeData` : normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result
- `ScaleData`: applies a linear transformation (‘scaling’) prior to dimensional reduction techniques like principal component analysis (PCA), e.g. standardization to mean of 0 and variance of 1 across cells
- `FindVariableFeatures`: calculates highly variable genes and focuses on these for downstream analysis to highlight biological signal in single-cell datasets. Highly variable genes are those genes that exhibit high cell-to-cell variation in the dataset (i.e, highly expressed in some cells but lowly expressed in others). 

A newer modeling framework for normalization and variance stabilization, called 'SCTransform`, replaces `NormalizeData`, `ScaleData`, and `FindVariableFeatures` and is proposed to improve common downstream analytical tasks such as variable gene selection, dimensional reduction, and differential expression. 

After running `SCTransform`, transformed data will be available in the `SCT` assay. During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage,

```r
library(sctransform)

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc
```

`pbmc[["SCT"]]$scale.data`: contains the residuals (normalized values), and is used directly as input to PCA. Please note that this matrix is non-sparse, and can therefore take up a lot of memory if stored for all genes. To save memory, we store these values only for highly variable genes (3,000 by default), by setting the return.only.var.genes = TRUE by default in the SCTransform() function call.


### 1.4 Dimension reduction and visualization
#### 1.4.1 Perform linear dimensional reduction
Identifying the true dimensionality of a dataset – can be challenging. It is recommended to use multiple approaches. The first is more supervised, exploring principal components (PCs) to determine relevant sources of heterogeneity. 

So we first perform linear dimensional reduction--PCA--on the scaled data. By default, only the previously determined variable features (3,000 for `SCTransform`) are used as input. 

```r
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)   # pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```
```r
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
## PC_ 1 
## Positive:  FTL, LYZ, FTH1, CST3, S100A9 
## Negative:  MALAT1, RPS27A, CCL5, RPS6, LTB 
## PC_ 2 
## Positive:  NKG7, CCL5, GZMB, GNLY, GZMA 
## Negative:  HLA-DRA, CD74, CD79A, HLA-DPB1, HLA-DQA1 
## PC_ 3 
## Positive:  S100A8, S100A9, LYZ, FTL, RPS12 
## Negative:  CD74, HLA-DRA, CD79A, HLA-DPB1, HLA-DQA1 
## PC_ 4 
## Positive:  FCGR3A, LST1, FCER1G, AIF1, IFITM3 
## Negative:  S100A8, S100A9, LYZ, LGALS2, CD14 
## PC_ 5 
## Positive:  GNLY, GZMB, FGFBP2, FCGR3A, PRF1 
## Negative:  CCL5, GPX1, PPBP, PF4, SDPR 
```
Seurat provides several useful ways of visualizing both cells and features that define the PCA, including `VizDimReduction`, `DimPlot`, `DimHeatmap`, and 
```r
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/e4a29f22-d13f-4f92-8553-8b9a3df3fcbe)
```r
DimPlot(pbmc, reduction = "pca")
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/78990975-5f13-4092-b4f0-586743c6f095)

In particular `DimHeatmap` allows for easy exploration of the primary sources of heterogeneity in a dataset. Both cells and features are ordered according to their PCA scores. Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets.
```r
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/18761652-7cfe-4120-a728-0a28551515e5)

To determine more quantatitively the number of PCs to be used for downstream analysis, we can use `ElbowPlot` to see when the SD reaches the plateau. 
```r
ElbowPlot(pbmc, ndims=30)
```
When using sctransform, more PCs can be included as the sctransform workflow performs more effective normalization which removes technical effects from the data. Users are encouraged to repeat downstream analyses with a different number of PCs and see how the clustering changes.

#### 1.4.2 Perform non-linear dimensional reduction (UMAP/tSNE) and clustering
Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. 

For clustering, Seurat applies a graph-based clustering approach to construct a K-nearest neighbor (KNN) graph based on the euclidean distance in PCA space. This step is performed using the `FindNeighbors` function, and takes as input the previously defined dimensionality of the dataset (first 30 PCs). Next, `FindClusters` applies modularity optimization techniques such as the Louvain algorithm (default) or SLM to iteratively group cells together. The clusters can be found using the Idents() function.

Cells that are grouped together within graph-based clusters determined below should co-localize on these dimension reduction plots.

```r
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE)
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/bdb53b18-619c-45d5-bae5-e3ac1b3bc5b8)

## 6. Marker genes identification
Seurat can help you find markers that define clusters via differential expression (DE). By default, it identifies both positive and negative markers of a single cluster (specified in ident.1), compared to all other cells.  `FindAllMarkers` automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

In Seurat v5, we use the presto package (as described here and available for installation here), to dramatically improve the speed of DE analysis, particularly for large datasets. For users who are not using presto, `min.pct` and `logfc.threshold` parameters can be set to increase the speed of DE testing.

```r
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
##               p_val avg_log2FC pct.1 pct.2     p_val_adj
## RPS27  5.342505e-116  0.7425496  1.00 0.998 6.688283e-112
## S100A4 2.224400e-101 -2.4259874  0.63 0.889  2.784726e-97
## RPL32  1.462416e-100  0.5981817  1.00 0.999  1.830799e-96
## RPS6    1.021260e-95  0.6212855  1.00 1.000  1.278515e-91
## RPS12   1.170400e-94  0.6631844  1.00 1.000  1.465223e-90
```
```r
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
```
```r
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

## # A tibble: 6,420 × 7
## # Groups:   cluster [12]
##       p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   
##       <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
## 1 1.50e-117       1.61 0.979 0.63  1.88e-113 0       LTB    
## 2 1.20e-111       1.43 0.938 0.445 1.51e-107 0       IL32   
## 3 3.22e- 95       1.62 0.759 0.3   4.03e- 91 0       IL7R   
## 4 8.71e- 79       1.10 0.895 0.416 1.09e- 74 0       CD3D   
## 5 1.83e- 73       2.46 0.421 0.101 2.29e- 69 0       AQP3   
## 6 2.69e- 72       1.05 0.929 0.593 3.37e- 68 0       LDHB   
## 7 1.50e- 70       3.89 0.219 0.018 1.88e- 66 0       TNFRSF4
## 8 7.16e- 59       1.71 0.596 0.237 8.97e- 55 0       CD2    
## 9 1.14e- 56       2.44 0.292 0.057 1.42e- 52 0       CD40LG 
## 10 1.90e- 55       1.05 0.779 0.395 2.38e- 51 0       CD3E   
## # ℹ 6,410 more rows
## # ℹ Use `print(n = ...)` to see more rows
```

Seurat has several tests for identifying DEG which can be set with the `test.use` parameter. For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).

```r
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
```

We include several tools for visualizing marker expression. VlnPlot() (shows expression probability distributions across clusters), and FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring RidgePlot(), CellScatter(), and DotPlot() as additional methods to view your dataset.

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

These are now standard steps in the Seurat workflow for visualization and clustering

```r
# Visualize canonical marker genes as violin plots.
VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5"), pt.size = 0.2, ncol = 3)
# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(pbmc, features = c("CD8A", "GZMK", "CCL5"), pt.size = 0.2, ncol = 3)

## B cell cluster - TCL1A, FCER2
VlnPlot(pbmc, features = c("TCL1A", "FCER2"), pt.size = 0.2, ncol = 2)
FeaturePlot(pbmc, features = c("TCL1A", "FCER2"), pt.size = 0.2, ncol = 3)
```

