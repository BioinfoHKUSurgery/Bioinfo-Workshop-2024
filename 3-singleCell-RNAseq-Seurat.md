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

## 1. Create Seurat object
We next use the count matrix to create a Seurat object. 
```r
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~\Downloads\pbmc3k_filtered_gene_bc_matrices\filtered_gene_bc_matrices\hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# An object of class Seurat 
# 13714 features across 2700 samples within 1 assay 
# Active assay: RNA (13714 features, 0 variable features)
# 1 layer present: counts
```
```r
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
# 3 x 30 sparse Matrix of class "dgCMatrix"
#  [[ suppressing 30 column names ‘AAACATACAACCAC-1’, ‘AAACATTGAGCTAC-1’, ‘AAACATTGATCAGC-1’ ... ]]
#                                                                   
# CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
# TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
# MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .
```
The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.

## 2. Data processing - Quality control 
The steps below represent the selection and filtration of cells based on QC metrics. A few QC metrics commonly used by the community include:

- The number of unique genes detected in each cell.
- Low-quality cells or empty droplets will often have very few genes
- Cell doublets or multiplets may exhibit an aberrantly high gene count
- the total number of molecules detected within a cell (correlates strongly with unique genes)
- The percentage of reads that map to the mitochondrial genome
- Low-quality / dying cells often exhibit extensive mitochondrial contamination

We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features. We use the set of all genes starting with MT- as a set of mitochondrial genes
```r
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```
The number of unique genes and total molecules are automatically calculated during CreateSeuratObject()
You can find them stored in the object meta data
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

orig.ident:
nCount_RNA:
nFeature_RNA:
percent.mt:

In the example below, we visualize QC metrics, and use these to filter cells.
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


We filter cells that have unique feature counts over 2,500 or less than 200
We filter cells that have >5% mitochondrial counts
```r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```
## 3. Normalization and scaling (Replace with SCTransform)
After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. In Seurat v5, Normalized values are stored in pbmc[["RNA"]]$data.
For clarity, in this previous line of code (and in future commands), we provide the default values for certain parameters in the function call. However, this isn’t required and the same behavior can be achieved with:
```r
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
```
While this method of normalization is standard and widely used in scRNA-seq analysis, global-scaling relies on an assumption that each cell originally contains the same number of RNA molecules. We and others have developed alternative workflows for the single cell preprocessing that do not make these assumptions. For users who are interested, please check out our SCTransform() normalization workflow. The method is described in ourpaper, with a separate vignette using Seurat here. The use of SCTransform replaces the need to run NormalizeData, FindVariableFeatures, or ScaleData (described below.)


## 4. Dimension reduction visualization
We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

Our procedure in Seurat is described in detail here, and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures() function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.
```r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/bf6355d3-93e6-4e4f-b600-be18c8fed443)

Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:

Shifts the expression of each gene, so that the mean expression across cells is 0
Scales the expression of each gene, so that the variance across cells is 1
This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
The results of this are stored in pbmc[["RNA"]]$scale.data
By default, only variable features are scaled.
You can specify the features argument to scale additional features
```r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```
n Seurat, we also use the ScaleData() function to remove unwanted sources of variation from a single-cell dataset. For example, we could ‘regress out’ heterogeneity associated with (for example) cell cycle stage, or mitochondrial contamination i.e.:
```r
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
```
However, particularly for advanced users who would like to use this functionality, we strongly recommend the use of our new normalization workflow, SCTransform(). The method is described in our paper, with a separate vignette using Seurat here. As with ScaleData(), the function SCTransform() also includes a vars.to.regress parameter.
## 5. Clustering
## 6. Marker genes identification
## 7. Integration
