### “Omics Data Analysis”

### 1. Introduction to Omics data

Omics approaches generate large-scale molecular biology data, which are analyzed using bioinformatics tools and computational approaches to extract meaningful biological information. Omics has broad applications in various fields, including basic research, clinical diagnostics, drug discovery, and personalized medicine. It encompasses various disciplines, including genomics, transcriptomics, proteomics, metabolomics, and epigenomics. In this session, we will focus on transcriptomics and how bioinformatics approaches can be used to analyse and interpret RNA sequencing data.

### 2. Case study: RNA-seq data from an in-vitro model of proliferative vitreoretinopathy (PVR)

Retinal scarring is a common complication associated with various retinal diseases, including proliferative vitreoretinopathy (PVR). A recent study (<https://www.nature.com/articles/s41467-022-30474-6>) investigated the anti-scarring effects of a thermogel polymer (PEP) on both an in-vivo and in-vitro model of PVR by performing RNA sequencing and identifying the deferentially expressed genes and pathways. In this practical session, we will use R and RStudio to analyse the processed RNA-seq data downloaded from the gene expression omnibus (GEO).

### 2.1 Importing the read count matrix

Open the web link to the dataset on GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176513>

Here you can find a description of the experimental design, links to the raw data, and processed data which in this case is a read count matrix called GSE176513_gene_count_matrix.csv.gz. To find this, scroll down to the bottom and under "supplementary file" click download (ftp).

![](images/clipboard-2351814831.png)

Once downloaded, you can read the file into R using the `read.csv()` function:

```         
count.mat<-read.csv('/Users/pauldavidblakeley/Documents/GSE176513_gene_count_matrix.csv', header = T, sep = ',', row.names = 1)
```

We want to use the gene ids as rownames, so we use the `row.names = 1` option to use the first column as rownames.

### 2.2 Importing the meta data file

The meta data information about each experimental sample can usually be extracted from GEO by using the SRA Run Selector found at the bottom of the page.

However, to save time for this tutorial this file has already been downloaded and reformatted. It can be found on the github page by clicking .......

To read the downloaded file into R use the following:

```         
#read study design file
design=read.csv('/Users/pauldavidblakeley/Documents/Novogene_meta_invitro.csv', header = TRUE, sep = ",", row.names = 1)
```

### 2.3 Using BiomaRt to convert Ensembl IDs to gene names

BiomaRt is an R package that allows for programmatic access to the BioMart databases, which are maintained by the Ensembl project. BioMart databases contain a wealth of information on genes, transcripts, and other genomic features for a wide range of species.

To install BiomaRt and use it to convert Ensembl IDs to gene names, you can use the following code:

```         
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

#use the human biomart
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Extract gene annotations from biomart
t2g <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                            "chromosome_name", "gene_biotype"), mart = human)
```

### 2.3 DESeq2 for normalization and differential expression

We are mainly interested in the 24 hour post treatment samples. The following code will subset the count matrix and the metadata file to ensure they only contain the 24h time point.

```         
counts2<-count.mat[, grep("24h", colnames(count.mat))]
design<-design[grep('24h', rownames(design)),]
```

This tutorial will rely on DESeq2 to perform normalization and differential expression analysis. We shall use the DESeqDataSetFromMatrix from DESeq2 to create a DESeqDataSet object from a count matrix and sample metadata.

First install DESeq2 using BiocManager installer tool

```         
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2"")
```

Then run DESeqDataSetFromMatrix and remove genes with fewer than 2 read counts for increased efficiency and faster run time,

```         
design.fac<-data.frame(design$Treatment, design$Time)
design.fac$cell_time<-paste(design$Treatment, design$Time, sep = '_')

cds2=DESeqDataSetFromMatrix(countData=counts2, colData=DataFrame(design.fac), design= ~ design.Treatment)

cds2 <- cds2[ rowSums(counts(cds2)) > 1, ]
```

Next we will run a normaization method provided in DESeq2 called variance stabilizing transformation. The VST is a popular approach for normalizing count data and reducing the dependency of the variance on the mean, which is often observed in RNA-seq data.

```         
#vst
dseq.vst<-vst(cds2, blind = TRUE, fitType = "local")
vst.mat<-assay(dseq.vst)
vst.mat2<-t(vst.mat)
```

### 2.4 Sample correlation heatmap and PCA

Before we proceed with differential expression significance testing provided by DESeq2, we would first like to check the samples for any outliers, check the within- and between- group variances, and to see if the PEP treatment is having an effect on the gene expression. We will use the VST normalized values perform cl

First calculate euclidean distance between samples. Other distance metrics such as Pearson correlation may also be used instead.

```         
sampleDists <- dist(vst.mat2, method = 'euclidean')
sampleDistMatrix <- as.matrix(sampleDists)
```

Then use the pheatmap package to cluster the samples and plot a heatmap

```         
install.packages('pheatmap')
install.packages('RColorBrewer')
library('pheatmap')
library('RColorBrewer')
library('ggplot2')


rownames(sampleDistMatrix) <- dseq.vst$cell_time
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, clustering_method = 'ward',
         col=colors)
```

This could produce the following heatmap including dendrogram:

![](images/clipboard-1522008039.png)

DESeq2 has a `plotPCA()` function which we use to generate a PCA based on the top 8000 most variable genes.

```         
data <- plotPCA(dseq.vst, intgroup=c("design.Treatment", "design.Time"), returnData=TRUE,ntop=8000)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2,label=design.Treatment,color=design.Treatment), size=1) +
  geom_point(size=3) +
  geom_text(size=2.5,nudge_y = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
```

![](images/clipboard-2982271965.png)

### 2.5 DE testing

Next we use the `DESeq` function to fit a negative binomial model and estimate dispersion and differential expression of genes between the two conditions **PolymerTNT_24h** vs **TNT_24h**

```         
dseq2.out<-DESeq(cds2)
```

The `results` function returns the results as a table containing log2 fold change, p-values, adjusted p-values and other statistics.

```         
res.CD <- as.data.frame(results(dseq2.out, contrast = c("design.Treatment", "PolymerTNT_24h", "TNT_24h")))
```

It should return a table like this

```         
                  baseMean log2FoldChange      lfcSE        stat        pvalue          padj
ENSG00000135454 7.867856e+01    2.258361003 0.24294700  9.29569420  1.462487e-20  6.311785e-19
ENSG00000221866 3.751892e+01   -0.250353584 0.28150129 -0.88935148  3.738142e-01  5.710085e-01
ENSG00000225131 9.244403e+00    0.133531409 0.60906066  0.21924156  8.264619e-01  9.107814e-01
ENSG00000287029 1.260023e+00    1.919200590 1.82661194  1.05068873  2.934016e-01            NA
ENSG00000252652 3.324833e-01   -1.991873684 3.27209542 -0.60874560  5.426931e-01            NA
ENSG00000138175 6.988709e+02   -0.055581295 0.07925023 -0.70133918  4.830914e-01  6.682258e-01
ENSG00000115977 1.192017e+03    0.015010685 0.07917980  0.18957719  8.496405e-01  9.233423e-01
```

We also want to add the average normalized expression values to the results table and add a column for gene names and chromosome location:

```         
#add average expression column
base.means<-sapply( levels(dseq2.out$design.Treatment), function(lvl) rowMeans( counts(dseq2.out,normalized=TRUE)[,dseq2.out$design.Treatment == lvl, drop=F] ) )
res.CD<-merge(as.data.frame(res.CD), base.means, by.x='row.names', by.y='row.names')

#add gene information column
res.CD<-merge(as.data.frame(res.CD), t2g, by.x='Row.names', by.y='ensembl_gene_id')
```

Finally to write the DESeq2 results to a CSV file we can use:

```         
#write deseq results
res.CD<-as.data.frame(res.CD)
write.csv(res.CD, file='deseq2_PVR_CvsD.csv')
```

### 2.6. Using volcano plots to display DE results

Volcano plots are commonly used in RNA-seq data analysis to visualize the statistical significance (usually represented by p-values or adjusted p-values) and the magnitude of differential gene expression between two conditions. In R, you can create volcano plots using various packages, such as ggplot2, base R, or EnhancedVolcano.

Here we use the volcano plot to show the highly DE genes between Poly(PEP) and control.

First we use `ifelse` to specify only genes with p-value \< 0.05 should be coloured, the rest remain grey

```         
res.CD$col<-ifelse(res.CD$padj<0.05, 'red4', 'grey30')
```

And to convert Pval into -log2Pval:

```         
res.CD$minl2pval<--log2(res.CD$pvalue)
```

Now ready to generate the volcano plot using ggplot2

```         
ggplot(data = res.CD, aes(y=minl2pval, x=log2FoldChange, col=col)) + 
  geom_point(colour=res.CD$col, size=1.5) + 
  geom_text(aes(log2FoldChange, minl2pval), size =1.75, col='red4', nudge_y = 10, label = ifelse(res.CD$padj<0.05 & (res.CD$log2FoldChange< -0.2 | res.CD$log2FoldChange>0.2), res.CD$external_gene_name, '')) +
  ylim(-0,500) + coord_cartesian(clip = 'off', expand = T, ylim = c(0,600)) + xlim(-5,5) + xlab('log2fc') + ylab('-log2(p-value)') +
  theme(text = element_text(size=5)) + theme_minimal()
```

![](images/clipboard-2010409566.png)

From looking at the top-right and top-left, you can identify the most highly differential expressed genes. We observe that some genes such as `GSR, ME1, SLC6A6` are components of the NRF2 signalling pathway, which id a critical cellular defense mechanism against oxidative stress. We want to see whether the NRF2 pathway is indeed upregulated after poly(PEP) treatment, so we generate a volcano plot highlighting all members of this signalling pathway.

First we download the list of NRF2 pathway genes from the github page and load into R

```         
nrf_genes<-read.csv('/Users/pauldavidblakeley/Downloads/NCOMMS-21-46079B_zip-2/Source data 2_RNA-Seq Data/nrf2_pathway/nrf.txt', header = F)
nrf_genes<-as.vector(nrf_genes$V1)
```

Then generate the volcano plot by specifying with `ifelse` that we only want to colour and label NRF2 pathways genes. Using the `alpha` option we can fade out the non NRF2 pathway genes to place more emphasis on the genes of interest.

```         

res.CD$col<-ifelse(res.CD$external_gene_name %in% nrf_genes, 'blueviolet', 'grey30')
res.CD$minl2pval<--log2(res.CD$pvalue)

ggplot(data = res.CD, aes(y=minl2pval, x=log2FoldChange, col=col)) + 
  geom_point(colour=res.CD$col, size=2.5, alpha=ifelse(res.CD$col == 'grey30', 0.11, 0.84)) + 
  geom_text(aes(log2FoldChange, minl2pval), size =3, col='black', nudge_y = 3, alpha=0.9, label = ifelse(res.CD$col=='blueviolet', res.CD$external_gene_name, '')) +
  coord_cartesian(clip = 'off', expand = T, ylim = c(0,400))+ xlim(-8,8) + xlab('log2fc') + ylab('-log2(p-value)') +
  theme(text = element_text(size=20)) + geom_hline(yintercept = 1.30103, linetype="dashed")+ theme_minimal()
```

![](images/clipboard-21520051.png)
