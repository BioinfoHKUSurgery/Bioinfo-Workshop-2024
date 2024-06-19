# 6-Machine learning practical session

## Sparse Partial Least Squares Regression Discriminant Analysis (sPLS-DA)

Sparse Partial Least Squares Discriminant Analysis (sPLS-DA) is a machine learning algorithm used for classification and dimensionality reduction. sPLS-DA is particularly useful for analyzing high-dimensional data, such as gene expression data, proteomics or metabolomics experiments as it includes a regulariztion compomnent to identify the most important variables to be included in the model

In sPLS-DA, the input data matrix is decomposed into a small number of latent variables (linear combinations of the original variables) which maximally seperate the data classes when used in the classification model.

### Installation

First lets install the mixOmics R package using BiocManager and the load the package.

```
BiocManager::install('mixOmics')
library(mixOmics)
```

### The Small Round Blue Cell Tumours (SRBCT) dataset

The Small Round Blue Cell Tumours (SRBCT) dataset was generating by studying the expression levels of a set of genes to test for markers of certain tumour types (Khan et al, 2001). SRBCTs are a group of neoplasms with a characteristic appearance (blue staining) most often observed in children.

The SRBCT dataset contains the following variables:
  
```srbct$gene``` (continuous matrix): 63 rows and 2308 columns. The expression levels of the 2308 genes across the 63 tested subjects.

```srbct$class``` (categorical vector): class vector of length 63. Contains the tumour class of each individual. There are four classes which include Burkitt Lymphoma (BL), Ewing Sarcoma (EWS), neuroblastoma (NB) and rhabdomyosarcoma (RMS).

```srbct$gene.name``` (string vectors): two lists of length 2308. There is the Image.ID component containing an integer ID for each gene. The Gene.Description component briefly describes the nature of each gene.

To confirm the correct dataframe was extracted, the dimensions are checked. The distribution of class labels is also examined. It can be seen that these class labels are not balanced.

```
data(srbct) # extract the small round bull cell tumour data

X <- srbct$gene # use the gene expression data as the X matrix
Y <- srbct$class # use the class data as the Y matrix
```
```
dim(X) # check the dimensions of the X dataframe
```
```
## [1]   63 2308
```

```
summary(Y) # check the distribution of class labels
```
```
## EWS  BL  NB RMS 
##  23   8  12  20
```




As in most cases when generating models, exploring the data to determine the major sources of variation is a good first step. PCA is be used for this and centering and scaling is recommended to homogenize the variance across the genes. ncomp is set to an arbitrarily high number to understand the captured variance across cotheremponents.

```
# run pca method on data
pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE) 
plot(pca.srbct)  # barplot of the eigenvalues (explained variance per component)
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/4b7a72cd-d750-4964-b1ed-323f3e47c53d)


This barchart shows 2 components would be sufficient to explain a moderate proportion of the data's variance . Next, the data is projected onto these two components to attempt to observe sources of variation.

```
plotIndiv(pca.srbct, group = srbct$class, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') # onto the PCA subspace
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/31f988b9-2c5e-4220-b358-9b4a49be70b5)


It seems that different tumour types do not separate or cluster across the two Principal components of the data. There are clusters, but these are not explained by the class variable. It can be inferred then that the major source of variation is not attributed to tumour type, and the high variance within tumour types means that differential expression tests would be expected to return zero or very few significant genes. This also highlights one of the limitations of PCA as an unsupervised technique that maximizes variance instead of separation between class labels.

### Building the initial sPLS-DA model

```
srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later
```

In machine learning, an initial model is usually generated which is then fine tuned by optomizing the model parmeters. Here, a PLS-DA model is fitted with ten components to evaluate the performance and the number of components necessary for the final model. 

A sample dimension reduction plot, including confidence ellipses, is shown below. This plot shows much better clustering of samples according to the tumour type when compared to the PCA output.

```
# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2, 
          group = srbct$class, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/0e6c32af-3d1a-415e-bdca-a92bde2be235)


The background.predict() function can also be utilised to depict the separation of class labels as seen in the sample dimension reduction plot. This plot provides intuition on how new unseen samples would be classified according to the model generated by sPLS-DA.

```
# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = srbct$class, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")
          
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/b753c921-0771-4b87-b3dd-d455b4ca6052)


### Optomizing the parameters of sPLS-DA

### Selecting the number of components

### The ncomp Parameter

The number of components to use is a crucial decision and is dictated by the performance of the PLS-DA model – i.e. its ability to correctly classify novel samples. The perf() function is used for this exactly. This is done with repeated cross-validation. Based on the output of this function, the optimal number of components to use can be identified.

### Cross validation


Five-fold cross validation was repeated ten times using the BER of max.dist as the performance measure – where minimisation was optimal. For larger datasets with numerous samples, at least 10 folds is recommended. 3 or 5 folds is appropriate for smaller datasets and those with minimal samples should use Leave-One-Out (LOO) validation. The cpus parameter allows for the use of parallelisation of computation as tuning can take a long time on low-to-mid range processors.

The overall error rate (OER) and balanced error rate (BER) for the three different distance metrics (explained further below) across the first ten components are shown below

```
# undergo performance evaluation in order to tune the number of components to use
perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

```     
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/c14cfb04-34b3-4800-aff9-75d696298cd4)


From this, it seems three components are appropriate as the error for each distance metric decreases by very incremental amounts after this. Components beyond the third are likely to provide negligible returns to the classification accuracy. A more empirical way to select this number is through the $choice.ncomp component of the perf() output object. It runs t-tests for a significant different in mean error rate across components. Using the max.dist metric, this suggests that the optimal number of components is 4. When to use each distance metric is explained further below.

```
perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()
```
```
#        max.dist centroids.dist mahalanobis.dist
#overall        4              8                7
#BER            4              8                7
### Selecting the number of variables
```

### The keepX Parameter

In order to determine the number of variables used to construct each latent component, the tune.splsda() function is utilised. This is performed iteratively, such that components are tuned one at a time. Through this function, the classification error rate can be extracted and averaged across folds and repeats. For larger datasets, an appropriate number of repeats would be around 50-100.

```
# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 10, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 2) # allow for paralleliation to decrease runtime
```

Note: If this step is taking too long on your laptop, please download the pretrained model ```tune.splsda.srbct.RDS``` which is avaialable in Data folder.
This can be loaded intoo R using ```tune.splsda.srbct<-readRDS('tune.splsda.RDS')```

The output of the tuning is shown below. The diamond indicates the optimal number of variables to keep for a given component, selected by which keepX value achieves the lowest classification error rate as determined with a one-sided t−test. The error bars indicate the standard deviation across the repeated, cross-validated folds.

```
plot(tune.splsda.srbct, col = color.jet(4)) # plot output of variable number tuning
```

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/74516724-9cb3-4a5a-8d30-292df8f37e3d)


This also aids in further tuning the number of components. While the tuning of component number (through perf()) yielded an optimal value of 4, conflicting results can be seen after the use of tune.splsda(), such that the optimal value is claimed to be 3. After the optmisation of the selected features, the fourth component seemingly minimises the BER negligibly, and so 3 components are used in the final model.

```
tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()
## [1] 3
```

The exact quantity of features to use for each component can also be extracted from this object:
```
tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()
## comp1 comp2 comp3 comp4 
##     9   260    30    10
```

These values are stored to form the final, optimised model.
```
optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]
```

### Final Model

Using all the tuned parameters from above, the final sPLS-DA model can be formed.

```
# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)
```

### Plots

### Sample dimension reduction plots

Sample plots will be used to show the distribution of the data in the latent space based on the final sPLS-DA model. The plots seen below  show a much cleaner seperation between tumour classes than what was observed in the earlier unoptimised PLS-DA model (without feature or component selection). The plots below show the dimension reduction plots for the first and second components (a), and the first and third components (b). The difference between (a) and (b) is indicative of the fact that different genes discriminate the samples differently. Genes which contributed to the third component separated the RMS and NB classes much better than those which contributed to the second. All three components were well suited to separate the BL class as it does not overlap any other cluster in either plot.

(a) first and second components
```
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = srbct$class, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/8eb0dee0-8bac-48ee-80dc-f30b4a33759c)


(b) first and third components
```
plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = srbct$class, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')
```
![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/6e209324-9ea4-4e1d-9904-0d3bcdc5f548)


We then generate a heatmap to depict the expression levels of each gene (selected for component construction) for every sample. Euclidean distance with a complete agglomeration method were used to yield this heatmap. It can be seen that certain sets of genes selected by the sPLS-DA had homogeneous expression for different classes. For example, nearly half of the genes had high expression with the EWS (blue) tumour.

```
# set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Tumour Type", # legend title
            cex = 0.7) # legend size

# generate the heatmap, using the legend and colouring rows by each sample's class
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), 
           legend = legend)
```

### Feature importance plots

The stability of a given feature, or gene in this case, is defined as the proportion of cross validation folds (across repeats) where it was selected for to be used for a given component. Stability values can be extracted via perf.splsda.srbct$features$stable and ploted (see below). Those with the highest stability are likely to be much more “important” for a given component. The genes used for the first component had consistently lower stability than the other two. This can be explained as there are various combinations of genes that are discriminative on component 1, whereas the number of combinations decreases as component 2 is formed.

```
# form new perf() object which utilises the final model
perf.splsda.srbct <- perf(final.splsda, 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.srbct$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.srbct$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.srbct$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
     
```

Another feature plot to be used is the correlation circle plot shown below. By considering both the correlation circle plot and the sample plot, a group of genes with a positive correlation with component 1 (EH domain, proteasome etc.) are observed to be associated with the BL samples. Two groups of genes are either positively or negatively correlated with component 2. These genes are likely to characterize either the NB and RMS tumor classes, or the EWS tumor class respectively.

```
var.name.short <- substr(srbct$gene.name[, 2], 1, 40) # form simplified gene names truncated to 40 charachters

plotVar(final.splsda, comp = c(1,2), var.names = list(var.name.short), cex = 2.5) # generate correlation circle plot

```

Even though we truncated the gene names, some are still dificult to visualize with 2 components. We can use instead use 3 components and make an interactive 3D plot using the ```rgl``` package which allows us to rotate and zoom to better visualize the genes

```
install.packages("rgl")
library("rgl")

plotVar(final.splsda, comp = c(1,2,3), var.names = list(var.name.short), cex = 0.6, style = '3d', label.axes.box = 'axes') #3d plot
```

### Prediction of tumor classs for unseen samples

When undergoing prediction, the sPLS-DA data must first be segmented into training and testing, such that there are novel samples to evaluate performance on. Otherwise, it runs the risk of “overfitting”, resulting in inflated predictive ability scores.

```
train <- sample(1:nrow(X), 50) # randomly select 50 samples in training
test <- setdiff(1:nrow(X), train) # rest is part of the test set

# store matrices into training and test set:
X.train <- X[train, ]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]
```

A model is then trained on the training data. Note that the previously calculated optimal.keepX values are used here. In real scenarios, the training model should be tuned itself. It is crucial that when tuning the training model, it is done in the absence of the testing data. This also reduces likelihood of overfitting.

```
# train the model
train.splsda.srbct <- splsda(X.train, Y.train, ncomp = optimal.ncomp, keepX = optimal.keepX)
```

The model is then applied on the test set using a specific distance metric. In this case, the Mahalanobis distance was used (arbitrarily).

```
# use the model on the Xtest set
predict.splsda.srbct <- predict(train.splsda.srbct, X.test, 
                                dist = "mahalanobis.dist")
```

To evaluate the predictive performance, confusion matrices can be used. Directly below is such a matrix for a model using just the first two components. Only one misclassification were made – one sample was claimed to belong to the NB class when it truthfully belonged to the RMS class.

```
# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda.srbct$class$mahalanobis.dist[,2]

table(factor(predict.comp2, levels = levels(Y)), Y.test)

##      Y.test
##       EWS BL NB RMS
##   EWS   5  0  0   0
##   BL    0  1  0   0
##   NB    1  0  4   1
##   RMS   0  0  0   1

```
Below is the equivalent matrix for the model using all three components. It can be seen that the classification accuracy increased.

```
# evaluate the prediction accuracy for the first three components
predict.comp3 <- predict.splsda.srbct$class$mahalanobis.dist[,3]
table(factor(predict.comp3, levels = levels(Y)), Y.test)

##      Y.test
##       EWS BL NB RMS
##   EWS   6  0  0   0
##   BL    0  1  0   0
##   NB    0  0  4   0
##   RMS   0  0  0   2
```

### Performance Plots

AUC plots can be used for performance evaluation. AUC scores are calculated from training cross-validation sets and averaged in the ```perf()``` function (```perf.plsda.srbct$auc``` and ```perf.plsda.srbct$auc.all``` for one vs. one class or one vs. all classes respectively).

AUROC plots for models containing one component and three components can be seen below and suggest that the sPLS-DA model can distinguish BL subjects from the other groups with a high true positive and low false positive rate, while the model is less well able to distinguish samples from other classes on component 1. The model including all three components has a perfect classification accuracy.

Note that if print = TRUE (as is by default), numerical output including AUC and a Wilcoxon test p-value for each ‘one vs. other’ class comparisons that are performed per component will be printed.

```
auc.splsda = auroc(final.splsda, roc.comp = 1, print = FALSE) # AUROC for the first component
auc.splsda = auroc(final.splsda, roc.comp = 3, print = FALSE) # AUROC for all three components
```

## Tasks

Using the ```plotVar``` function and ```rgl``` can you identify the 8 genes with high positive correlation with NB and RMS tumors (component 2)?

Is there any way you can make an improved model with a higher AUC by tuning the parameters in a different way ? 

*(hint: in the help tab of R Studio type ```tune.splsda``` to find a list of parameter  that can be modified )*


