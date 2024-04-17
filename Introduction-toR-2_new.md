### “Introduction to R: data manipulation and visualization”

### 1. Data manipulation and control structures

R’s rich ecosystem of packages, statistical analysis capabilities,
reproducibility features, and active community make it a powerful tool
for data manipulation and visualization. Its versatility and flexibility
make it suitable for a wide range of data analysis tasks, from
exploratory data analysis to advanced statistical modeling.

R offers a wide range of packages and functions for data manipulation,
transformation, and cleaning. Apart from the multitude of base-r
functions such as `apply` and loop constructs, The `dplyr` and
`tidyverse` packages also provide a powerful set of tools for filtering,
sorting, aggregating, and joining datasets. These packages make it easy
to perform complex data manipulations efficiently.

### 1. For Loops

In R, you can use a `for` loop to iterate over a sequence of values or
elements and perform a set of operations repeatedly.

    # Iterate over a sequence of numbers
    for (i in 1:5) { 
      # Print the current value of i 
      print(i) 
    } 

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5


You can also use a for loop to iterate over elements of a vector, list,
or any other iterable object. Here’s an example:

    # Iterate over a vector samples

    samples <- c("control-1", "control-2", "treatment-1", "treatment-2")
    for (s in samples) { 
       #print each sample name
       print(s)
    }

    ## [1] "control-1"
    ## [1] "control-2"
    ## [1] "treatment-1"
    ## [1] "treatment-2"

#### 1.1. Nested For loops

To illustrate a nested loop example on clinical data using a data frame
in R, let’s assume you have a dataset with clinical information about
patients, including their ID, blood pressure, and cholesterol levels.
You want to perform a specific operation on each patient’s data using a
nested loop structure. Here’s an example:

    # Sample clinical data

    clinical_data <- data.frame(id = c(1, 2, 3),
                                blood_pressure = c(120, 130, 140),
                                cholesterol = c(180, 200, 220) )

    print(clinical_data)

    ##   id blood_pressure cholesterol
    ## 1  1            120         180
    ## 2  2            130         200
    ## 3  3            140         220

    #Nested loop to iterate over each patient's data

    for (i in 1:nrow(clinical_data)) {
      # Print patient ID
      cat("Patient ID:", clinical_data$id[i], "\n")
        # Inner loop for specific operations on each patient's data
        for (j in 2:ncol(clinical_data)) {
          # Perform a specific operation (e.g., print blood pressure and cholesterol)
          cat(colnames(clinical_data)[j], ":", clinical_data[i, j], "\n")
        }
      cat("\n")
    }

    ## Patient ID: 1 
    ## blood_pressure : 120 
    ## cholesterol : 180 
    ## 
    ## Patient ID: 2 
    ## blood_pressure : 130 
    ## cholesterol : 200 
    ## 
    ## Patient ID: 3 
    ## blood_pressure : 140 
    ## cholesterol : 220


#### 1.2. Conditional IF/ELSE statements

Next we will use IF/ELSE statements to check if any patient has high
blood pressure and to prescibe a specific medication in the case of high
BP.

    # For loop with a conditional statement to perform a specific operation on each patient's data
    for (i in 1:nrow(clinical_data)) {
      
      # Check if blood pressure is above 
      if (clinical_data$blood_pressure[i] > 130) {
        
        # If blood pressure is above 130, print patient ID and blood pressure level 
        cat("Patient ID:",     clinical_data$id[i], "\n")
        cat("Cholesterol level:", clinical_data$cholesterol[i], "\n") 
        cat("High blood pressure:", clinical_data$blood_pressure[i], "\n")
        cat("Prescription:", "Zestril", "\n")
      
        } else {
        
          # If blood pressure is not above 130, print patient ID and cholesterol level cat("Patient ID:", clinical_data$id[i], "\n")
        cat("Patient ID:", clinical_data$id[i], "\n")
        cat("Cholesterol level:", clinical_data$cholesterol[i], "\n")
        
        }
      
      cat("\n")
    }

    ## Patient ID: 1 
    ## Cholesterol level: 180 
    ## 
    ## Patient ID: 2 
    ## Cholesterol level: 200 
    ## 
    ## Patient ID: 3 
    ## Cholesterol level: 220 
    ## High blood pressure: 140 
    ## Prescription: Zestril


#### 1.3. Apply() functions

The `apply()` function in R allows you to apply a specified function to
each row or column of a data frame or matrix. It usually runs much
faster than for loops and is more concise to code, but may limited
flexibility when dealing with more complicated data manipulation
procedures.

Here we will apply the R mean() function to calculate the average
cholesterol and blood pressure for the patients

    # Apply function to calculate the mean for each column column_means

    column_means <- apply(clinical_data, 2, mean)
    print(column_means)

    ##             id blood_pressure    cholesterol 
    ##              2            130            200

    #> id blood_pressure cholesterol 
    #> 2 130 200

#### 2. Data Visualization in R

One of the main advantages of R is the extensive data visualization
functionality. ggplot2 is the most widely used package and follows a
layered approach, where you can add different layers to a plot to build
up the desired visualization, making it very flexible and easy to modify
plots

#### 2.1 Basic scatterplot in ggplo2

    #example dataframe showing blood pressure increases with patient age
    df <- data.frame( age = c(45, 56, 34, 67, 51, 42, 59, 38, 48),
                      BP = c(120, 130, 125, 145, 135, 115, 140, 130, 125))
                
    #intall and load ggplot2
    #install.packages('ggplot2')
    library('ggplot2')

    #generate scatterplot
    ggplot(data = df, aes(x = age, y = BP)) + 
      geom_point(size=4, col='blue') + 
      labs(x = "age", y = "BP") + 
      ggtitle("Scatter Plot: age vs BP")

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/494321e7-126f-40ab-8e6c-ded7889cd490)


##### 2.2 Basic barplot in ggplo2

For barplots use the `ggplot()` function to set up the plot and the
`geom_bar()` function to add the bars.


    # Create a dataframe

    df <- data.frame( condition = c("Condition A", "Condition B", "Condition C"),
                      count = c(20, 30, 15))

    # Create the bar plot
    ggplot(data = df, aes(x = condition, y = count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(x = "Condition", y = "Count") +
      ggtitle("Bar Plot")

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/3376f6bd-d66e-4dfc-8cfd-ab06aeb2de94)

\#####. 2.3. Basic boxplots

ggplot2 can also be used to generate boxplots. In this example, we
create a dataframe named `df` with two columns: `group` and `value`. The
`group` column represents the different groups in the clinical study
(e.g., treatment groups or control groups), and the `value` column
represents the corresponding measurements or observations.

    # Create a dataframe

    df <- data.frame( group = rep(c("Group A", "Group B", "Group C"), each = 50),
                      value = c(rnorm(50, mean = 10, sd = 2), 
                      rnorm(50, mean = 15, sd = 3), rnorm(50, mean = 12, sd = 1.5)) )

    #generate boxplots
    ggplot(data = df, aes(x = group, y = value)) + 
      geom_boxplot(fill = "steelblue", color = "black") +
      labs(x = "Group", y = "Value") +
      ggtitle("Boxplot of treatment groups")

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/2c4de8b3-dd38-4182-be56-141e0a195ae2)

The `geom_boxplot()` function adds the boxplots to the plot. By default,
the boxplots show the median, interquartile range, and whiskers. You can
customize the appearance of the boxplots by using additional arguments
in `geom_boxplot()`, such as `fill` for the color of the boxes and
`color` for the color of the outlines.

\#####. 2.4. Adding trendlines

Using geom\_smooth() it is relatively straightforward to add trendlines
to line plot. Here a linear model fit is used, but others such as local
polynormal regression fitting (loess) can be used.

    patients <- c("Patient 1", "Patient 2", "Patient 3", "Patient 4", "Patient 5")
    values1 <- c(10, 25, 30, 40, 50)
    values2 <- c(20, 39, 25, 43, 60)

    #Create a scatter plot with trendline using ggplot2 

    data <- data.frame(patients, values1, values2)

    ggplot(data, aes(x = values1, y = values2)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("TPM Gene A") +
    ylab("TPM Gene B") +
    ggtitle("Gene Expression Scatter Plot with Trendline")

    ## `geom_smooth()` using formula = 'y ~ x'

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/f8cb4626-5937-48b8-9218-4c186a5cf869)


We can also add standard error layer to the graph

    ggplot(data, aes(x = values1, y = values2)) +
      geom_point() +
      geom_smooth(method = "lm", se=TRUE) +
      xlab("TPM Gene A") +
      ylab("TPM Gene B") +
      ggtitle("Gene Expression Scatter Plot with Trendline and SE")

    ## `geom_smooth()` using formula = 'y ~ x'

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/9c8711b3-4ccd-440f-9864-2c673a468086)

\#####. 2.5. Advanced Boxplots

scale\_fill\_brewer() from the ColorBrewer package contains predefined
color schemes which can be used to colour the boxplots plots. The
theme() function is used for changing plot aesthetics such as gridlines
and background with predefined styles. In the case below the
‘minimal\_theme’ was used

    library('reshape2')
    library('readr')
    library('dplyr')

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    num_genes <- 100
    num_samples <- 10

    #Generate random gene expression data
    gene_expression <- matrix(rnorm(num_genes * num_samples), nrow = num_genes, ncol = num_samples)

    #Convert gene expression data to a data frame
    gene_expression_df <- as.data.frame(gene_expression)
    gene_expression_melted <- melt(gene_expression_df)

    ## No id variables; using all as measure variables

    #Create colored boxplots using ggplot2

    boxplot_colored <- ggplot(gene_expression_melted,
    aes(x = variable, y = value, fill=variable)) +
    geom_boxplot() +
    scale_fill_brewer(palette = 'RdPu', n=10) +
    labs(x = "Samples", y = "Gene Expression", title = "Colored Boxplots of Gene Expression Data") +
    theme_minimal()

    boxplot_colored

    ## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette RdPu is 9
    ## Returning the palette you asked for with that many colors

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/039ff0ea-dde1-40b6-83a5-113829aeb13d)

\#####. 2.6. Advanced plotting: Dot plots for scRNA-seq data.

Gene expression dotplots are often used to visualize Single cell RNA-seq
data. Colour intensity = log2(count+1) Dot size = proportion of cells
expressing gene

Apart from ggplot2, other packages are used to generate the dotplot:
Cowplot for improved aesthetics Viridis for colour schemes dplyr for
modifying columns

    library('dplyr')
    library('cowplot')

    gene_cluster <- read_tsv('https://github.com/davemcg/davemcg.github.io/raw/master/content/post/scRNA_dotplot_data.tsv.gz')

    ## Rows: 888 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): Gene, cluster, Group
    ## dbl (3): cell_ct, cell_exp_ct, count
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    markers <- gene_cluster$Gene %>% unique()

    gene_cluster %>% filter(Gene %in% markers) %>%
    mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>%
    filter(count > 0, `% Expressing` > 1) %>%
    ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) +
    geom_point() +
    cowplot::theme_cowplot() +
    theme(axis.line = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') +
    theme(axis.ticks = element_blank()) +
    scale_color_gradientn(colours = viridis:: viridis (20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 12')

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165875740/5e6180df-d00f-43e3-93e1-a0f2ca9e7901)
