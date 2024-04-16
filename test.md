### "Introduction to R: data manipulation and visualization"

### 1. Data manipulation

R's rich ecosystem of packages, statistical analysis capabilities, reproducibility features, and active community make it a powerful tool for data manipulation and visualization. Its versatility and flexibility make it suitable for a wide range of data analysis tasks, from exploratory data analysis to advanced statistical modeling.

R offers a wide range of packages and functions for data manipulation, transformation, and cleaning. Apart from the multitude of base-r functions such as `apply` and loop constructs, The `dplyr` and `tidyverse` packages also provide a powerful set of tools for filtering, sorting, aggregating, and joining datasets. These packages make it easy to perform complex data manipulations efficiently.

### 1. For Loops

In R, you can use a `for` loop to iterate over a sequence of values or elements and perform a set of operations repeatedly.

```{r}

# Iterate over a sequence of numbers
for (i in 1:5) { 
  # Print the current value of i 
  print(i) 
} 

#>[1] 1
#>[1] 2 
#>[1] 3 
#>[1] 4 
#>[1] 5

```

You can also use a for loop to iterate over elements of a vector, list, or any other iterable object. Here's an example:

```{r}

# Iterate over a vector samples

samples <- c("control-1", "control-2", "treatment-1", "treatment-2")
for (s in samples) { 
   #print each sample name
   print(s)
}

#>[1] "control-1" 
#>[1] "control-2"
#>[1] "treatment-1"
#>[1] "treatment-2"


```

#### 1.2. Nested For loops

To illustrate a nested loop example on clinical data using a data frame in R, let's assume you have a dataset with clinical information about patients, including their ID, blood pressure, and cholesterol levels. You want to perform a specific operation on each patient's data using a nested loop structure. Here's an example:

```{r}

# Sample clinical data clinical_data

clinical_data <- data.frame(id = c(1, 2, 3),
                            blood_pressure = c(120, 130, 140),
                            cholesterol = c(180, 200, 220) )

print(clinical_data)

#>  id blood_pressure cholesterol
#>  1            120         180
#>  2            130         200
#>  3            140         220


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

#> Patient ID: 1 #\>blood_pressure : 120 #\>cholesterol : 180
#> Patient ID: 2 #\>blood_pressure : 130 #\>cholesterol : 200
#> Patient ID: 3 #\>blood_pressure : 140 #\>cholesterol : 220
```

#### 1.2. Conditional IF/ELSE statements

Next we will use IF/ELSE statements to check if any patient has high blood pressure and to prescibe a specific medication in the case of high BP.

```{r}

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


#> Patient ID: 1
#> Cholesterol level: 180

#> Patient ID: 2 
#> Cholesterol level: 200

#> Patient ID: 3 #\>Cholesterol level: 220
#> High blood pressure: 140
#> Prescription: Zestril 

```

#### 1.3. Apply() functions

The `apply()` function in R allows you to apply a specified function to each row or column of a data frame or matrix. It usually runs much faster than for loops and is more concise to code, but may limited flexibility when dealing with more complicated data manipulation procedures.

Here we will apply the R mean() function to calculate the average cholesterol and blood pressure for the patients

```{r}

# Apply function to calculate the mean for each column column_means

column_means <- apply(clinical_data, 2, mean)
print(column_means)

#> id blood_pressure cholesterol 
#> 2 130 200
```

#### 2. Data Visualization in R

One of the main advantages of R is the extensive data visualization functionality. ggplot2 is the most widely used package and follows a layered approach, where you can add different layers to a plot to build up the desired visualization, making it very flexible and easy to modify plots

#### 2.1 Basic scatterplot in ggplo2

```{r}
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
```

```{r echo=FALSE}
ggplot(data = df, aes(x = age, y = BP)) + 
  geom_point(size=4, col='blue') + 
  labs(x = "age", y = "BP") + 
  ggtitle("Scatter Plot: age vs BP")
```
