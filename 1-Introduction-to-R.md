# R programming
`R` is a free, open-source programming language used for statistical computing and graphical presentation while `RStudio` is an integrated development environment (IDE) featureing tools for plotting, viewing history, debugging and managing your workspace for R and Python. 

In the era of data science, there are several pros and cons of using R programming in analyzing clinical or biological data:

<img src="https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/3d928712-a632-43e1-b7da-0bcdfc445d96" width=600 >


## 1. Installation
For this and the coming tutorials, you will need to install
- R: [https://cloud.r-project.org/](https://cloud.r-project.org)
- RStudio: [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)

After installation, launch RStudio (and, therefore R at the same time). You will see 4 default primary panes in the user interface:
- **Source** : where you view or edit your source codes (such as `.R`, `.rmd`)
- **Console** : where you run the codes interactively
- **Environments** : containing the `Environment`, `History`, `Connections`, `Tutorial` tabs, etc
- **Output** : containing the `Files`, `Plots`, `Packages`, `Help`, `Viewer`, and `Presentation` tabs

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/5e51192a-ef2a-44b4-9653-19ebfa423f23)


## 2. Very basics of R
The workspace is your current R working environment which includes the user-defined objects (e.g. vectors, data frames, lists, and functions). At the end of an R session, user can save an image of the current workspace that can be reloaded the next time when R is started under the `Environment` of the Environments pane (top right). Here are some standard commands for managing your workspace.

```
getwd() # print the current working directory - cwd
setwd("C:/Users/Dr. Clara Tang/Documents")  # change to targeted directory
ls()    # list the objects in the current workspace
```

### 2.1 Data Types
Generally for all programming lanuages, we use variables to store data/information. Variables are reserved memory locations to store values. Different types of variables have different memory usage, which may affect the speed of data processing. 

In R, a variable is referred to an `object`. There are a few commonly used data types predefined in the built-in environment, including
1. Numeric
2. Integers
3. Logical
4. Characters

You can use the class() and typeof() functions to check the class and data type of any variable.

#### 2.1.1 Numeric (or double)
The `numeric` is for numeric values, as the most common and the default data type. 
```
x <- c(1, 2, 4, 16)
x
#>[1]  1  2  4 16
class(x)
#> [1] "numeric"
typeof(x)
#> [1] "double"
```

#### 2.1.2 Integer
The `integer` is another data type used for the set of all integers. You can use the capital ‘L’ notation as a suffix to specify a particular value as the integer data type. Also, you can convert a value into an integer type using the as.integer() function.
```
y <- c(1L, 2L, 4L, 16L)
y
#> [1] 1 2 4 16
class(y)
#> [1] "integer"
typeof(y)
#> [1] "integer"
```
```
# Assign a integer value to y
y <- 20
y_in <- as.integer(20)
 
# is y an integer?
print(is.integer(y))
#> [1] FALSE

# is y_in an integer?
print(is.integer(y_in))
#> [1] TRUE
```

#### 2.1.3 Logical
The `logical` data type takes either a value of `true` or `false`. A logical value is often generated when comparing variables.

```
2>1
z <- c(TRUE, TRUE, FALSE, FALSE)
z
#> [1]  TRUE  TRUE  FALSE FALSE
typeof(z)
#> [1] "logical"
```

#### 2.1.4 Character
The `character` is a data type where you have all the alphabets and special characters. It stores character values or strings. Strings in R can contain alphabets, numbers, and symbols. The character type is usually denoted by wrapping the value inside single or double inverted commas.
```
w <- c("su", "rg", "er", "y")
w
#> [1] "su", "rg", "er", "y"
typeof(w)
#> [1] "character"
```

### 2.2 Data Structures

## 3. Data management
### 3.1 Data input : read in .csv table

## 4. Statistical analysis

