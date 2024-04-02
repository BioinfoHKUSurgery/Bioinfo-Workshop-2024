# R programming
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
5. Factors

You can use the class() and typeof() functions to check the class and data type of any variable.

#### 2.1.1 Numeric (or double)
The `numeric` is for numeric values, as the most common and the default data type. 
```r
x <- c(1, 2, 4, 16)
x
#> [1]  1  2  4 16
class(x)
#> [1] "numeric"
typeof(x)
#> [1] "double"
```

#### 2.1.2 Integer
The `integer` is another data type used for the set of all integers. You can use the capital ‘L’ notation as a suffix to specify a particular value as the integer data type. Also, you can convert a value into an integer type using the as.integer() function.
```r
y <- c(1L, 2L, 4L, 16L)
y
#> [1] 1 2 4 16
class(y)
#> [1] "integer"
typeof(y)
#> [1] "integer"
```
```r
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

```r
2>1
z <- c(TRUE, TRUE, FALSE, FALSE)
z
#> [1]  TRUE  TRUE  FALSE FALSE
typeof(z)
#> [1] "logical"
```

#### 2.1.4 Character
The `character` is a data type where you have all the alphabets and special characters. It stores character values or strings. Strings in R can contain alphabets, numbers, and symbols. The character type is usually denoted by wrapping the value inside single or double inverted commas.
```r
w <- c("su", "rg", "er", "y")
w
#> [1] "su", "rg", "er", "y"
typeof(w)
#> [1] "character"
```

#### 2.1.5 Factor ***
`Factor` is a special case of `character` data type, which is often used to represent categorical data. Categories of each factor data type are known as `LEVELS`. E.g. for gender,

```r
gender <- factor(c("female", "male", "male", "female", "male"))
gender
#> [1] female male male   female male  
#> Levels: female male
```
To know the different levels of a factor variable, use `levels()`:
```r
levels(gender)
#> [1] "female" "male"
```

### 2.2 Data Structures
Data structure involves how the data is organised, accessed and modified. Using an appropriate data structure may largely improve computing efficiency.

#### 2.2.1 Vector
`Vector` is a basic data structure in R. It contains elements in the same data type (no matter double, integer, character, etc). You can check data type by using typeof() function and length of the vector by length() function.
```r
x
length(x)
#> [1] 4
```
```r
x <- rep(3, 5)    # repeat 3 for 5 times
x <- 1:12         # integer from 1 to 12
```

#### 2.2.2 Matrix
`Matrix` is a two-dimensional data structure. It is in principle built based on vector but has more convenient built-in functions for computation. It has rows and columns, both of which can also have names. To check the dimensions, you can use the `dim()` function.

```r
A <- matrix(1:12, nrow=3)
A
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    4    7   10
#> [2,]    2    5    8   11
#> [3,]    3    6    9   12

dim(A)
#> [1] 3 4

B <- matrix(1:12, nrow=3, byrow=TRUE)
B
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    2    3    4
#> [2,]    5    6    7    8
#> [3,]    9   10   11   12

colnames(A) <- c("C1","C2","C3","C4")
rownames(A) <- c("R1","R2","R3")
A
#>    C1 C2 C3 C4
#> R1  1  4  7 10
#> R2  2  5  8 11
#> R3  3  6  9 12
```

To index vector, you can use `logical` or `integer` (starting from 1), or the element name if it has. We can also use negative integers to return all elements except those specified. But we cannot mix positive and negative integers while indexing and real numbers if used are truncated to integers.

```r
x <- 1:12

x[3]
#> [1] 3

x[2:5]
#> [1] 2 3 4 5

x[c(2, 5, 6)]                   # index with integer
#> [1] 2 5 6

x[c(TRUE, FALSE, FALSE, TRUE)]  # index with logical value
#> [1]  1  4  5  8  9 12
```
```r
A
#>    C1 C2 C3 C4
#> R1  1  4  7 10
#> R2  2  5  8 11
#> R3  3  6  9 12

A[1, 2]
#> [1] 4

A[1, "C2"]
#> [1] 4

A[1, c(2, 3)]
#> C2 C3 
#>  4  7

A[1:2, c(2, 3)]
#>    C2 C3
#> R1  4  7
#> R2  5  8

A[-1, -1]
```
You can also modify values of an element of the vector/matrix by index.
```
1.2.2.2 Modify values
A[1, 2:4] <- c(-3, -5, 20)
A
#>    C1 C2 C3 C4
#> R1  1 -3 -5 20
#> R2  2  5  8 11
#> R3  3  6  9 12
```

#### 2.2.3 List
Different from `vector` that has all elements in the same data type, the list data structure can have components of mixed data types. We can use `str()` function to view the structure of a list (or any object).
```r
x <- list(2.5, TRUE, 1:3)
x
#> [[1]]
#> [1] 2.5
#> 
#> [[2]]
#> [1] TRUE
#> 
#> [[3]]
#> [1] 1 2 3

str(x)
#> List of 3
#>  $ : num 2.5
#>  $ : logi TRUE
#>  $ : int [1:3] 1 2 3
We can also have a name for each element:

x <- list("a" = 2.5, "b" = TRUE, "c" = 1:3)
```
Different from `vector` and `matrix`, for a `list`, you need to use double-layer square brackets, either by numeric index or name. Alternatively, you can also use `$` symbol with the name.
```r
x[[3]]
#> [1] 1 2 3

x[["c"]]
#> [1] 1 2 3

x$c
#> [1] 1 2 3
```

#### 2.2.4 Data Frame
`Data frame` is widely used for rectangular data, where each column has the same data type (`vector`) but different columns can have different data types (like Excel). It can be treated as a special type of list: A list of vectors with the **same length**.

```r
df <- data.frame("SampleID" = 1:5, "Age" = 15:19, 
                 "Name" = c("John","Peter","Paul","Mary","Harry"))
df
#   SampleID Age  Name
# 1        1  15  John
# 2        2  16 Peter
# 3        3  17  Paul
# 4        4  18  Mary
# 5        5  19 Harry
```

## 3. Data management
### 3.1 Data input : read in .csv table

## 4. Statistical analysis

