# R programming
`R` is a free, open-source programming language used for statistical computing and graphical presentation while `RStudio` is an integrated development environment (IDE) featureing tools for plotting, viewing history, debugging and managing your workspace for R and Python. 

In the era of data science, there are several pros and cons of using R programming in analyzing clinical or biological data:

<img src="https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/3d928712-a632-43e1-b7da-0bcdfc445d96" width=600 >


## 1. Installation
For this and the coming tutorials, you will need to install
- R: [https://cloud.r-project.org/](https://cloud.r-project.org)
- RStudio: [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)

After installation, launch RStudio (and, therefore R at the same time). You will see 4 default primary panes in the user interface:
- **Source** : where you document your source codes
- **Console** : where you run the codes interactively
- **Environments** : containing the `Environment`, `History`, `Connections`, `Tutorial` tabs, etc
- **Output** : containing the `Files`, `Plots`, `Packages`, `Help`, `Viewer`, and `Presentation` tabs

![image](https://github.com/BioinfoHKUSurgery/Bioinfo-Workshop-2024/assets/165180561/5e51192a-ef2a-44b4-9653-19ebfa423f23)

## 2. Data types
In R programming language, there are a few commonly used data types predefined in the built-in environment. In total, R has five basic data types:

1. Numeric
2. Integers
3. Logical
4. Characters
5. Complex

You can use the class() and typeof() functions to check the class and data type of any variable. 

```{R}
# Example
x <- "Surgery"
typeof(x)
```
