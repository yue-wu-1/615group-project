
# Qlocalstat - A Package for Assessing the Same Population Assumption in Two-Sample MR using Linkage Disequilibrium Data

```{r}
library(Qlocalstat)
```

<!-- badges: start -->
<!-- badges: end -->

Qlocalstat is a score-statistic based method that takes GWAS summary statistics (estimate and standard error) for two traits (exposure, x and outcome, y) and an LD matrix for the exposure for a set of nearby variants (locus), then calculates a score statistic that determines if the heterogeneity of summary statistics in the locus is higher than expected (which is reflective of population mismatch). 

## Installation

You can install Qlocalstat from Github via:

``` r
devtools::install_github("yue-wu-1/615group-project")
library(Qlocalstat)
```

or by downloading the tarball in the repository and installing in R using the command:

```{r example}
install.packages("Qlocalstat_1.0.tar.gz",repos = NULL)
```

This package contains four functions: `Qstat`, the main function that calculates the score statistic and compares it to the appropriate $\chi^2$ null distribution, `PDcheck`, which assesses whether a correlation matrix is positive definite, `eigen.pinv`, a function called by `Qstat` that calculates the eigenvalue-based pseudoinverse if it is asked for, and `matrix_thin`, a function that thins the LD matrix ahead of time for better stability. To see more details about these functions, use `help().`

```{r example}
help(Qstat)
help(PDcheck)
help(eigen.pinv)
help(matrix_thin)
```
