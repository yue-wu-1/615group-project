
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
```

or by downloading the tarball in the repository via:

```{r example}
install.packages("Qlocalstat_1.0.tar.gz",repos = NULL)
```
