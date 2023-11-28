\name{Qlocalstat-package}
\alias{Qlocalstat-package}
\alias{Qlocalstat}
\docType{package}
\title{
\packageTitle{Qlocalstat}
}
\description{
\packageDescription{Qlocalstat}
}
\details{
% The DESCRIPTION file:
% \packageDESCRIPTION{Qlocalstat}
% \packageIndices{Qlocalstat}
The \code{Qstat} function in this package is the core function; it calculates a score statistic to test the same-population assumption in two-sample MR in a locus of variants based on GWAS summary statistics (estimate and standard error) for an exposure and outcome dataset and the LD matrix of the variants in the exposure data. This function also includes a variety of options for processing these inputs (in particular the LD matrix) to ensure numerical stability in calculation of this score statistic.

This package also includes two other functions - \code{PDcheck} and \code{eigen.pinv} - the former checks whether an input LD matrix is positive definite and suggests workarouns in the \code{Qstat} function if it is not; the latter is called within \code{Qstat} as part of the \code{"eigen"} argument for \code{SVD_thresh} to calculate an eigenvalue-based SVD pseudoinverse
}
\author{
\packageAuthor{Qlocalstat}

Maintainer: \packageMaintainer{Qlocalstat}
}
\references{
Zhu X, Stephens M. \href{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/issue-3/Bayesian-large-scale-multiple-regression-with-summary-statistics-from-genome/10.1214/17-AOAS1046.full}{BAYESIAN LARGE-SCALE MULTIPLE REGRESSION WITH SUMMARY STATISTICS FROM GENOME-WIDE ASSOCIATION STUDIES.} Ann Appl Stat. 2017;11(3):1561-1592. doi: 10.1214/17-aoas1046. Epub 2017 Oct 5. PMID: 29399241; PMCID: PMC5796536.

\href{https://www.cse.unr.edu/~bebis/MathMethods/SVD/lecture.pdf}{Truncated SVD approximation of matrix inverse}

\href{https://web.eecs.umich.edu/~aey/recent/regular.pdf}{Regularized Matrix Computations} (for Tikhonov regularization and truncated SVD)
}

\keyword{ package }
\seealso{
}
\examples{
  R <- matrix(c(1, 0.964356, 0.964356, 1), nrow = 2, ncol = 2)
  bx <- c(0.4724217, 0.4681516)
  by <- c(0.1840862, 0.1632529)
  
  Qstat(center = "mean", bx = bx, by = by,
                          se_bx = rep(0.05, 2), se_by = rep(0.05,2), R = R,
                          weak_filter = TRUE, weak_thresh = 2,
                          SVD = TRUE, SVD_thresh = "eigen", reg = 1e-16)
}
