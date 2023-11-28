\name{Qstat}
\alias{Qstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Score Statistic to Assess the Same-Population Assumption}
\description{
Calculates a score statistic to test the same-population assumption in two-sample MR in a locus of variants based on GWAS summary statistics (estimate and standard error) for an exposure and outcome dataset and the LD matrix of the variants in the exposure data.
}
\usage{
Qstat(center = NULL, bx, by, se_bx, se_by, R,
weak_filter = FALSE, weak_thresh = NULL,
SVD = FALSE, SVD_thresh = NA, reg = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{center}{The choice of center for the score statistic - this can be a user specified number, or one of three characters for preset options: \code{"index"} for the Wald ratio of the index variant (the variant that is most significantly associated with the exposure), \code{"mean"} for the mean Wald ratio of the variants in the locus, or \code{"median"} for the median Wald ratio of variants in the locus.}
  \item{bx}{The GWAS effect estimates for the exposure data in ancestry 1.}
  \item{by}{The GWAS effect estimates for the outcome data in ancestry 2.}
  \item{se_bx}{The standard errors of the exposure GWAS effect estimates.}
  \item{se_by}{The standard errors of the outcome GWAS effect estimates.}
  \item{R}{The LD matrix of the locus variants in the ancestry of the exposure data. Make sure ahead of time that this matrix will result in a SPD result using \code{PDcheck}, or instability may result!}
  \item{weak_filter}{A logical indicating whether to filter out weak instruments in the locus based on the magnitude of the z-score of their associations with the exposure.}
  \item{weak_thresh}{The cutoff for z-score magnitude if filtering out weak instruments.}
  \item{SVD}{A logical indicating whether or not to approximate the matrix inverse in the score statistic using the truncated SVD pseudoinverse.}
  \item{SVD_thresh}{Either a user-specified number indicating the cutoff to filter out small singular values for the threshold-based pseudoinverse, or the character "eigen" to calculate the eigenvalue-based pseudoinverse.}
  \item{reg}{A positive number indicating the constant for use for Tikhonov regularization - setting this parameter too large may bias your results!}
}
\details{
The LD matrix should be symmetric, and the dimensions of the LD matrix and the lengths of each of the summary statistics vectors should match with each other.
}
\value{
Returns a list with three components
\item{Qstat}{The score statistic for a locus}
\item{pval}{The corresponding p-value for the score statistic}
\item{df}{The degrees of freedom for the ch-squared distribution of the score statistic under the null hypothesis}
}
\references{
Zhu X, Stephens M. \href{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/issue-3/Bayesian-large-scale-multiple-regression-with-summary-statistics-from-genome/10.1214/17-AOAS1046.full}{BAYESIAN LARGE-SCALE MULTIPLE REGRESSION WITH SUMMARY STATISTICS FROM GENOME-WIDE ASSOCIATION STUDIES.} Ann Appl Stat. 2017;11(3):1561-1592. doi: 10.1214/17-aoas1046. Epub 2017 Oct 5. PMID: 29399241; PMCID: PMC5796536.

\href{https://www.cse.unr.edu/~bebis/MathMethods/SVD/lecture.pdf}{Truncated SVD approximation of matrix inverse}

\href{https://web.eecs.umich.edu/~aey/recent/regular.pdf}{Regularized Matrix Computations} (for Tikhonov regularization and truncated SVD)
}
\author{
Li, Jack; Wu, Yue; Wang; Zhuoyu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{PDcheck}}, \code{\link{eigen.pinv}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
