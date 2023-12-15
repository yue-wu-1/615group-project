\name{matrix_thin}
\alias{matrix_thin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Thinning an LD matrix by removing perfectly (or nearly perfectly) correlated entries}
\description{
When the sample size for calculating a correlation matrix is small, there will occasionally be perfectly (or nearly perfectly) correlated variants in an LD matrix, represented by entries that are exactly or nearly exactly 1 on the off-diagonals of the matrix. Inclusion of too many of these perfectly correlated variants can destabilize the score statistic. To address this, this function "thins" out such a matrix by removing these perfectly correlated variants.
}
\usage{
matrix_thin(matrix, R2 = 0.99999, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{The LD matrix to be thinned}
  \item{R2}{our choice of correlation threshold for thinning. The default argument is set to 0.99999 to allow for some error in numeric precision in the calculation}
  \item{verbose}{whether to print information about what indices were removed (if any)}
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
  cor_mat <- matrix(c(1, 0.95, 1, 0.95, 1, 0.95, 1, 0.95, 1), nrow = 3)
  matrix_thin(cor_mat)
}
