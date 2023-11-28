\name{PDcheck}
\alias{PDcheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checking Positive Definiteness of the Score Statistic Matrix
}
\description{
Calculates the \code{Omega} of the score statistic based on an LD matrix and summary statistics, checks whether this matrix is positive definiteness, then, if the matrix is not positive definite, prints recommendations on which parameters to use for the truncated SVD pseudoinverse.
}
\usage{
PDcheck(R, se_bx, se_by, center)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{The LD matrix of the locus variants in the ancestry of the exposure data - make sure it's symmetric!}
  \item{se_bx}{The standard errors of the exposure GWAS effect estimates.}
  \item{se_by}{The standard errors of the outcome GWAS effect estimates.}
  \item{center}{The choice of center for the score statistic - unlike the "center" of \code{Qstat}, this only accepts numbers.}
}
\details{
This method calculates the eigenvalues of the matrix \code{Omega} - if they are not all positive (i.e. the matrix is not PD), then it determines which eigenvalues are negative and prints messages describing how to handle these negative eigenvalues when using \code{Qstat} (either a range of good SVD thresholds if using a threshold-based pseudoinverse, or whether to use the eigenvalue-based pseudoinverse).
}
\value{
The recommendations and check for positive definiteness are printed as messages - the function either returns 0 if the matrix is positive definite, or the negative eigenvalue with the largest magnitude if the matrix is not positive definite.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Li, Jack; Wu, Yue; Wang; Zhuoyu
}
\note{
Should be run on the LD matrix and summary statistics ahead of time before running \code{Qstat}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Qstat}}
}
\examples{
R <- matrix(c(1, 0.990459, 0.875089, 0.882848,
                    0.990459, 1, 0.884442, 0.874105,
                    0.875089, 0.884442, 1, 0.991415,
                    0.882848, 0.874105, 0.991415, 1), nrow = 4, ncol = 4)

center <- 0.1809437

PDcheck(R = R, se_bx = rep(0.05, 4), se_by = rep(0.05, 4), center = center)
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
