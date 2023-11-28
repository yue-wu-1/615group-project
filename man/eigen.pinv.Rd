\name{eigen.pinv}
\alias{eigen.pinv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Eigenvalue-based Pseudoinverse}
\description{
Calculates the eigenvalue-based pseudoinverse of a matrix, \code{Omega}.
}
\usage{
eigen.pinv(Omega)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Omega}{The matrix to calculate the pseudoinverse - this should be a symmetric matrix.}
}
\details{
  The eigenvalue-based pseudoinverse relies on the fact that correlation matrices are SPD (meaning that all eigenvalues should be positive), while estimates of correlation matrices may not be SPD, but are still guaranteed to be symmetric (meaning that the absolute values of eigenvalues and singular values are the same, but may differ in sign). This pseudoinverse removes the contribution of singular values that correspond to negative eigenvalues (which are normally impossible), which are highly likely to be inaccurately calculated.
}
\value{
Returns a matrix with the same dimensions as \code{Omega}.
}

\author{
Li, Jack; Wu, Yue; Wang; Zhuoyu
}
\note{
This function is called within the \code{Qstat} function; there isn't usually a reason to run it on its own.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
  Omega <- matrix(c(1, 0.990459, 0.875089, 0.882848,
                    0.990459, 1, 0.884442, 0.874105,
                    0.875089, 0.884442, 1, 0.991415,
                    0.882848, 0.874105, 0.991415, 1), nrow = 4, ncol = 4)
  
  eigen.pinv(Omega)
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
