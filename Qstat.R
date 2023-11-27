#' Locus-Based Test to Assess The Same-Population Assumption
#'
#' @param center Specify a center of choice from "index" "mean" "median"
#' "index" picks the most significant SNP as the representative
#' "mean" picks the mean of the ratios
#' "median" picks the median of the ratios
#' @param bx The beta of the exposure traits of the ancestry 1
#' @param by The beta of the outcome traits of the ancestry 2
#' @param se_bx The standard error of the beta of the exposure traits of the ancestry 1
#' @param se_by The standard error of the beta of the outcome traits of the ancestry 2
#' @param R The LD matrix
#' @param weak_filter filters out weak instruments below a certain squared z-score threshold
#' @param weak_thresh thresh of filter
#' @param SVD uses SVD to approximate the matrix inverse to get rid of the singular values below a certain threshold
#' @param SVD_thresh for inverse approximation by SVD
#' @param reg regularization parameter (if the threshold is set too high, your results will be biased!)
#' @return A list containing Qstat pvla (p value) and df (degree of fredom)
#' @examples 6
#' @importFrom stats median pchisq
#' @export
Qstat <- function(center = NULL, bx, by, se_bx, se_by, R,
                  weak_filter = FALSE, weak_thresh = NULL,
                  SVD = FALSE, SVD_thresh = NA, reg = 0) {
  ##An NA dropping step that is always performed - getting rid of any NA results
  ratios <- by/bx
  if (any(is.na(ratios))) {
    ix <- which(is.na(ratios))
    R <- R[-ix, -ix]
    ratios <- ratios[-ix]
    se_by <- se_by[-ix]
    se_bx <- se_bx[-ix]
    bx <- bx[-ix]
    by <- by[-ix]
  }

  ##filtering out weak instruments below a certain squared z-score threshold
  if (weak_filter) {
    absZ <- abs(bx/se_bx)
    ix <- which(absZ < weak_thresh)
    if (length(ix) != 0) { #if the length of ix is not non-zero, do this
      R <- R[-ix, -ix]
      ratios <- ratios[-ix]
      se_by <- se_by[-ix]
      se_bx <- se_bx[-ix]
      bx <- bx[-ix]
      by <- by[-ix]
    }
  }
  if (length(by) == 0) {
    print("You've removed everything!")
    return(list(Qstat = NA, pval = NA, df = NA))
  }


  ##The Q-statistic depends on the choice of center
  #while the user can specify this, we will also give them three simple choices: "index", "mean", "median"
  if (center == "index") { #picking the most significant SNP as the representative
    absZ <- abs(bx/se_bx)
    ind_var <- which.max(absZ) #the index variant
    center <- by[ind_var]/bx[ind_var]
  }

  if (center == "mean") { ##the mean of the ratios
    center <- mean(ratios)
  }

  if (center == "median") { ##the median of the ratios
    center <- median(ratios)
  }

  #Calculating Omega, and applying Tikhonov regularization if desired.
  #Beware of over-applying regularization: if the threshold is set too high, your results will be biased!
  Omega <- diag(se_bx) %*% R %*% diag(se_bx) + 1/center^2*diag(se_by) %*% R %*% diag(se_by)
  Omega <- Omega + reg*diag(dim(Omega)[1])

  if (SVD) {
    #using SVD to approximate the matrix inverse to get rid of the singular values below a certain threshold
    svdO <- svd(Omega)
    if (is.na(SVD_thresh)) {
      stop("filtering threshold needed for inverse approximation by SVD!")
    }
    dinv <- ifelse(svdO$d > SVD_thresh, 1/svdO$d, 0)
    invOmega <- svdO$v %*% diag(dinv) %*% t(svdO$u)
    Qstat <- as.numeric(t(by/center - bx) %*% invOmega %*% (by/center - bx))
    pval <- pchisq(Qstat, df = length(bx) - 1, lower.tail = FALSE)
    return(list(Qstat = Qstat, pval = pval, df = length(bx) - 1))
  } else {
    #no SVD approximation
    Qstat <- as.numeric(t(by/center - bx) %*% solve(Omega) %*% (by/center - bx))
    pval <- pchisq(Qstat, df = length(bx) - 1, lower.tail = FALSE)
    return(list(Qstat = Qstat, pval = pval, df = length(bx) - 1))
  }
}

