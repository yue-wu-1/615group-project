#' A Score Statistic to Assess the Same-Population Assumption
#'
#' @param center The choice of center for the score statistic - this can be a user specified number, or one of three characters for preset options: \code{"index"} for the Wald ratio of the index variant (the variant that is most significantly associated with the exposure), \code{"mean"} for the mean Wald ratio of the variants in the locus, or \code{"median"} for the median Wald ratio of variants in the locus.
#' @param bx The GWAS effect estimates for the exposure data in ancestry 1.
#' @param by The GWAS effect estimates for the outcome data in ancestry 2.
#' @param se_bx The standard errors of the exposure GWAS effect estimates.
#' @param se_by The standard errors of the outcome GWAS effect estimates.
#' @param R The LD matrix of the locus variants in the ancestry of the exposure data. Make sure ahead of time that this matrix will result in a SPD result using \code{PDcheck}, or instability may result!
#' @param weak_filter A logical indicating whether to filter out weak instruments in the locus based on the magnitude of the z-score of their associations with the exposure.
#' @param weak_thresh The cutoff for z-score magnitude if filtering out weak instruments.
#' @param SVD A logical indicating whether or not to approximate the matrix inverse in the score statistic using the truncated SVD pseudoinverse.
#' @param SVD_thresh Either a user-specified number indicating the cutoff to filter out small singular values for the threshold-based pseudoinverse, or the character "eigen" to calculate the eigenvalue-based pseudoinverse.
#' @param reg A positive number indicating the constant for use for Tikhonov regularization - setting this parameter too large may bias your results!
#' @return Returns a list with three components
#' \item{Qstat}{The score statistic for a locus}
#' \item{pval}{The corresponding p-value for the score statistic}
#' \item{df}{The degrees of freedom for the ch-squared distribution of the score statistic under the null hypothesis}
#' @examples   R <- matrix(c(1, 0.964356, 0.964356, 1), nrow = 2, ncol = 2)
#' bx <- c(0.4724217, 0.4681516)
#' by <- c(0.1840862, 0.1632529)
#'
#' Qstat(center = "mean", bx = bx, by = by,
#'       se_bx = rep(0.05, 2), se_by = rep(0.05,2), R = R,
#'       weak_filter = TRUE, weak_thresh = 2,
#'       SVD = TRUE, SVD_thresh = "eigen", reg = 1e-16)
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

  if (length(by) <= 1) {
    print("Not enough variants remain!")
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

    if (is.numeric(SVD_thresh)) {
      if (min(svdO$d) > SVD_thresh) {
        message(paste("NOTE: All singular values are above the selected threshold of", SVD_thresh))
        message(paste("The smallest singular value is", min(svdO$d)))
      }
      dinv <- ifelse(svdO$d > SVD_thresh, 1/svdO$d, 0)
      invOmega <- svdO$v %*% diag(dinv) %*% t(svdO$u)
      Qstat <- as.numeric(t(by/center - bx) %*% invOmega %*% (by/center - bx))
      if (Qstat < 0) warning("Q-statistic is less than zero, consider using a stricter SVD threshold!")
      pval <- pchisq(Qstat, df = length(bx) - 1, lower.tail = FALSE)
      return(list(Qstat = Qstat, pval = pval, df = length(bx) - 1))
    }

    if (SVD_thresh == "eigen") {
      #removing singular values that have a negative eigenvalue
      invOmega <- eigen.pinv(Omega)
      Qstat <- as.numeric(t(by/center - bx) %*% invOmega %*% (by/center - bx))
      pval <- pchisq(Qstat, df = length(bx) - 1, lower.tail = FALSE)
      return(list(Qstat = Qstat, pval = pval, df = length(bx) - 1))
    } else {
      stop("Invalid choice of SVD_thresh - either choose a number or use 'eigen'!")
    }
  } else {
    #no SVD approximation
    Qstat <- as.numeric(t(by/center - bx) %*% solve(Omega) %*% (by/center - bx))
    if (Qstat < 0) warning("Q-statistic is less than zero, consider using the pseudoinverse!")
    pval <- pchisq(Qstat, df = length(bx) - 1, lower.tail = FALSE)
    return(list(Qstat = Qstat, pval = pval, df = length(bx) - 1))
  }
}

#' @export
PDcheck <- function(R, se_bx, se_by, center) {
  #checks if the Omega generated by a certain choice of center of LD matrix, standard errors, and center is positive definite by checking eigenvalues
  #If it is not PD, points out which singular values have corresponding negative eigenvalues
  #then recommends a course of action for the pseudoinverse based on this.
  Omega <- diag(se_bx) %*% R %*% diag(se_bx) + 1/center^2*diag(se_by) %*% R %*% diag(se_by)
  eigenO <- eigen(Omega)

  if (any(eigenO$values < 0)) {
    message(cat("This matrix is not positive definite - negative eigenvalue indices include",
                which(eigenO$values < 0),
                "with eigenvalues",
                eigenO$values[which(eigenO$values < 0)], sep = " "))

    if (max(abs(eigenO$values[which(eigenO$values < 0)])) < min(eigenO$values[which(eigenO$values > 0)])) {
      #this is an easy case where the negative eigenvalues are also the smallest singular values.
      message(cat("We recommend using the pseudoinverse, with a SVD threshold of above",
                  max(abs(eigenO$values[which(eigenO$values < 0)])),
                  "and below",
                  min(eigenO$values[which(eigenO$values > 0)]), sep = " "))
    } else {
      #a more complicated situation where there isn't complete separation between the non-negative and negative eigenvalues
      #(that is, the negative eigenvalues aren't also the smallest singular values)
      #in this case, we can either just ignore the lack of perfect separation, or perform an eigenvalue based removal
      message(cat("Negative eigenvalues are not also the smallest singular values:",
                  "You can either use a SVD threshold of",
                  max(abs(eigenO$values[which(eigenO$values < 0)])),
                  "or use the 'eigen' argument for SVD_thresh", sep = " "))
    }
    return(max(abs(eigenO$values[which(eigenO$values < 0)])))
  } else {
    message("Matrix is positive definite!")
    return(0)
  }
}

#' @export
eigen.pinv <- function(Omega) {
  #an eigenvalue based pseudoinverse
  eigenO <- eigen(Omega)
  svdO <- svd(Omega)

  if (any(eigenO$values < 0)) {
    neg_eigen <- eigenO$values[eigenO$values < 0]
    dinv <- vector()
    for (j in 1:length(svdO$d)) {
      if (any(abs(svdO$d[j] + neg_eigen) < 1e-15)) {
        dinv[j] <- 0 #remove the singular value because it equals a negative eigenvalue
      } else {
        dinv[j] <- 1/svdO$d[j] #keep the singular value, and invert it
      }
    }
  } else { #the matrix is positive definite, no worries.
    dinv <- 1/svdO$d
  }
  return(svdO$v %*% diag(dinv) %*% t(svdO$u))
}

#' @export
matrix_thin <- function(matrix, R2 = 0.99999, verbose = FALSE) {
  n <- dim(matrix)[1]
  rm_index <- vector()
  for (i in 1:(n - 1)) {
    lds <- matrix[i, (i+1):n]^2
    redundant <- which(lds > R2) + i
    rm_index <- c(rm_index, redundant)
    rm_index <- unique(rm_index)
  }

  if (length(rm_index) >= 1) {
    if (verbose) {message(paste("Removed entries", rm_index))}
  } else {
    if (verbose){message("No variants removed by thinning!")}
  }

  matrix_thinned <- matrix[-rm_index, -rm_index]
  if (length(matrix_thinned) == 1) {
    warning("All but one variant removed - the score statistic needs at least two variants.")
  }

  return(matrix_thinned)
}
