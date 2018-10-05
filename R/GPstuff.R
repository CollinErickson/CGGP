
#' Calculate correlation matrix for two sets of points in one dimension
#' 
#' Note that this isn't the correlation between two vectors.
#' It is for two
#'
#' @param x1 Vector of coordinates from same dimension
#' @param x2 Vector of coordinates from same dimension
#' @param theta Is this log theta???
#'
#' @return Matrix
#' @export
#'
#' @examples
#' CorrMat(c(0,.2,.4),c(.1,.3,.5), 10)
CorrMat <- function(x1, x2, .., logtheta, theta) {#browser()
  if (missing(theta)) { theta <- exp(logtheta)}
  thetasqrt3 <- theta*sqrt(3)
  d1 = length(x1)
  d2 = length(x2)
  # Would be easier to use outer
  diffmat = abs(t(matrix(rep(x1, each = d2), nrow = d2)) - (matrix(rep(x2, each =
                                                                         d1), nrow = d1)))
  # This is close to Matern 3/2 but missing sqrt(3) in both places
  C = (1 + diffmat / thetasqrt3) * exp(-diffmat / thetasqrt3)
  return(C)
}



#' Calculate derivative of CorrMat with respect to theta
#'
#' @param x1 First vector of 1D points
#' @param x2 Second vector of 1D points
#' @param theta Actually log theta, terrible notation
#'
#' @return Matrix
#' @export
#'
#' @examples
#' dCorrMat(c(0,.2,.4),c(.1,.3,.5), 10)
dCorrMat <- function(x1, x2, ..., logtheta, theta) {
  if (missing(theta)) {theta <- exp(logtheta)}
  thetasqrt3 <- theta*sqrt(3)
  d1 = length(x1)
  d2 = length(x2)
  diffmat = abs(t(matrix(rep(x1, each = d2), nrow = d2)) - (matrix(rep(x2, each =
                                                                         d1), nrow = d1)))
  C = (diffmat / thetasqrt3) ^ 2 * exp(-diffmat / thetasqrt3)
  return(C)
}
