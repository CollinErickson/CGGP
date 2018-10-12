
#' Calculate correlation matrix for two sets of points in one dimension
#' 
#' Note that this isn't the correlation between two vectors.
#' It is for two
#'
#' @param x1 Vector of coordinates from same dimension
#' @param x2 Vector of coordinates from same dimension
#' @param logtheta Log of correlation parameter
#' @param theta Correlation parameter
#' @param ... Don't use, just forces theta to be named
#'
#' @return Matrix
#' @export
#'
#' @examples
#' CorrMat(c(0,.2,.4),c(.1,.3,.5), theta=.1)
CorrMat <- function(x1, x2, ..., logtheta, theta) {#browser()
  if (missing(theta)) { theta <- exp(logtheta)}
  if (any(theta<0)) stop("Theta < 0 in CorrMat")
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



#' Calculate correlation of points with themselves
#' 
#' Use to get diagonal of correlation matrix.
#'
#' @param x1 Vector of 1D points
#' @param theta Correlation parameters
#' @param logtheta Log of correlation parameters
#' @param nugget Nugget to add to each correlation.
#' @param ... Don't use, just forces theta to be named
#'
#' @return Correlations of x1 with itself
#' @export
#'
#' @examples
#' diag_corrMat(c(0,.2,.4,.6,.8,1),theta=.1, nugget=0) # Should be all ones
diag_corrMat <- function(x1, ..., logtheta, theta, nugget) {
  if (missing(theta)) { theta <- exp(logtheta)}
  if (any(theta<0)) stop("Theta < 0 in CorrMat")
  thetasqrt3 <- theta*sqrt(3)
  d1 = length(x1)
  diffmat = abs(x1-x1) # just zeros, rep(0, d1)
  #C = 10^(-10)*(diffmat<10^(-6))+(1 + diffmat / exp(theta)) * exp(-diffmat / exp(theta))
  C = nugget + (1 + diffmat / thetasqrt3) * exp(-diffmat / thetasqrt3) # Just ones
  return(C)
}


#' Calculate derivative of CorrMat with respect to logtheta
#'
#' @param x1 First vector of 1D points
#' @param x2 Second vector of 1D points
#' @param logtheta Log correlation parameter
#' @param theta Correlation parameter
#' @param ... Don't use, just forces theta to be named
#'
#' @return Matrix
#' @export
#'
#' @examples
#' dCorrMat(c(0,.2,.4),c(.1,.3,.5), theta=.1)
dCorrMat <- function(x1, x2, ..., logtheta, theta) {
  if (missing(theta)) {theta <- exp(logtheta)}
  thetasqrt3 <- theta*sqrt(3)
  d1 = length(x1)
  d2 = length(x2)
  diffmat = abs(t(matrix(rep(x1, each = d2), nrow = d2)) - (matrix(rep(x2, each =
                                                                         d1), nrow = d1)))
  dC = (diffmat / thetasqrt3) ^ 2 * exp(-diffmat / thetasqrt3)
  return(dC)
}
