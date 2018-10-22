# to change correlation function: write lik, dlik, and IMSE

#' Gaussian error function
#'
#' @param x Points to evaluate at
#'
#' @importFrom stats integrate optim pnorm
#' @importFrom Rcpp evalCpp
#' @return ERF of x
#' @export
#'
#' @examples
#' erf(1)
erf <- function(x) {
  2*pnorm(sqrt(2)*x) - 1
}

#' Gaussian correlation matrix
#'
#' @param x1 Vector of 1D points
#' @param x2 Vector of 1D points
#' @param theta Correlation parameter in 1D, one number
#' @param logtheta Log of correlation parameters
#' @param nugget Nugget to add to each correlation.
#'
#' @return Correlation matrix
#' @export
#'
#' @examples
#' gausscorr(c(0,.5,1), c(.1,.4,.5,.9), theta=.1, nugget=0)
gausscorr <- function(x1, x2, theta, logtheta, nugget) {
  if (missing(theta)) { theta <- exp(logtheta)}
  outer(x1, x2, function(a,b) {exp(-(a-b)^2/theta)})
}

#' Derivative of Gaussian correlation matrix
#'
#' @inheritParams gausscorr
#'
#' @return Derivative of Gaussian correlation matrix in 1D w.r.t. theta
#' @export
#'
#' @examples
#' dgausscorr(c(0,.5,1), c(.1,.4,.5,.9), theta=.1, nugget=0)
dgausscorr <- function(x1, x2, theta, logtheta, nugget) {
  if (missing(theta)) { theta <- exp(logtheta)}
  outer(x1, x2, function(a,b) {exp(-(a-b)^2/theta) * (a-b)^2/theta^2})
}

#' Gaussian correlation matrix
#'
#' @inheritParams gausscorr
#'
#' @return Correlation matrix
#' @export
#'
#' @examples
#' gausscorr(c(0,.5,1), c(.1,.4,.5,.9), theta=.1, nugget=0)
diag_gausscorr <- function(x1, theta, logtheta, nugget) {
  rep(1+nugget, length(x1))
}

#' Calculate MSE for Gaussian correlation
#'
#' @param xl Vector of evaluated points in 1D
#' @param theta Correlation parameter, one number
#'
#' @return MSE
#' @export
#'
#' @examples
#' th = .1
#' xl <- c(0,.5,.9)
#' MSE_calc.out <- MSEgauss(xl=xl, theta=th)
MSEgauss <- function(xl, theta) {
  S = gausscorr(xl, xl, theta)
  #t = exp(theta)
  n = length(xl)
  Ci = solve(S)
  
  
  A = matrix(rep(xl, each = n), nrow = n)
  a = pmin(A, t(A))
  b = pmax(A, t(A))
  
  # Calculated using symbolic integration. Checked in tests/testmse.R compared to numerical integration.
  out0 <- exp(-(a-b)^2/(2*theta)) * (erf((a+b)/sqrt(2*theta)) - erf((a+b-2)/sqrt(2*theta))) # Make matrix
  out1 <- .5*sqrt(pi/2)*sqrt(theta)*out0 # Just scaling
  #MSE = 1 - sum(diag(out1 %*% Ci))
  MSE = 1 - sum((out1 * Ci))
  MSE
}