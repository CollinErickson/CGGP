# This was replaced when changing to integral form.

#' Calculate MSE over single dimension for Matern 3/2
#' 
#' Equation found using
#'
#' @param xl Vector of points in 1D
#' @param logtheta Log of correlation parameters.
#' @param theta Correlation parameters
#' @param nugget Nugget to add to diagonal of correlation matrix.
#' @param ... Don't use, just forces theta to be named
#'
#' @return MSE value
#' @export
#'
#' @examples
#' MSE_calc(xl=c(0,.5,.9), theta=1, nugget=.001)
MSE_calc_Mat32Exact <- function(xl, ..., logtheta, theta, nugget) {
  if (missing(theta)) {theta <- exp(logtheta)}
  S = CorrMat(xl, xl, theta=theta)
  diag(S) <- diag(S) + nugget
  #t = exp(theta)
  t = theta * sqrt(3)
  n = length(xl)
  Ci = solve(S)
  
  
  A = matrix(rep(xl, each = n), nrow = n)
  a = pmin(A, t(A))
  b = pmax(A, t(A))
  
  t2 = 1.0 / t
  t3 = a + b - 2.0
  t4 = t2 * t3
  t5 = exp(t4)
  t6 = a - b
  t7 = t2 * t6
  t8 = exp(t7)
  t9 = t ^ 2
  t10 = a * t2 * 2.0
  t11 = exp(t10)
  t12 = a * b * 2.0
  
  out1 = t5 * (-3.0 / 2.0) + a * t5 * (3.0 / 4.0) - a * t8 * (3.0 / 4.0) +
    b * t5 * (3.0 / 4.0) + b * t8 * (3.0 / 4.0) - t * t5 * (5.0 / 4.0) - t2 *
    t5 * (1.0 / 2.0) + t * t8 * (5.0 / 4.0) + a * t2 * t5 * (1.0 / 2.0) + b *
    t2 * t5 * (1.0 / 2.0) - t2 * exp(-t2 * (a + b)) * (t9 * 5.0 + t12 + a *
                                                         t * 3.0 + b * t * 3.0 - t9 * t11 * 5.0 + a * t * t11 * 3.0 - b * t * t11 *
                                                         3.0) * (1.0 / 4.0) - a * b * t2 * t5 * (1.0 / 2.0) - 1.0 / t ^ 2 * t6 *
    t8 * (t9 * 6.0 - t12 - a * t * 6.0 + b * t * 6.0 + a ^ 2 + b ^ 2) * (1.0 /
                                                                           6.0)
  
  MSE = diag_corrMat(.5, theta=theta, nugget=nugget) - sum(diag(out1 %*% Ci))
  
  MSE
}


if (F) {
  xl <- seq(0,1,l=33)
  # Make sure values are same
  MSE_calc(xl=xl, theta=.1, nugget=0)
  MSE_calc_Mat32Exact(xl=xl, theta=.1, nugget=0)
  # Check timing
  microbenchmark::microbenchmark(
    MSE_calc(xl=xl, theta=.1, nugget=0),
    MSE_calc_Mat32Exact(xl=xl, theta=.1, nugget=0)
    
  )
  # They take similar time, but for large n the approximate one is faster
}