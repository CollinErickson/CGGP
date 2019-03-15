CorrMat <- function(x1, x2, theta) {
  d1 = length(x1)
  d2 = length(x2)
  diffmat = abs(t(matrix(rep(x1, each = d2), nrow = d2)) - (matrix(rep(x2, each =
                                                                         d1), nrow = d1)))
  C = 10^(-10)*(diffmat<10^(-6))+(1 + diffmat / exp(theta)) * exp(-diffmat / exp(theta))
  return(C)
}

diag_corrMat <- function(x1, theta) {
  d1 = length(x1)
  diffmat = abs(x1-x1)
  C = 10^(-10)*(diffmat<10^(-6))+(1 + diffmat / exp(theta)) * exp(-diffmat / exp(theta))
  return(C)
}

dCorrMat <- function(x1, x2, theta) {
  d1 = length(x1)
  d2 = length(x2)
  diffmat = abs(t(matrix(rep(x1, each = d2), nrow = d2)) - (matrix(rep(x2, each =
                                                                         d1), nrow = d1)))
  C = (diffmat / exp(theta)) ^ 2 * exp(-diffmat / exp(theta))
  return(C)
}
