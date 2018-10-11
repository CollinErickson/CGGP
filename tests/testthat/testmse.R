test_that("MSE_calc for Matern 3/2 works", {
  # These all assume sigma2 is 1
  lth = -1
  xl <- c(0,.5,.9)
  nugget <- 0
  MSE_calc.out <- MSE_calc(xl=xl, logtheta=lth, nugget=nugget)
  
  S = CorrMat(xl, xl, logtheta = lth)
  Ci = solve(S)
  matern32 <- function(x,lth) {(1+abs(x)/exp(lth)/sqrt(3))*exp(-abs(x)/exp(lth)/sqrt(3))}
  integrand <- Vectorize(function(x) {1+nugget-matern32(x-xl,lth=lth) %*% Ci %*% matern32(x-xl,lth=lth)})
  curve(integrand); abline(v=xl)
  integrate.out <- integrate(integrand, lower=0, upper=1)
  expect_equal(MSE_calc.out, integrate.out$value, tol=1e-4)
})

test_that("MSE_calc for Matern 3/2 works with nugget", {
  lth = -1
  xl <- c(0,.5,.9)
  nugget <- 1e-3
  MSE_calc.out <- MSE_calc(xl=xl, logtheta=lth, nugget=nugget)
  
  S = CorrMat(xl, xl, logtheta = lth)
  diag(S) <- diag(S) + nugget
  Ci = solve(S)
  matern32 <- function(x,lth) {(1+abs(x)/exp(lth)/sqrt(3))*exp(-abs(x)/exp(lth)/sqrt(3))}
  integrand <- Vectorize(function(x) {1+nugget-matern32(x-xl,lth=lth) %*% Ci %*% matern32(x-xl,lth=lth)})
  curve(integrand); abline(v=xl)
  integrate.out <- integrate(integrand, lower=0, upper=1)
  expect_equal(MSE_calc.out, integrate.out$value, tol=1e-4)
})


test_that("MSE_calc for Gauss works", {
  th = .1
  xl <- c(0,.5,.9)
  MSE_calc.out <- MSEgauss(xl=xl, theta=th)
  
  S = gausscorr(xl, xl, th)
  Ci = solve(S)
  gauss <- function(x,th) {exp(-(x)^2/th)}
  integrand <- Vectorize(function(x) {1-gauss(x-xl,th=th) %*% Ci %*% gauss(x-xl,th=th)})
  curve(integrand); abline(v=xl)
  integrate.out <- integrate(integrand, lower=0, upper=1)
  expect_equal(MSE_calc.out, integrate.out$value, tol=1e-4)
})
