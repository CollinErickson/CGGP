test_that("Correlation CorrMatCauchy works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- c(.05,.9,-.3)
  th <- runif(3,-1,1)
  
  # First check return_numpara is right
  expect_equal(SGGP_internal_CorrMatCauchy(return_numpara=TRUE), 3)
  
  # Now check correlation
  cauchy1 <- SGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th)
  expect_is(cauchy1, "matrix")
  expect_equal(dim(cauchy1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  cauchyfunc <- function(a,b,theta) {
    expLS <- exp(3*theta[1])
    expHE <- exp(3*theta[2])
    alpha = 2*exp(3*theta[3]+4)/(1+exp(3*theta[3]+4))
    diffmat <- outer(a, b, Vectorize(function(aa,bb) abs(aa-bb)))
    h = diffmat/expLS
    halpha = h^alpha
    pow = -expHE/alpha
    (1+halpha)^pow
  }
  cauchy2 <- cauchyfunc(x1, x2, theta=th)
  expect_equal(cauchy1, cauchy2)
  
  # Now check that dC is actually grad of C
  cauchy_C_dC <- SGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:3) {
    thd <- c(0,0,0)
    thd[i] <- eps
    numdC <- (SGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th+thd) -
              SGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # numdC <- (-SGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th+2*thd) +
    #           8*  SGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th+thd) + 
    #           -8*  SGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th-thd) +
    #             SGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th-2*thd)) / (12*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, cauchy_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})
