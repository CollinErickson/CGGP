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

test_that("Correlation CorrMatCauchySQT works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- c(.05,.9,-.3)
  th <- runif(3,-1,1)
  
  # First check return_numpara is right
  expect_equal(SGGP_internal_CorrMatCauchySQT(return_numpara=TRUE), 3)
  
  # Get error when you give in theta of wrong size
  expect_error(SGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  cauchy1 <- SGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta=th)
  expect_is(cauchy1, "matrix")
  expect_equal(dim(cauchy1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  cauchyfunc <- function(x1,x2,theta) {
    
    expTILT = exp((theta[3]))
    expLS = exp(3*(theta[1]))
    x1t = (x1+10^(-2))^expTILT
    x2t = (x2+10^(-2))^expTILT
    x1ts = x1t/expLS
    x2ts = x2t/expLS
    
    diffmat =abs(outer(x1ts,x2ts,'-')); 
    expHE = exp(3*(theta[2]))
    h = diffmat
    alpha = 2*exp(5)/(1+exp(5))
    halpha = h^alpha
    pow = -expHE/alpha
    
    (1+halpha)^pow
  }
  cauchy2 <- cauchyfunc(x1, x2, theta=th)
  expect_equal(cauchy1, cauchy2)
  
  # Now check that dC is actually grad of C
  cauchy_C_dC <- SGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:3) {
    thd <- c(0,0,0)
    thd[i] <- eps
    numdC <- (SGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta=th+thd) -
                SGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, cauchy_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatCauchySQ works", {
  x1 <- runif(5)
  x2 <- runif(4)
  # th <- c(.05,.9)
  th <- runif(2,-1,1)
  
  # First check return_numpara is right
  expect_equal(SGGP_internal_CorrMatCauchySQ(return_numpara=TRUE), 2)
  
  # Get error when you give in theta of wrong size
  expect_error(SGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta = c(.1,.1,.4)))
  
  # Now check correlation
  cauchy1 <- SGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta=th)
  expect_is(cauchy1, "matrix")
  expect_equal(dim(cauchy1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  cauchyfunc <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-')); 
    
    expLS = exp(3*theta[1])
    expHE = exp(3*theta[2])
    h = diffmat/expLS
    alpha = 2*exp(0+6)/(1+exp(0+6))
    halpha = h^alpha
    pow = -expHE/alpha
    
    (1-10^(-10))*(1+halpha)^pow+10^(-10)*(diffmat<10^(-4))
  }
  cauchy2 <- cauchyfunc(x1, x2, theta=th)
  expect_equal(cauchy1, cauchy2)
  
  # Now check that dC is actually grad of C
  cauchy_C_dC <- SGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:2) {
    thd <- c(0,0)
    thd[i] <- eps
    numdC <- (SGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta=th+thd) -
                SGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, cauchy_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatGaussian works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(SGGP_internal_CorrMatGaussian(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(SGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- SGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  gaussianfunc <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-')); 
    diffmat2 <- diffmat^2
    expLS = exp(3*theta[1])
    h = diffmat2/expLS
    C = (1-10^(-10))*exp(-h) + 10^(-10)*(diffmat<10^(-4))
    C
  }
  corr2 <- gaussianfunc(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- SGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (SGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta=th+thd) -
                SGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatMatern32 works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(SGGP_internal_CorrMatMatern32(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(SGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- SGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  matern32func <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    C = (1-10^(-10))*(1+sqrt(3)*h)*exp(-sqrt(3)*h) + 10^(-10)*(diffmat<10^(-4))
    C
  }
  corr2 <- matern32func(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- SGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (SGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta=th+thd) -
                SGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})



test_that("Correlation CorrMatMatern52 works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(SGGP_internal_CorrMatMatern52(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(SGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- SGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  matern52func <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    C = (1-10^(-10))*(1+sqrt(5)*h+5/3*h^2)*exp(-sqrt(5)*h) + 10^(-10)*(diffmat<10^(-4))
    C
  }
  corr2 <- matern52func(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- SGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (SGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta=th+thd) -
                SGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatPowerExp works", {
  nparam <- 2
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(nparam,-1,1)
  
  # First check return_numpara is right
  expect_equal(SGGP_internal_CorrMatPowerExp(return_numpara=TRUE), nparam)
  
  # Get error when you give in theta of wrong size
  expect_error(SGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta = rep(0, nparam-1)))
  expect_error(SGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta = rep(0, nparam+1)))
  
  # Now check correlation
  corr1 <- SGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  PowerExpfunc <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    alpha <- 1 + (theta[2]+1)/2
    h = diffmat/expLS
    C = (1-10^(-10))*exp(-(h)^alpha) + 10^(-10)*(diffmat<10^(-4))
    C
  }
  corr2 <- PowerExpfunc(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- SGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:nparam) {
    thd <- rep(0, nparam)
    thd[i] <- eps
    numdC <- (SGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta=th+thd) -
                SGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})
