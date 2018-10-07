test_that("lik matches exact calculation", {
  # Double check lik is correct using full equations.
  # Can only be done for small sample sizes this way.
  
  # Create covariance function for vectors
  CorrMat32Vecs <- function(u,v, theta) {
    d <- abs(u-v)
    prod((1+d/theta/sqrt(3)) * exp(-d/theta/sqrt(3)))
  }
  CorrMat32Full <- function(X, theta) {
    n <- nrow(X)
    outer(1:n, 1:n, Vectorize(function(i,j) {CorrMat32Vecs(X[i,],X[j,],theta)}))
  }
  
  # Use borehole
  borehole <- function(x) {
    rw <- x[, 1] * (0.15 - 0.05) + 0.05
    r <-  x[, 2] * (50000 - 100) + 100
    Tu <- x[, 3] * (115600 - 63070) + 63070
    Hu <- x[, 4] * (1110 - 990) + 990
    Tl <- x[, 5] * (116 - 63.1) + 63.1
    Hl <- x[, 6] * (820 - 700) + 700
    L <-  x[, 7] * (1680 - 1120) + 1120
    Kw <- x[, 8] * (12045 - 9855) + 9855
    
    m1 <- 2 * pi * Tu * (Hu - Hl)
    m2 <- log(r / rw)
    m3 <- 1 + 2 * L * Tu / (m2 * rw ^ 2 * Kw) + Tu / Tl
    return(m1 / m2 / m3)
  }
  d = 8
  testf<-function (x) {  return(borehole(x))} 
  
  # Create SG to test on
  SG = SGcreate(rep(0, d), rep(1, d),21) #create the design.  it has so many entries because i am sloppy
  Y = testf(SG$design) #the design is $design, simple enough, right?
  logtheta <- rep(0,8)
  lk <- lik(logtheta = logtheta, SG = SG, y = Y)
  expect_equal(lk,
               log(sum(Y * solve(CorrMat32Full(SG$design, theta=exp(logtheta)), Y))/length(Y)) + 
                 log(det(CorrMat32Full(SG$design, theta=exp(logtheta)))) / length(Y) +
                 sum(logtheta^2)/length(Y)
  )
  
  SG = SGcreate(rep(0, d), rep(1, d),31) #create the design.  it has so many entries because i am sloppy
  Y = testf(SG$design) #the design is $design, simple enough, right?
  logtheta <- -(1:8)/10
  lk <- lik(logtheta = logtheta, SG = SG, y = Y)
  expect_equal(lk,
               log(sum(Y * solve(CorrMat32Full(SG$design, theta=exp(logtheta)), Y))/length(Y)) + 
                 log(det(CorrMat32Full(SG$design, theta=exp(logtheta)))) / length(Y) +
                 sum(logtheta^2)/length(Y)
  )
  
  SG = SGcreate(rep(0, d), rep(1, d),41) #create the design.  it has so many entries because i am sloppy
  Y = testf(SG$design) #the design is $design, simple enough, right?
  logtheta <- rep(.1,8)
  lk <- lik(logtheta = logtheta, SG = SG, y = Y)
  expect_equal(lk,
               log(sum(Y * solve(CorrMat32Full(SG$design, theta=exp(logtheta)), Y))/length(Y)) + 
                 log(det(CorrMat32Full(SG$design, theta=exp(logtheta)))) / length(Y) +
                 sum(logtheta^2)/length(Y)
  )
  
  
})