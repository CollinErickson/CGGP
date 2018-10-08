test_that("Prediction matches exact on small samples", {
  
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
  
  # Create covariance function for vectors
  CorrMat32Vecs <- function(u,v, theta) {
    d <- abs(u-v)
    prod((1+d/theta/sqrt(3)) * exp(-d/theta/sqrt(3)))
  }
  CorrMat32Full <- function(X, theta) {
    n <- nrow(X)
    outer(1:n, 1:n, Vectorize(function(i,j) {CorrMat32Vecs(X[i,],X[j,],theta)}))
  }
  CorrMat32Full2 <- function(X, U, theta) {
    outer(1:nrow(X), 1:nrow(U), Vectorize(function(i,j) {CorrMat32Vecs(X[i,],U[j,],theta)}))
  }
  
  SG = SGcreate(rep(0, d), rep(1, d),101) #create the design.  it has so many entries because i am sloppy
  Y = testf(SG$design) #the design is $design, simple enough, right?
  logtheta <- rep(0,8)
  logtheta <- logthetaMLE(SG=SG,y=Y)
  
  n <- 50
  xp <- matrix(runif(d*10),n,8)
  SGpred <- SGGPpred(xp=xp, SG=SG, y=Y, logtheta = logtheta)
  
  my <- mean(Y)
  dy <- Y - my
  Sig <- CorrMat32Full(SG$design, theta=exp(logtheta))
  s <- CorrMat32Full2(xp, SG$design, theta = exp(logtheta))
  expred <- my + s %*% solve(Sig, dy)
  
  plot(expred, SGpred$mean); abline(a=0,b=1, col=2)
  expect_equal(SGpred$mean, expred)
  
})  
  
