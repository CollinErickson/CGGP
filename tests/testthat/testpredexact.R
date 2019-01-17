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
  CorrMatCauchySQTPts <- function(x1, x2, theta) {
    
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
  CorrMatCauchySQTVecs <- function(x1,x2, theta) {
    prod(sapply(1:length(x1), function(i) {
      CorrMatCauchySQTPts(x1[i], x2[i], theta[(1+3*i-3):(3*i)])
    }))
  }
  CorrMatCauchySQTFull <- function(X, theta) {
    n <- nrow(X)
    outer(1:n, 1:n, Vectorize(function(i,j) {CorrMatCauchySQTVecs(X[i,],X[j,],theta)}))
  }
  CorrMatCauchySQTFull2 <- function(X, U, theta) {
    outer(1:nrow(X), 1:nrow(U), Vectorize(function(i,j) {CorrMatCauchySQTVecs(X[i,],U[j,],theta)}))
  }
  
  SG = SGGPcreate(d=d,31) # Need size to be small to avoid computationally singular in solve
  Y = testf(SG$design) #the design is $design, simple enough, right?
  SG <- SGGPfit(SG=SG, Y=Y)
  
  n <- 50
  xp <- matrix(runif(d*n),n,d)
  SGpred <- SGGPpred(xp=xp, SG=SG)
  
  my <- mean(Y)
  dy <- Y - my
  Sig <- CorrMatCauchySQTFull(SG$design, theta=SG$thetaMAP)
  s <- CorrMatCauchySQTFull2(xp, SG$design, theta = SG$thetaMAP)
  expred <- my + s %*% solve(Sig, dy)
  
  # Check mean predictions
  # plot(expred, SGpred$mean); abline(a=0,b=1, col=2)
  expect_equal(SGpred$mean, expred)
  
  # Test var predictions
  # Calculating s2 like this doesn't work since we use MAP
  # s2 <- c(t(Y) %*%solve(Sig, Y) / length(Y))
  # Just use the MAP value
  s2 <- SG$sigma2MAP[1,1]
  exvar <- s2 * (1 - colSums(t(s) * solve(Sig, t(s))))
  print(1/SGpred$var* exvar)
  plot(exvar, SGpred$var); abline(a=0,b=1, col=2)
  expect_equal(SGpred$var, exvar)
})

test_that("predMV works", {
  SG <- SGGPcreate(d=3, batchsize=100)
  f1 <- function(x){x[1]+x[2]^2}
  f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
  y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
  y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
  y <- cbind(y1, y2)
  SG <- SGGPfit(SG, Y=y)
  yMVpred <- SGGPpred(SG$design, SG=SG)$mean
  expect_equal(yMVpred[,1], y1, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  
  # Doesn't work since there's no way to update Y without updating parameters too.
  # xpred <- matrix(runif(100*3),100,3)
  # SG1 <- SGGPfit(SG, Y=y1)
  # SG2 <- SGGPfit(SG, Y=y2)
  # y1pred <- SGGPpred(xpred, SG=SG1)$mean
  # y2pred <- SGGPpred(xpred, SG=SG2)$mean
  # yMVpred <- SGGPpred(xpred, SG=SG)$mean
  # expect_equal(yMVpred[,1], c(y1pred), tol=1e-2)
  # expect_equal(yMVpred[,2], c(y2pred), tol=1e-2)
})

test_that("Supplemented works", {
  d <- 3
  SG <- SGGPcreate(d=d, batchsize=100)
  f1 <- function(x){x[1]+sin(2*pi*x[1]) + x[2]^2}
  y1 <- apply(SG$design, 1, f1)
  SG <- SGGPfit(SG, Y=y1)
  
  # Add supplemental data
  nsup <- 20
  xsup <- matrix(runif(d*nsup), nsup, d)
  ysup <- apply(xsup, 1, f1)
  # SG$supplemented <- TRUE
  # SG$Xs <- xsup
  # SG$Ys <- ysup
  
  # # Get error when not fit
  # expect_error(SGGPpred(xp=xsup, SG=SG))
  
  # Should work after fitting
  SG <- SGGPfit(SG, Y=y1, Xs=xsup, Ys=ysup)
  
  # Predictions should match values at supplemented points
  expect_equal(c(SGGPpred(xp=xsup, SG=SG)$me), ysup, tol=1e-4)
  
  # # Predict at points
  # n <- 50
  # xp <- matrix(runif(d*10),n,d)
  # SGpred <- SGGPpred(xp=xp, SG=SG)
})

test_that("supplemental with MV output works", {
  
  SG <- SGGPcreate(d=3, batchsize=100)
  f1 <- function(x){x[1]+x[2]^2}
  f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
  y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
  y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
  y <- cbind(y1, y2)
  
  xsup <- matrix(runif(3*30), ncol=3)
  ysup1 <- apply(xsup, 1, f1)
  ysup2 <- apply(xsup, 1, f2)
  ysup <- cbind(ysup1, ysup2)
  
  SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup)
  yMVpred <- SGGPpred(SG$design, SG=SG)$mean
  expect_equal(yMVpred[,1], y1, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  
})
