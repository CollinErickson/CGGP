# -----------------------------------------------
# Test different output parameter dimensions
# -----------------------------------------------

# Use 3 dim output. 3rd func is LinComb of first two, so PCA should have 2 dim
f1 <- function(x){x[1]+x[2]^2 + cos(x[3]^2*2*pi*4) - 3.3}
f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
f3 <- function(x) {f1(x) + .3*f2(x)}
f4 <- function(x) {-1.1*f1(x) + .8*f2(x)}
f5 <- function(x) {.2*f1(x)}
f <- function(x) {
  if (is.matrix(x)) {return(t(apply(x, 1, f)))}
  c(f1(x), f2(x), f3(x), f4(x), f5(x))#, rep(f5(x),20))
}
d <- 3
outd <- 5
outd_pca <- 2

nsup <- 20
xsup <- matrix(runif(nsup*d), nsup, d)
ysup <- f(xsup)
eps.sup <- 1e-2 # Use difference accuracy for supp data preds


ntest <- 20
xtest <- matrix(runif(ntest*d), ntest, d)
ytest <- f(xtest)
eps.test <- 1e-2 # Use difference accuracy for testp data preds


test_that("1. MV output, PCA, 1opd", {
  
  # First check MV with PCA
  # Keep sample size small to pred on Xsup will be closer to exact
  SG <- SGGPcreate(d=d, batchsize=30)
  expect_is(SG, "SGGP")
  y <- f(SG$design)
  expect_error(SG <- SGGPfit(SG, Y=y), NA) # No error
  expect_length(SG$thetaMAP, d*SG$numpara)
  expect_true(!is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$Y) == outd)
  expect_true(ncol(SG$y) == outd_pca)
  yMVpred <- SGGPpred(SG$design, SGGP=SG)$mean
  expect_equal(yMVpred, y, 1e-4)
  
  # Check that append works without error, don't save it
  for (sel.method in c("UCB", "TS", "Greedy")) {
    expect_error(SGGPappend(SG, 30, sel.method), NA)
  }
  
  # Add supplemental data
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- SGGPpred(SG, xsup)$me
  expect_equal(ysuppred, ysup, eps.sup)
  
  
  # Add supplemental data, but fewer rows than 2*d so it runs other code in fit
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup[1:2,], Ys=ysup[1:2,]), NA) # No error
  ysuppred <- SGGPpred(SG, xsup[1:2,])$me
  expect_equal(ysuppred, ysup[1:2,], tol=eps.sup)
  
  # print(c((mean(abs(ysuppred[,1] - ysup1[1:2]))),
  # (mean(abs(ysuppred[,2] - ysup2[1:2]))),
  # (mean(abs(ysuppred[,3] - ysup3[1:2])))))
})


test_that("2. MV output, NO PCA, 1opd", {
  
  
  # First check MV with PCA
  SG <- SGGPcreate(d=d, batchsize=30)
  expect_is(SG, "SGGP")
  y <- f(SG$design)
  expect_error(SG <- SGGPfit(SG, Y=y, use_PCA = F), NA) # No error
  expect_length(SG$thetaMAP, d*SG$numpara)
  expect_true(!is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$Y) == outd)
  expect_true(ncol(SG$y) == outd)
  yMVpred <- SGGPpred(SG, SG$design)$mean
  expect_equal(yMVpred, y, 1e-4)
  
  # Check that append works without error, don't save it
  for (sel.method in c("UCB", "TS", "Greedy")) {
    expect_error(SGGPappend(SG, 30, sel.method, NA))
  }
  
  # Add supplemental data
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- SGGPpred(SG, xsup)$me
  expect_equal(ysuppred, ysup, eps.sup)
})


test_that("3. MV output, PCA, separate opd", {
  # First check MV with PCA
  SG <- SGGPcreate(d=d, batchsize=30)
  expect_is(SG, "SGGP")
  y <- f(SG$design)
  expect_error(SG <- SGGPfit(SG, Y=y, use_PCA = T, separateoutputparameterdimensions = T), NA) # No error
  expect_length(SG$thetaMAP, d*SG$numpara*outd_pca)
  expect_true(is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$thetaMAP) == outd_pca)
  expect_true(ncol(SG$Y) == outd)
  expect_true(ncol(SG$y) == outd_pca)
  yMVpred <- SGGPpred(SG, SG$design)$mean
  expect_equal(yMVpred, y, 1e-4)
  
  # Check that append works without error, don't save it
  for (sel.method in c("UCB", "TS", "Greedy")) {
    expect_error(SGGPappend(SG, 30, sel.method), NA)
  }
  
  # Append new data, use RIMSEperpoint
  SG <- SGGPappend(SG, 100, RIMSEperpoint = TRUE, selectionmethod = "UCB")
  y <- f(SG$design)
  expect_error(SG <- SGGPfit(SG, Y=y), NA) # No error
  
  # Add supplemental data
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- SGGPpred(SG, xsup)
  expect_equal(ysuppred$mean, ysup, eps.sup)
  expect_true(all(!is.na(ysuppred$var)))
  expect_true(all((ysuppred$var>0)))
})


test_that("4. MV output, NO PCA, separate opd", {
  
  
  # First check MV with PCA
  SG <- SGGPcreate(d=3, batchsize=30)
  expect_is(SG, "SGGP")
  y <- f(SG$design)
  
  # Fit a model to only each individual output, will compare thetaMAP later
  seed <- sample(1:10000, 1)
  set.seed(seed)
  SG1 <- SGGPfit(SG, Y=y[,1])
  SG2 <- SGGPfit(SG, Y=y[,2])
  SG3 <- SGGPfit(SG, Y=y[,3])
  SG4 <- SGGPfit(SG, Y=y[,4])
  SG5 <- SGGPfit(SG, Y=y[,5])
  
  # Now fit all
  set.seed(seed)
  expect_error(SG <- SGGPfit(SG, Y=y, use_PCA = F, separateoutputparameterdimensions = T), NA) # No error
  expect_length(SG$thetaMAP, d*SG$numpara*outd)
  expect_true(is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$thetaMAP) == outd)
  expect_true(ncol(SG$Y) == outd)
  expect_true(ncol(SG$y) == outd)
  
  # Check that thetaMAP for each dimension matches when model is only fit to that dimension.
  expect_equal(SG$thetaMAP[,1], SG1$thetaMAP)
  expect_equal(SG$thetaMAP[,2], SG2$thetaMAP)
  expect_equal(SG$thetaMAP[,3], SG3$thetaMAP)
  expect_equal(SG$thetaMAP[,4], SG4$thetaMAP)
  expect_equal(SG$thetaMAP[,5], SG5$thetaMAP)
  # Check that predictions on these match
  expect_equal(c(SGGPpred(SG, xtest)$me[,1]), c(SGGPpred(SG1, xtest)$me))
  expect_equal(c(SGGPpred(SG, xtest)$me[,2]), c(SGGPpred(SG2, xtest)$me))
  expect_equal(c(SGGPpred(SG, xtest)$me[,3]), c(SGGPpred(SG3, xtest)$me))
  expect_equal(c(SGGPpred(SG, xtest)$me[,4]), c(SGGPpred(SG4, xtest)$me))
  expect_equal(c(SGGPpred(SG, xtest)$me[,5]), c(SGGPpred(SG5, xtest)$me))
  
  # Now check predictions
  yMVpred <- SGGPpred(SG$design, SGGP=SG)$mean
  expect_equal(yMVpred, y, 1e-4)
  
  # Check that append works without error, don't save it
  for (sel.method in c("UCB", "TS", "Greedy")) {
    expect_error(SGGPappend(SG, 30, sel.method), NA)
  }
  
  # Add supplemental data
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- SGGPpred(SG, xsup)$me
  expect_equal(ysuppred, ysup, eps.sup)
})

