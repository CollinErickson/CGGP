# -----------------------------------------------
# Test different output parameter dimensions
# -----------------------------------------------

# Use 3 dim output. 3rd func is LinComb of first two, so PCA should have 2 dim
f1 <- function(x){x[1]+x[2]^2 + cos(x[3]^2*2*pi*4) - 3.3}
f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
f3 <- function(x) {f1(x) + .3*f2(x)}

nsup <- 20
xsup <- matrix(runif(nsup*3), nsup, 3)
ysup1 <- apply(xsup, 1, f1)
ysup2 <- apply(xsup, 1, f2)
ysup3 <- apply(xsup, 1, f3)
ysup <- cbind(ysup1, ysup2, ysup3)
eps.sup <- 1e-2 # Use difference accuracy for supp data preds


ntest <- 20
xtest <- matrix(runif(ntest*3), ntest, 3)
ytest1 <- apply(xtest, 1, f1)
ytest2 <- apply(xtest, 1, f2)
ytest3 <- apply(xtest, 1, f3)
ytest <- cbind(ytest1, ytest2, ytest3)
eps.test <- 1e-2 # Use difference accuracy for testp data preds


test_that("1. MV output, PCA, 1opd", {
  
  # First check MV with PCA
  # Keep sample size small to pred on Xsup will be closer to exact
  SG <- SGGPcreate(d=3, batchsize=30)
  expect_is(SG, "SGGP")
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y3 <- apply(SG$design, 1, f3)
  y <- cbind(y1, y2, y3)
  expect_error(SG <- SGGPfit(SG, Y=y), NA) # No error
  expect_length(SG$thetaMAP, 3*3)
  expect_true(!is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$Y) == 3)
  expect_true(ncol(SG$y) == 2)
  yMVpred <- SGGPpred(SG$design, SG=SG)$mean
  expect_equal(yMVpred[,1], y1, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  
  # Add supplemental data
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- SGGPpred(xsup, SG)$me
  expect_equal(ysuppred[,1], ysup1, eps.sup)
  expect_equal(ysuppred[,2], ysup2, eps.sup)
  expect_equal(ysuppred[,3], ysup3, eps.sup)
  
  
  # Add supplemental data, but fewer rows than 2*d so it runs other code in fit
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup[1:2,], Ys=ysup[1:2,]), NA) # No error
  ysuppred <- SGGPpred(xsup[1:2,], SG)$me
  expect_equal(ysuppred[,1], ysup1[1:2], eps.sup)
  expect_equal(ysuppred[,2], ysup2[1:2], eps.sup)
  expect_equal(ysuppred[,3], ysup3[1:2], eps.sup)
  
  # print(c((mean(abs(ysuppred[,1] - ysup1[1:2]))),
  # (mean(abs(ysuppred[,2] - ysup2[1:2]))),
  # (mean(abs(ysuppred[,3] - ysup3[1:2])))))
})


test_that("2. MV output, NO PCA, 1opd", {
  
  
  # First check MV with PCA
  SG <- SGGPcreate(d=3, batchsize=30)
  expect_is(SG, "SGGP")
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y3 <- apply(SG$design, 1, f3)
  y <- cbind(y1, y2, y3)
  expect_error(SG <- SGGPfit(SG, Y=y, use_PCA = F), NA) # No error
  expect_length(SG$thetaMAP, 3*3)
  expect_true(!is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$Y) == 3)
  expect_true(ncol(SG$y) == 3)
  yMVpred <- SGGPpred(SG$design, SG=SG)$mean
  expect_equal(yMVpred[,1], y1, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  
  # Add supplemental data
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- SGGPpred(xsup, SG)$me
  expect_equal(ysuppred[,1], ysup1, eps.sup)
  expect_equal(ysuppred[,2], ysup2, eps.sup)
  expect_equal(ysuppred[,3], ysup3, eps.sup)
})


test_that("3. MV output, PCA, separate opd", {
  
  # First check MV with PCA
  SG <- SGGPcreate(d=3, batchsize=100)
  expect_is(SG, "SGGP")
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y3 <- apply(SG$design, 1, f3)
  y <- cbind(y1, y2, y3)
  expect_error(SG <- SGGPfit(SG, Y=y, use_PCA = T, separateoutputparameterdimensions = T), NA) # No error
  expect_length(SG$thetaMAP, 3*3*2)
  expect_true(is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$thetaMAP) == 2)
  expect_true(ncol(SG$Y) == 3)
  expect_true(ncol(SG$y) == 2)
  yMVpred <- SGGPpred(SG$design, SG=SG)$mean
  expect_equal(yMVpred[,1], y1, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  
  # Append new data, use RIMSEperpoint
  SG <- SGGPappend(SG, 100, RIMSEperpoint = TRUE, selectionmethod = "UCB")
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y3 <- apply(SG$design, 1, f3)
  y <- cbind(y1, y2, y3)
  expect_error(SG <- SGGPfit(SG, Y=y), NA) # No error
  
  # Add supplemental data
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- SGGPpred(xsup, SG)$me
  expect_equal(ysuppred[,1], ysup1, eps.sup)
  expect_equal(ysuppred[,2], ysup2, eps.sup)
  expect_equal(ysuppred[,3], ysup3, eps.sup)
})


test_that("4. MV output, NO PCA, separate opd", {
  
  
  # First check MV with PCA
  SG <- SGGPcreate(d=3, batchsize=100)
  expect_is(SG, "SGGP")
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y3 <- apply(SG$design, 1, f3)
  y <- cbind(y1, y2, y3)
  
  # Fit a model to only each individual output, will compare thetaMAP later
  seed <- sample(1:10000, 1)
  set.seed(seed)
  SG1 <- SGGPfit(SG, Y=y1)
  SG2 <- SGGPfit(SG, Y=y2)
  SG3 <- SGGPfit(SG, Y=y3)
  
  # Now fit all
  set.seed(seed)
  expect_error(SG <- SGGPfit(SG, Y=y, use_PCA = F, separateoutputparameterdimensions = T), NA) # No error
  expect_length(SG$thetaMAP, 3*3*3)
  expect_true(is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$thetaMAP) == 3)
  expect_true(ncol(SG$Y) == 3)
  expect_true(ncol(SG$y) == 3)
  
  # Check that thetaMAP for each dimension matches when model is only fit to that dimension.
  expect_equal(SG$thetaMAP[,1], SG1$thetaMAP)
  expect_equal(SG$thetaMAP[,2], SG2$thetaMAP)
  expect_equal(SG$thetaMAP[,3], SG3$thetaMAP)
  # Check that predictions on these match
  expect_equal(c(SGGPpred(xtest, SG)$me[,1]), c(SGGPpred(xtest, SG1)$me))
  expect_equal(c(SGGPpred(xtest, SG)$me[,2]), c(SGGPpred(xtest, SG2)$me))
  expect_equal(c(SGGPpred(xtest, SG)$me[,3]), c(SGGPpred(xtest, SG3)$me))
  
  # Now check predictions
  yMVpred <- SGGPpred(SG$design, SG=SG)$mean
  expect_equal(yMVpred[,1], y1, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  
  # Add supplemental data
  expect_error(SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- SGGPpred(xsup, SG)$me
  expect_equal(ysuppred[,1], ysup1, eps.sup)
  expect_equal(ysuppred[,2], ysup2, eps.sup)
  expect_equal(ysuppred[,3], ysup3, eps.sup)
})
