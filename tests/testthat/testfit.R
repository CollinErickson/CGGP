# SGGPfit is used in many other tests, so not as many here

test_that("SGGPfit works with Laplace approx", {
  d <- 3
  SG <- SGGPcreate(d=d, batchsize=30)
  f <- function(x){x[1]+x[2]^2}
  y <- apply(SG$design, 1, f)
  SG <- SGGPfit(SG, Y=y)
  
  # Check that samples from Laplace are near max
  if (F) {
    stripchart(as.data.frame(t(SG$thetaPostSamples)))
    stripchart(as.data.frame(t(SG$thetaMAP)), add=T, col=2, pch=19)
  }
  expect_true(all(abs(rowMeans(SG$thetaPostSamples)- SG$thetaMAP) < apply(SG$thetaPostSamples, 1, sd)))
  
  neglogpost <- SGGP_internal_gneglogpost(rep(.5,length(SG$thetaMAP)), SG, y)
  expect_is(neglogpost, "matrix")
  expect_length(neglogpost, length(SG$thetaMAP))
  
  
  neglogpost2 <- SGGP_internal_gneglogpost(rep(.5,length(SG$thetaMAP)), SG, y, return_lik = T)
  expect_is(neglogpost2, "list")
  expect_length(neglogpost2[[1]], 1)
  expect_length(neglogpost2[[2]], length(SG$thetaMAP))
  
  # Check neglogpost grad matches gneglogpost
  theta <- SG$thetaMAP / 2 # Don't want values near -1 or +1
  epsval <- 1e-4
  thetagrad <- SGGP_internal_gneglogpost(theta, SG, SG$y)
  for (i in 1:length(SG$thetaMAP)) {
    eps <- rep(0, length(SG$thetaMAP))
    eps[i] = eps[i] + epsval
    numgrad <- (SGGP_internal_neglogpost(theta + eps, SG, SG$y) - 
                  SGGP_internal_neglogpost(theta - eps, SG, SG$y)) / (2*epsval)
    expect_equal(thetagrad[i], numgrad, tol=1e-4)
  }
  
  # Works with supplementary data
  nsup <- 30
  xsup <- matrix(runif(nsup*3), nsup, 3)
  ysup <- apply(xsup, 1, f)
  SG <- SGGPappend(SG, 20)
  ynew <- apply(SG$design_unevaluated, 1, f)
  expect_error(SG <- SGGPfit(SG, Ynew=ynew, Xs=xsup, Ys=ysup), NA)
  
  # Check neglogpost grad matches gneglogpost for all HandlingSuppData options
  theta <- SG$thetaMAP / 2 # Don't want values near -1 or +1
  epsval <- 1e-4
  for (handling in c("Correct", "Only", "Ignore", "Mixture", "MarginalValidation", "FullValidation")) {
    thetagrad <- SGGP_internal_gneglogpost(theta, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)
    numgrad <- rep(0, length(SG$thetaMAP))
    for (i in 1:length(SG$thetaMAP)) {
      eps <- rep(0, length(SG$thetaMAP))
      eps[i] = eps[i] + epsval
      # numgrad <- (SGGP_internal_neglogpost(theta + eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) - 
      #               SGGP_internal_neglogpost(theta - eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)) / (2*epsval)
      numgrad[i] <- (-SGGP_internal_neglogpost(theta + 2*eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) + 
                    8*SGGP_internal_neglogpost(theta + eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) - 
                    8*SGGP_internal_neglogpost(theta - eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) + 
                    SGGP_internal_neglogpost(theta - 2*eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)) / (12*epsval)
      # 
      # print(numgrad)
    }
    expect_equal(c(thetagrad), numgrad, tol=1e-2, info = handling)
  }
  # numDeriv::grad(function(th) {SGGP_internal_neglogpost(th, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)}, theta)
})

test_that("SGGPfit works with Ynew - scalar output", {
  
  
  SG <- SGGPcreate(d=3, batchsize=20)
  f <- function(x){x[1]+x[2]^2}
  y <- apply(SG$design, 1, f)
  # Can't give in Y and Ynew
  expect_error(SG <- SGGPfit(SG, Y=y, Ynew=y))
  # Works with just Ynew
  SG <- SGGPfit(SG, Ynew=y)
  expect_true(all(!is.na(SG$thetaPostSamples)))
  
  # After append, only give in Ynew
  SG <- SGGPappend(SG, 50)
  ynew <- apply(SG$design_unevaluated, 1, f)
  # Again error if give in Y and Ynew
  expect_error(SGGPfit(SG, Y=c(y, ynew), Ynew=ynew))
  # Get error when Ynew is wrong size
  expect_error(SGGPfit(SG, Ynew=ynew[1:(length(ynew)-1)]))
  # Works fine
  expect_is(SG <- SGGPfit(SG, Ynew = ynew), "SGGP")
})

test_that("SGGPfit works with Ynew - vector output", {
  
  
  SG <- SGGPcreate(d=3, batchsize=20)
  f1 <- function(x){x[1]+x[2]^2}
  f2 <- function(x){x[2]^1.3+sin(2*pi*x[3])}
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y <- cbind(y1, y2)
  # Can't give in Y and Ynew
  expect_error(SG <- SGGPfit(SG, Y=y, Ynew=y))
  # Other errors when giving ynew and wrong number of rows
  expect_error(SGGPfit(SG, Ynew=y[1:5,]))
  expect_error(SGGPfit(SG, Ynew=y1[1:5]))
  
  # Works with just Ynew
  SG <- SGGPfit(SG, Ynew=y)
  
  # After append, only give in Ynew
  SG <- SGGPappend(SG, 50)
  ynew1 <- apply(SG$design_unevaluated, 1, f1)
  ynew2 <- apply(SG$design_unevaluated, 1, f2)
  ynew <- cbind(ynew1, ynew2)
  # Again error if give in Y and Ynew
  expect_error(SGGPfit(SG, Y=rbind(y, ynew), Ynew=ynew))
  # Get error when Ynew is wrong size
  expect_error(SGGPfit(SG, Ynew=ynew[1:(nrow(ynew)-1),]))
  # Get error if you give in vector
  expect_error(SGGPfit(SG, Ynew=ynew1))
  # Works fine
  expect_is(SG <- SGGPfit(SG, Ynew = ynew), "SGGP")
  
  
})

test_that("Using MCMC approx works", {
  f1 <- function(x){x[1]+x[2]^2}
  f2 <- function(x){x[2]^1.3+sin(2*pi*x[3])}
  
  # First 1D output
  SG <- SGGPcreate(d=3, batchsize=20)
  y1 <- apply(SG$design, 1, f1)
  SG$numPostSamples <- 5
  
  expect_error(SG <- SGGPfit(SG, Ynew=y1, laplaceapprox = F), NA)
  expect_equal(dim(SG$thetaPostSamples), c(3*SG$numpara,5))
  
  # All loglikelihoods should be finite for samples
  expect_true(all(is.finite(apply(SG$thetaPostSamples, 2, SGGP_internal_neglogpost, SG=SG, y=SG$y))))
  
  # Now with multiple output
  SG <- SGGPcreate(d=3, batchsize=20)
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y <- cbind(y1, y2)
  
  # speed up test by reducing number of thetaPostSamples
  SG$numPostSamples <- 3
  
  SG1 <- SGGPfit(SG, Ynew=y, laplaceapprox = F)
  SG2 <- SGGPfit(SG, Ynew=y, laplaceapprox = T)
  expect_equal(dim(SG1$thetaPostSamples), dim(SG2$thetaPostSamples))
  
  SG1 <- SGGPfit(SG, Ynew=y, laplaceapprox = F, separateoutputparameterdimensions = T)
  SG2 <- SGGPfit(SG, Ynew=y, laplaceapprox = T, separateoutputparameterdimensions = T)
  expect_equal(dim(SG1$thetaPostSamples), dim(SG2$thetaPostSamples))
  if (F) {
    boxplot(SG1$thetaPostSamples[1,,1], SG2$thetaPostSamples[1,,1])
  }
})

test_that("postvarmatcalc", {
  # These don't work unless x1 and x2 have same length.
  # Maybe only makes sense when x1 == x2?
  x1 <- runif(5)
  x2 <- x1
  o1 <- SGGP_internal_postvarmatcalc(x1, x2,
                               xo=c(.11), theta=c(.1,.2,.3),
                               CorrMat=SGGP_internal_CorrMatCauchySQT,
                               returndPVMC=F, returndiagonly=F)
  o2 <- SGGP_internal_postvarmatcalc(x1, x2,
                                     xo=c(.11), theta=c(.1,.2,.3),
                                     CorrMat=SGGP_internal_CorrMatCauchySQT,
                                     returndPVMC=F, returndiagonly=T)
  o3 <- SGGP_internal_postvarmatcalc(x1, x2,
                                     xo=c(.11), theta=c(.1,.2,.3),
                                     CorrMat=SGGP_internal_CorrMatCauchySQT,
                                     returndPVMC=T, returndiagonly=F)
  o4 <- SGGP_internal_postvarmatcalc(x1, x2,
                                     xo=c(.11), theta=c(.1,.2,.3),
                                     CorrMat=SGGP_internal_CorrMatCauchySQT,
                                     returndPVMC=T, returndiagonly=T)
  expect_equal(o1, o3$Sigma_mat)
  expect_equal(diag(o1), o2)
  expect_equal(o4$Sigma_mat, diag(o3$Sigma_mat))
  expect_equal(o4$dSigma_mat[,1], diag(o3$dSigma_mat[,1:5]))
  expect_equal(o4$dSigma_mat[,2], diag(o3$dSigma_mat[,6:10]))
  expect_equal(o4$dSigma_mat[,3], diag(o3$dSigma_mat[,11:15]))
})