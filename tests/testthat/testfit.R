# SGGPfit is used in many other tests, so not as many here

test_that("SGGPfit works", {
  
  SG <- SGGPcreate(d=3, batchsize=20)
  f <- function(x){x[1]+x[2]^2}
  y <- apply(SG$design, 1, f)
  SG <- SGGPfit(SG, Y=y)
  
  neglogpost <- SGGP_internal_gneglogpost(rep(.5,9), SG, y)
  expect_is(neglogpost, "matrix")
  expect_length(neglogpost, 9)
  
  
  neglogpost2 <- SGGP_internal_gneglogpost(rep(.5,9), SG, y, return_lik = T)
  expect_is(neglogpost2, "list")
  expect_length(neglogpost2[[1]], 1)
  expect_length(neglogpost2[[2]], 9)
  
  # # Works with supplementary data
  # nsup <- 10
  # xsup <- matrix(runif(nsup*3), nsup, 3)
  # ysup <- apply(xsup, 1, f)
  # SG <- SGGPappend(SG, 20)
  # ynew <- apply(SG$design_unevaluated, 1, f)
  # expect_error(SG <- SGGPfit(SG, Ynew=ynew, Xs=xsup, Ys=ysup), NA)
  
})

test_that("SGGPfit works with Ynew - scalar output", {
  
  
  SG <- SGGPcreate(d=3, batchsize=20)
  f <- function(x){x[1]+x[2]^2}
  y <- apply(SG$design, 1, f)
  # Can't give in Y and Ynew
  expect_error(SG <- SGGPfit(SG, Y=y, Ynew=y))
  # Works with just Ynew
  SG <- SGGPfit(SG, Ynew=y)
  
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

test_that("Not using LaPlace approx works", {
  f1 <- function(x){x[1]+x[2]^2}
  f2 <- function(x){x[2]^1.3+sin(2*pi*x[3])}
  
  # First 1D output
  SG <- SGGPcreate(d=3, batchsize=20)
  y1 <- apply(SG$design, 1, f1)
  
  expect_error(SG <- SGGPfit(SG, Ynew=y1, laplaceapprox = F), NA)
  
  # Now with multiple output
  SG <- SGGPcreate(d=3, batchsize=20)
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y <- cbind(y1, y2)
  
  # speed up test by reducing number of thetaPostSamples
  SG$numPostSamples <- 10
  
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
