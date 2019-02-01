test_that("calcpw", {
  SG <- SGGPcreate(d=3, batchsize=100)
  y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  pw       <- SGGP_internal_calcpw(      SG=SG, y=y, theta=SG$thetaMAP)
  pwanddpw <- SGGP_internal_calcpwanddpw(SG=SG, y=y, theta=SG$thetaMAP)
  expect_equal(pw, pwanddpw$pw)
  
  # Check dpw with numerical grad
  eps <- 1e-4
  for (i in 1:9) {
    delta <- rep(0,9)
    delta[i] <- eps
    num_dpw <- (SGGP_internal_calcpw(SG=SG, y=y, theta=SG$thetaMAP + delta) - 
                  SGGP_internal_calcpw(SG=SG, y=y, theta=SG$thetaMAP - delta)) / 
      (2*eps)
    expect_equal(num_dpw, pwanddpw$dpw[,i], 1e-4)
  }
  
  nr <- nrow(SG$design)
  nb <- max(SG$uo[1:SG$uoCOUNT,])
  # Get list back when ask for return_lS
  pw_with_lS <- SGGP_internal_calcpw(SG=SG, y=y, theta=SG$thetaMAP, return_lS = T)
  expect_is(pw_with_lS, "list")
  expect_length(pw_with_lS, 2)
  expect_equal(names(pw_with_lS), c("pw", "lS"))
  expect_is(pw_with_lS$pw, "numeric")
  expect_length(pw_with_lS$pw, nr)
  expect_is(pw_with_lS$lS, "matrix")
  expect_equal(dim(pw_with_lS$lS), c(nb,3))
  rm(pw_with_lS)
  # Other case when y is matrix
  pw_with_lS2 <- SGGP_internal_calcpw(SG=SG, y=cbind(y,y), theta=SG$thetaMAP, return_lS = T)
  expect_is(pw_with_lS2, "list")
  expect_length(pw_with_lS2, 2)
  expect_equal(names(pw_with_lS2), c("pw", "lS"))
  expect_is(pw_with_lS2$pw, "matrix")
  expect_equal(dim(pw_with_lS2$pw), c(nr,2))
  expect_is(pw_with_lS2$lS, "matrix")
  expect_equal(dim(pw_with_lS2$lS), c(nb,3))
  rm(pw_with_lS2)
  
  # Same but now with dw and dpw
  # Get list back when ask for return_lS
  pw_with_lS <- SGGP_internal_calcpwanddpw(SG=SG, y=y, theta=SG$thetaMAP, return_lS = T)
  expect_is(pw_with_lS, "list")
  expect_length(pw_with_lS, 4)
  expect_equal(names(pw_with_lS), c("pw", "dpw", "lS", "dlS"))
  expect_is(pw_with_lS$pw, "numeric")
  expect_length(pw_with_lS$pw, nr)
  expect_is(pw_with_lS$dpw, "matrix")
  expect_equal(dim(pw_with_lS$dpw), c(nr,9))
  expect_is(pw_with_lS$lS, "matrix")
  expect_equal(dim(pw_with_lS$lS), c(nb,3))
  expect_is(pw_with_lS$dlS, "matrix")
  expect_equal(dim(pw_with_lS$dlS), c(nb,9))
  rm(pw_with_lS)
  # Other case when y is matrix
  pw_with_lS <- SGGP_internal_calcpwanddpw(SG=SG, y=cbind(y,y), theta=SG$thetaMAP, return_lS = T)
  expect_is(pw_with_lS, "list")
  expect_length(pw_with_lS, 4)
  expect_equal(names(pw_with_lS), c("pw", "dpw", "lS", "dlS"))
  expect_is(pw_with_lS$pw, "numeric")
  expect_length(pw_with_lS$pw, (nr*2))
  expect_is(pw_with_lS$dpw, "matrix")
  expect_equal(dim(pw_with_lS$dpw), c(nr*2,9))
  expect_is(pw_with_lS$lS, "matrix")
  expect_equal(dim(pw_with_lS$lS), c(nb,3))
  expect_is(pw_with_lS$dlS, "matrix")
  expect_equal(dim(pw_with_lS$dlS), c(nb,9))
  rm(pw_with_lS)
})

test_that("calcpw chol error", {
  
  SG <- SGGPcreate(d=3, batchsize=100, corr = 'gauss')
  y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  SG <- SGGPfit(SG, y)
  # # Force error when theta is small
  # expect_error(SGGP_internal_calcpw(SG, y, -400))
  # Not an error, but Inf. The error is caught, so no real error is returned
  expect_equal(SGGP_internal_calcpw(SG, y, -400), Inf)
})

