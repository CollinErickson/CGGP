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
    num_dpw <- (SGGP_internal_calcpw(      SG=SG, y=y, theta=SG$thetaMAP + delta) - 
                  SGGP_internal_calcpw(      SG=SG, y=y, theta=SG$thetaMAP - delta)) / 
      (2*eps)
    expect_equal(num_dpw, pwanddpw$dpw[,i], 1e-4)
  }
})