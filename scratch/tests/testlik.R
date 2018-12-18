test_that("Likelihood matches dLLH for Matern 3/2", {
  
  SG <- SGcreate(d=3, batchsize=100)
  y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  th <- c(.1,.1,.1)
  f <- lik(th, SG=SG, y=y)
  expect_length(f, 1)
  g <- glik(th, SG=SG, y=y)
  expect_is(g, "numeric")
  expect_equal(length(g), 3)
  theps <- 1e-5
  expect_equal(g[1], (lik(th+c(theps/2,0,0), SG=SG, y=y) - lik(th-c(theps/2,0,0), SG=SG, y=y)) / theps, tolerance=1e-4)
  expect_equal(g[2], (lik(th+c(0,theps/2,0), SG=SG, y=y) - lik(th-c(0,theps/2,0), SG=SG, y=y)) / theps, tolerance=1e-4)
  expect_equal(g[3], (lik(th+c(0,0,theps/2), SG=SG, y=y) - lik(th-c(0,0,theps/2), SG=SG, y=y)) / theps, tolerance=1e-4)
  
})

test_that("More lik tests", {
  SG <- SGcreate(d=3, batchsize=100)
  y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  
  # Check inf for large
  expect_equal(lik(logtheta=c(5,6,7), SG=SG, y=y), Inf)
  
  out <- logthetaMLE(SG=SG, y=y)
  expect_is(out, "SGGP")
  
  # Check it actually is minimum
  expect_lt(lik(logtheta=out$logtheta, SG=SG, y = y),
            lik(logtheta=c(0,0,0), SG=SG, y = y)) # Perturbation is worse
  
  out2 <- logthetaMLE(SG=SG, y=y, return_optim = T)
  expect_is(out2, "list")
  expect_equal(length(out2), 5)
})

test_that("MV Likelihood matches dLLH for Matern 3/2", {
  
  SG <- SGcreate(d=3, batchsize=100)
  y1 <- apply(SG$design, 1, function(x){exp(x[1])+x[2]^.2+rnorm(1,0,.01)})
  y2 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  y <- cbind(y1, y2)
  th <- c(.1,.2,.3)
  f <- likMV(th, SG=SG, y=y)
  expect_length(f, 1)
  g <- glikMV(th, SG=SG, y=y)
  expect_is(g, "numeric")
  expect_equal(length(g), 3)
  theps <- 1e-5
  expect_equal(g[1], (likMV(th+c(theps/2,0,0), SG=SG, y=y) - likMV(th-c(theps/2,0,0), SG=SG, y=y)) / theps, tolerance=1e-4)
  expect_equal(g[2], (likMV(th+c(0,theps/2,0), SG=SG, y=y) - likMV(th-c(0,theps/2,0), SG=SG, y=y)) / theps, tolerance=1e-4)
  expect_equal(g[3], (likMV(th+c(0,0,theps/2), SG=SG, y=y) - likMV(th-c(0,0,theps/2), SG=SG, y=y)) / theps, tolerance=1e-4)
  
})


test_that("More MV lik tests", {
  SG <- SGcreate(d=3, batchsize=100)
  y1 <- apply(SG$design, 1, function(x){exp(x[1])+x[2]^.2+rnorm(1,0,.01)})
  y2 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  y <- cbind(y1, y2)
  
  # Check inf for large
  expect_equal(likMV(logtheta=c(5,6,7), SG=SG, y=y), Inf)
  
  # Test logthetaMLEMV
  
  # Returns correct type
  out <- logthetaMLEMV(SG=SG, y=y)
  expect_is(out, "SGGP")
  
  # Check it actually is minimum
  expect_lt(likMV(logtheta=out$logtheta, SG=SG, yMV = y),
                   likMV(logtheta=out$logtheta+.02, SG=SG, yMV = y)) # Perturbation is worse
  
  # It can return just optim output
  out2 <- logthetaMLEMV(SG=SG, y=y, return_optim = T)
  expect_is(out2, "list")
  expect_equal(length(out2), 5)
})