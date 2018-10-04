test_that("Likelihood matches dLLH for Matern 3/2", {
  
  SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
  y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  th <- c(.1,.1,.1)
  f <- lik(th, SG=SG, y=y)
  expect_equal(dim(f), c(1,1))
  g <- glik(th, SG=SG, y=y)
  expect_equal(dim(g), c(1,3))
  theps <- 1e-5
  expect_equal(g[1], (lik(th+c(theps/2,0,0), SG=SG, y=y) - lik(th-c(theps/2,0,0), SG=SG, y=y))[1,1] / theps, tolerance=1e-4)
  expect_equal(g[2], (lik(th+c(0,theps/2,0), SG=SG, y=y) - lik(th-c(0,theps/2,0), SG=SG, y=y))[1,1] / theps, tolerance=1e-4)
  expect_equal(g[3], (lik(th+c(0,0,theps/2), SG=SG, y=y) - lik(th-c(0,0,theps/2), SG=SG, y=y))[1,1] / theps, tolerance=1e-4)
})
