test_that("Validation works", {
  set.seed(0)
  
  SG <- SGcreate(d=3, batchsize=100)
  f1 <- function(x){x[1]+x[2]^2+rnorm(1,0,.01)}
  y <- apply(SG$design, 1, f1)
  Xval <- matrix(runif(3*100), ncol=3)
  Yval <- apply(Xval, 1, f1)
  val <- validation(c(.1,.1,.1), SG=SG, y=y, xval=Xval, yval=Yval)
  expect_is(val, "numeric")
  expect_equal(length(val), 1)
  
  # Make sure logthetaVALID returns best
  SG <- logthetaVALID(SG=SG, y=y, xval=Xval, yval=Yval)
  expect_lt(validation(SG$logtheta, SG=SG, y=y, xval=Xval, yval=Yval), 
            validation(SG$logtheta+c(.1,.1,.1), SG=SG, y=y, xval=Xval, yval=Yval))
  
  # Check gvalidation is grad of validation
  th <- c(.1,.2,.3)
  g <- gvalidation(logtheta = th, SG=SG, y=y, xval=Xval, yval=Yval)
  theps <- 1e-2
  expect_equal(g[1],
               (validation(th+c(theps/2,0,0), SG=SG, y=y, xval=Xval, yval=Yval) - 
                  validation(th-c(theps/2,0,0), SG=SG, y=y, xval=Xval, yval=Yval)) / theps,
               tolerance=1e-4)
  expect_equal(g[2],
               (validation(th+c(0,theps/2,0), SG=SG, y=y, xval=Xval, yval=Yval) - 
                  validation(th-c(0,theps/2,0), SG=SG, y=y, xval=Xval, yval=Yval)) / theps,
               tolerance=1e-4)
  expect_equal(g[3],
               (validation(th+c(0,0,theps/2), SG=SG, y=y, xval=Xval, yval=Yval) - 
                  validation(th-c(0,0,theps/2), SG=SG, y=y, xval=Xval, yval=Yval)) / theps,
               tolerance=1e-4)
  
  # Make sure big logtheta gives error
  expect_equal(validation(logtheta=c(4.1,.1,.1), SG=SG, y=y, xval=Xval, yval=Yval), Inf)
  
  # Make sure return_optim works
  expect_is(SG, "SGGP")
  expect_is(SG, "list")
  optout <- logthetaVALID(SG=SG, y=y, xval=Xval, yval=Yval, return_optim=T)
  expect_is(optout, "list")
  expect_equal(length(optout), 5)
})