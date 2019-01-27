test_that("Plots work", {
  
  SG <- SGGPcreate(d=3, batchsize=20, corr="m52")
  f <- function(x){x[1]+x[2]*log(x[1]+.1)+sin(2*pi*4*x[2]^2) + sqrt(x[1])*cos(2*pi*5*x[3])}
  y <- apply(SG$design, 1, f)
  # SG <- SGGPfit(SG, Y=y)
  expect_error(SG <- SGGPfit(SG, Y=y), NA)
  
  # Heat map
  pheat <- SGGPheat(SG)
  expect_is(pheat, "gg")
  expect_is(pheat, "ggplot")
  
  # Histogram
  phist <- SGGPhist(SG)
  expect_is(phist, "gg")
  expect_is(phist, "ggplot")
  # # Don't actually want warning for log scale, but getting it, so expect it
  # expect_warning(print(phist))
  
  # Blockplot
  expect_error(bp <- SGGPblockplot(SG), NA)
  expect_is(bp, "ggplot")
  
  # Validation plot
  Xval <- matrix(runif(3*100), ncol=3)
  Yval <- apply(Xval, 1, f)
  pval <- SGGPvalplot(SGGP=SG, Xval=Xval, Yval=Yval)
  expect_is(pval, "gg")
  expect_is(pval, "ggplot")
  expect_error(SGGPvalplot(SGGP=SG, Xval=Xval, Yval=Yval, plot_with = 'base'), NA)
  
  # Val stats
  vstats <- valstats(c(0,1,2), c(.01,.01,.01), c(0,1.1,1.9))
  expect_is(vstats, "data.frame")
  expect_equal(nrow(vstats), 1)
  rm(vstats)
  
  # SGGP Val stats
  vstats <- SGGPvalstats(SG, Xval, Yval)
  expect_is(vstats, "data.frame")
  expect_equal(nrow(vstats), 1)
  
  # Multiple outputs
  SG2 <- SGGPcreate(d=3, batchsize=100)
  f1 <- function(x){x[1]+x[2]^2/(x[1]+1) +sin(x[1]*x[3])}
  f2 <- function(x){x[1]^1.3+.4*sin(2*pi*x[2])*x[1] + (x[1]+1)*exp(x[3])+10}
  y1 <- apply(SG2$design, 1, f1)#+rnorm(1,0,.01)
  y2 <- apply(SG2$design, 1, f2)#+rnorm(1,0,.01)
  y <- cbind(y1, y2)
  SG2 <- SGGPfit(SG2, Y=y)
  Xval <- matrix(runif(3*100), ncol=3)
  Yval <- cbind(apply(Xval, 1, f1),
                apply(Xval, 1, f2))
  sv <- SGGPvalstats(SG2, Xval, Yval)
  expect_is(sv, "data.frame")
  expect_equal(nrow(sv), 2)
  # expect_true(all(sv$RMSE<.1))
  expect_equal(nrow(SGGPvalstats(SG2, Xval, Yval, bydim=FALSE)), 1)
  expect_error(SGGPvalstats(SG2, Xval, Yval[,1]))
  
  # Corr plot
  expect_error(p <- SGGPcorrplot(), NA)
  expect_is(p, "ggplot")
  expect_error(p <- SGGPcorrplot(SGGP_internal_CorrMatCauchySQ, theta=c(-.9,.8,.7,-.8)), NA)
  expect_is(p, "ggplot")
  expect_error(p <- SGGPcorrplot(SG), NA)
  expect_is(p, "ggplot")
  expect_error(p <- SGGPcorrplot(SG2), NA)
  expect_is(p, "ggplot")
  expect_error(SGGPcorrplot(SG, plot_with = "base"), NA)
  rm(p)
  
  # Projection plot
  # These should work fine, return ggplot
  expect_error(pp1 <- SGGPprojectionplot(SG), NA)
  expect_is(pp1, "ggplot")
  expect_error(pp1 <- SGGPprojectionplot(SG2), NA)
  expect_is(pp1, "ggplot")
  expect_error(pp2 <- SGGPprojectionplot(SG, proj = c(0)), NA)
  expect_is(pp2, "ggplot")
  # Error if proj is not 1 or d dim
  expect_error(SGGPprojectionplot(SG, proj=c(1,1)))
  expect_error(SGGPprojectionplot(SG, proj=c(1,1,1,1)))
})
