test_that("Plots work", {
  
  SG <- SGGPcreate(d=3, batchsize=20, corr="m52")
  f <- function(x){x[1]+log(x[1]+.1)+sin(2*pi*4*x[2]^2) + cos(2*pi*5*x[3])}
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
  
  # Validation plot
  Xval <- matrix(runif(3*100), ncol=3)
  Yval <- apply(Xval, 1, f)
  pval <- SGGPvalplot(SGGP=SG, Xval=Xval, Yval=Yval)
  expect_is(pval, "gg")
  expect_is(pval, "ggplot")
  expect_error(SGGPvalplot(SGGP=SG, Xval=Xval, Yval=Yval, plot_with = 'base'), NA)
  
  # Val stats
  vstats <- SGGPvalstats(SG, Xval, Yval)
  expect_is(vstats, "data.frame")
  expect_equal(nrow(vstats), 1)
  
  # Corr plot
  p <- SGGP_internal_CorrPlot(SGGP_internal_CorrMatCauchySQ, theta=c(-.9,.8,.7,-.8))
  expect_is(p, "ggplot")
  p <- SGGP_internal_CorrPlot(SG)
  expect_is(p, "ggplot")
  expect_error(SGGP_internal_CorrPlot(SG, plot_with = "base"), NA)
})
