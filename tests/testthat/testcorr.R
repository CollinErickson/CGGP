test_that("Gaussian correlation works", {
  # expect_is for class checking 
  # expect_equal for numeric
  x1 <- runif(5)
  x2 <- runif(4)
  th <- .1
  gc <- gausscorr(x1, x2, th)
  expect_equal(gc[2,3], exp(-(x1[2]-x2[3])^2/th))
  expect_equal(dim(gc), c(5,4))
  
  dgc <- dgausscorr(x1, x2, th)
  expect_equal(dim(dgc), c(5,4))
  teps <- 1e-5
  expect_equal(dgc, (gausscorr(x1, x2,th+teps/2) - gausscorr(x1, x2,th-teps/2))/teps)
  rm(list=ls())
})

test_that("Matern 3/2 correlation works", {
  # Note that it doesn't have sqrt(3), it's absorbed into theta unfortunately
  # (1 + diffmat / exp(theta)) * exp(-diffmat / exp(theta))
  x1 <- runif(5)
  x2 <- runif(4)
  th <- .1
  lth <- log(th)
  mc <- CorrMat(x1, x2, lth)
  expect_equal(mc[2,3], (1+abs(x1[2]-x2[3])/th) * exp(-abs(x1[2]-x2[3])/th))
  expect_equal(dim(mc), c(5,4))
  
  dmc <- dCorrMat(x1, x2, lth)
  expect_equal(dim(dmc), c(5,4))
  teps <- 1e-5
  expect_equal(dmc, (CorrMat(x1, x2,lth+teps/2) - CorrMat(x1, x2,lth-teps/2))/teps)
})
