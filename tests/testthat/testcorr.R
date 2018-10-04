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
})