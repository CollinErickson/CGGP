test_that("print and predict work", {
  
  SG = SGGPcreate(3,201, corr="cauchy")
  tp <- capture.output(print(SG))
  expect_is(tp, "character")
  expect_gt(length(tp), 5)
})