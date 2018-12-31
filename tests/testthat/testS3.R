test_that("print and predict work", {
  
  SG = SGGPcreate(3,201)
  tp <- capture.output(print(SG))
  expect_is(tp, "character")
  expect_gt(length(tp), 5)
})