test_that("SGGPcreate works", {
  expect_error(SG <- SGGPcreate())
  expect_error(SG <- SGGPcreate(d=3))
  expect_error(SG <- SGGPcreate(batchsize=20))
  expect_error(SG <- SGGPcreate(d=1, batchsize=20))
  
  
  SG <- SGGPcreate(d=3, batchsize=20)
  expect_is(SG, "list")
  expect_is(SG, "SGGP")
})