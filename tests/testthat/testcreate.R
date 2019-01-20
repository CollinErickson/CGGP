test_that("SGGPcreate works", {
  expect_error(SG <- SGGPcreate())
  expect_error(SG <- SGGPcreate(d=3))
  expect_error(SG <- SGGPcreate(batchsize=20))
  expect_error(SG <- SGGPcreate(d=1, batchsize=20))
  
  
  SG <- SGGPcreate(d=3, batchsize=20, corr = "MATern32")
  expect_is(SG, "list")
  expect_is(SG, "SGGP")
  
  # Create a big one that forces it to allocate more memory
  expect_error(SGGPcreate(d=4, 7000, corr="powerEXP"), NA)
})
