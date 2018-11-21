test_that("Append appends something", {
  
  SG <- SGcreate(d=3, batchsize=100)
  SG <- SGappend(theta=c(.1,.1,.1), SG=SG, batchsize=20)
  expect_gte(nrow(SG$design), 110)
  expect_lte(nrow(SG$design), 120)
})

test_that("Append is smart", {
  SG <- SGcreate(d=5, batchsize=10)
  # Should put more in first dimension
  SG <- SGappend(logtheta=c(-2,2,2,2,2), SG=SG, batchsize=100)
  
  # First dim should have gone into further level
  expect_gt(max(SG$uo[1:SG$ss,1]), max(SG$uo[1:SG$ss,2:5]))
  
  # And have higher average level
  expect_gt(mean(SG$uo[1:SG$ss,1]), max(apply(SG$uo[1:SG$ss,2:5], 2, mean)))
})