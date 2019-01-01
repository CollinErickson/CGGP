test_that("SGGPappend works", {
  SG <- SGGPcreate(d=3, batchsize=100)
  f <- function(x){x[1]+x[2]^2}
  y <- apply(SG$design, 1, f)
  SG <- SGGPfit(SG, Y=y)
  for (i in c("UC", "Greedy", "TS")) {
    lastN <- nrow(SG$design)
    SG <- SGGPappend(SGGP=SG, batchsize=20)
    expect_is(SG, "SGGP")
    expect_gt(nrow(SG$design), lastN)
    y <- apply(SG$design, 1, f)
    SG <- SGGPfit(SG, Y=y)
  }
})