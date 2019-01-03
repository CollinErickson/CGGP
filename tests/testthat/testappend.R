test_that("SGGPappend works", {
  SG <- SGGPcreate(d=3, batchsize=100)
  f <- function(x){x[1]+x[2]^2}
  y <- apply(SG$design, 1, f)
  SG <- SGGPfit(SG, Y=y)
  for (i in c("UCB", "Greedy", "TS")) {
    lastN <- nrow(SG$design)
    SG <- SGGPappend(SGGP=SG, batchsize=20, selectionmethod = i)
    expect_is(SG, "SGGP")
    expect_gt(nrow(SG$design), lastN)
    y <- apply(SG$design, 1, f)
    SG <- SGGPfit(SG, Y=y)
  }
  expect_error(SGGPappend(SG, 20, selectionmethod = "UC"))
})

test_that("SGGPappend works with large number", {
  # Start it small
  SG <- SGGPcreate(d=3, batchsize=20)
  f <- function(x){x[1]+log(x[1]+.1)+sin(2*pi*4*x[2]^2) + cos(2*pi*5*x[3])}
  y <- apply(SG$design, 1, f)
  SG <- SGGPfit(SG, Y=y)
  lastN <- nrow(SG$design)
  
  # Adding 2000 will force it to increase ML and add rows to uo, pila, pala, etc.
  expect_silent(SG <- SGGPappend(SGGP=SG, batchsize=2000))
  expect_is(SG, "SGGP")
  expect_gt(nrow(SG$design), lastN)
  # y <- apply(SG$design, 1, f)
  # SG <- SGGPfit(SG, Y=y)
  # expect_error(SGGPappend(SG, 20, selectionmethod = "UC"))
  
  # Do we want prediction to work before fit is run? If yes, below should not give an error, but it does now.
  # ypred <- SGGPpred(SG$design, SG)
  
  # Check that prediction is still exact after adding rows to uo, pila, w, etc.
  y <- apply(SG$design, 1, f)
  SG <- SGGPfit(SG, Y=y)
  ypred <- SGGPpred(SG$design, SG)
  expect_equal(y, c(ypred$mean), tol=1e-4)
})

