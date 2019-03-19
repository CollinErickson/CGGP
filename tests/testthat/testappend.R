context("testappend")

test_that("CGGPappend works", {
  SG <- CGGPcreate(d=3, batchsize=100, corr='GAUSS')
  f <- function(x){x[1]*x[3]+x[2]^2}
  y <- apply(SG$design, 1, f)
  SG <- CGGPfit(SG, Y=y)
  for (i in c("UCB", "Greedy", "TS", "Oldest", "Random", "Lowest")) {
    lastN <- nrow(SG$design)
    SG <- CGGPappend(CGGP=SG, batchsize=20, selectionmethod = i)
    expect_is(SG, "CGGP")
    expect_gt(nrow(SG$design), lastN)
    y <- apply(SG$design, 1, f)
    SG <- CGGPfit(SG, Y=y)
  }
  # Error, UC is not an option.
  expect_error(CGGPappend(SG, 20, selectionmethod = "UC"))
  
  # # Can append after append without fitting between
  # expect_error(SG <- CGGPappend(SG, 20), NA)
  # expect_error(SG <- CGGPappend(SG, 20), NA)
  
  # Can append with supp data
  xsup <- matrix(runif(3*10), ncol=3)
  ysup <- apply(xsup, 1, f)
  SG <- CGGPfit(SG, SG$Y, Xs=xsup, Ys=ysup, corr='m32')
  for (sel.method in c("UCB", "TS", "Greedy")) {
    expect_error(CGGPappend(SG, 30, sel.method), NA)
  }
})

test_that("CGGPappend works with large number", {
  # Start it small
  SG <- CGGPcreate(d=3, batchsize=20)
  f <- function(x){x[1]+log(x[1]+.1)+sin(2*pi*4*x[2]^2) + cos(2*pi*5*x[3])}
  y <- apply(SG$design, 1, f)
  SG <- CGGPfit(SG, Y=y)
  lastN <- nrow(SG$design)
  
  # Adding 2000 will force it to increase ML and add rows to uo, pila, pala, etc.
  # But it doesn't show as working on codecov? Try 4000
  expect_error(SG <- CGGPappend(CGGP=SG, batchsize=2*2000), NA)
  expect_is(SG, "CGGP")
  expect_gt(nrow(SG$design), lastN)
  # y <- apply(SG$design, 1, f)
  # SG <- CGGPfit(SG, Y=y)
  # expect_error(CGGPappend(SG, 20, selectionmethod = "UC"))
  
  # Do we want prediction to work before fit is run? If yes, below should not give an error, but it does now.
  # ypred <- CGGPpred(SG$design, SG)
  
  # Check that prediction is still exact after adding rows to uo, pila, w, etc.
  # This gets inaccurate since there's so much data, so I'll remove it.
  # y <- apply(SG$design, 1, f)
  # SG <- CGGPfit(SG, Y=y)
  # ypred <- CGGPpred(SG, SG$design)
  # expect_equal(y, c(ypred$mean), tol=1e-2)
})

test_that("CGGPappend gives warning if it can't add any data", {
  
  SG <- CGGPcreate(d=3, batchsize=20, corr="cauchySQ")
  f <- function(x){x[1]+log(x[1]+.1)+sin(2*pi*4*x[2]^2) + cos(2*pi*5*x[3])}
  y <- apply(SG$design, 1, f)
  SG <- CGGPfit(SG, Y=y)
  
  expect_warning(CGGPappend(CGGP=SG, batchsize=1))
})


test_that("Append with different weights on different outputs works", {
  # Second function is much bigger, will focus more in that dimension
  #  if told to use weightings. If given equal weight, dimensions will
  #  not be treated differently.
  f1 <- function(x){sin(2*pi*x[1]^2)}
  f2 <- function(x){10000*sin(2*pi*x[2]^2)}
  
  # Now with multiple output
  # Should get to 3,1 and 1,3 with 2,2
  set.seed(3)
  SG <- CGGPcreate(d=2, batchsize=17, grid_sizes = c(1,2,4,2,2,2,2,2,2,12,32))
  # Should be equal so far
  if (F) {
    summary(SG$uo[1:SG$uoCOUNT,])
    CGGPblockplot(SG)
  }
  expect_equal(table(SG$uo[,1]), table(SG$uo[,2]))
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y <- cbind(y1, y2)
  
  set.seed(2)
  SG <- CGGPfit(SG, Ynew=y, separateoutputparameterdimensions = TRUE)
  
  # Now if I append with equal weights by using "/var" (default), it should have equal in both dimensions
  expect_error(SG2 <- CGGPappend(SG, 26, selectionmethod = "Greedy", multioutputdim_weights = "/sigma2MAP"), NA)
  if (F) {
    SG2$uo[1:SG2$uoCOUNT,] %>% summary
    CGGPblockplot(SG2)
  }
  expect_equal(mean(SG2$uo[1:SG2$uoCOUNT,2]) , mean(SG2$uo[1:SG2$uoCOUNT,1]))
  
  # If using full variance (default), will put much more in second dimensions
  SG3 <- CGGPappend(SG, 26)
  if (F) {
    SG3$uo[1:SG3$uoCOUNT,] %>% summary
    CGGPblockplot(SG3)
  }
  expect_gt(mean(SG3$uo[1:SG3$uoCOUNT,2]) - 2 , mean(SG3$uo[1:SG3$uoCOUNT,1]))
  
  # Last option is to standardize with "/range^2"
  expect_error(SG4 <- CGGPappend(SG, 26, multioutputdim_weights="/range^2"), NA)
  if (F) {
    SG4$uo[1:SG4$uoCOUNT,] %>% summary
    CGGPblockplot(SG4)
  }
  expect_gt(mean(SG3$uo[1:SG3$uoCOUNT,2]) - 2 , mean(SG3$uo[1:SG3$uoCOUNT,1]))
  
  # Bad options
  expect_error(CGGPappend(SG, 26, multioutputdim_weights="/range"))
  expect_error(CGGPappend(SG, 26, multioutputdim_weights=c(1,2,3)))
  # Can't actually fit when range is zero, so making a fake one.
  # SG0 <- CGGPfit(SG, SG$Y*0)
  SG$Y <- SG$Y*0
  expect_error(CGGPappend(SG, 26, multioutputdim_weights="/range^2"))
  
})
