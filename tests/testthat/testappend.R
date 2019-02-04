test_that("SGGPappend works", {
  SG <- SGGPcreate(d=3, batchsize=100, corr='GAUSS')
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
  # Error, UC is not an option.
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
  expect_error(SG <- SGGPappend(SGGP=SG, batchsize=2000), NA)
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
  ypred <- SGGPpred(SG, SG$design)
  expect_equal(y, c(ypred$mean), tol=1e-2)
})

test_that("SGGPappend gives warning if it can't add any data", {
  
  SG <- SGGPcreate(d=3, batchsize=20, corr="cauchySQ")
  f <- function(x){x[1]+log(x[1]+.1)+sin(2*pi*4*x[2]^2) + cos(2*pi*5*x[3])}
  y <- apply(SG$design, 1, f)
  SG <- SGGPfit(SG, Y=y)
  
  expect_warning(SGGPappend(SGGP=SG, batchsize=1))
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
  SG <- SGGPcreate(d=2, batchsize=17, grid_sizes = c(1,2,4,2,2,2,2,2,2,12,32))
  # Should be equal so far
  if (F) {
    summary(SG$uo[1:SG$uoCOUNT,])
    SGGPblockplot(SG)
  }
  expect_equal(table(SG$uo[,1]), table(SG$uo[,2]))
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y <- cbind(y1, y2)
  
  set.seed(2)
  SG <- SGGPfit(SG, Ynew=y, use_PCA = FALSE, separateoutputparameterdimensions = TRUE)
  
  # Now if I append with equal weights by using "/var" (default), it should have equal in both dimensions
  SG2 <- SGGPappend(SG, 26, selectionmethod = "Greedy", multioutputdim_weights = "/sigma2MAP")
  if (F) {
    SG2$uo[1:SG2$uoCOUNT,] %>% summary
    SGGPblockplot(SG2)
  }
  expect_equal(mean(SG2$uo[1:SG2$uoCOUNT,2]) , mean(SG2$uo[1:SG2$uoCOUNT,1]))
  
  # If using full variance (default), will put much more in second dimensions
  SG3 <- SGGPappend(SG, 26)
  if (F) {
    SG3$uo[1:SG3$uoCOUNT,] %>% summary
    SGGPblockplot(SG3)
  }
  expect_gt(mean(SG3$uo[1:SG3$uoCOUNT,2]) - 2 , mean(SG3$uo[1:SG3$uoCOUNT,1]))
  
  # Last option is to standardize with "/range^2"
  SG4 <- SGGPappend(SG, 26, multioutputdim_weights="/range^2")
  if (F) {
    SG4$uo[1:SG4$uoCOUNT,] %>% summary
    SGGPblockplot(SG4)
  }
  expect_gt(mean(SG3$uo[1:SG3$uoCOUNT,2]) - 2 , mean(SG3$uo[1:SG3$uoCOUNT,1]))
})
