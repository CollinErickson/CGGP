test_that("1. Create, append, predict with only supp, scalar out", {
  d <- 3
  # f <- function(x){1*(cos(x[3]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + log(x[3]+.2))}
  f <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  
  set.seed(0)
  nsup <- 30
  xsup <- matrix(runif(nsup*d), nsup, d)
  ysup <- apply(xsup, 1, f)
  ntest <- 500
  xtest <- matrix(runif(ntest*d), ntest, d)
  ytest <- apply(xtest, 1, f)
  
  # Create with only supp
  expect_error(s1 <- SGGPcreate(d, 0, Xs=xsup, Ys=ysup, corr="CauchySQ"), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  
  # Predict with only supp on supp points
  expect_error(p1sup <- SGGPpred(s1, xsup), NA)
  expect_equal(c(p1sup$mean), ysup, tol=1e-6)
  expect_true(all(p1sup$var < 1e-8))
  
  # Predict with only supp on test points
  expect_error(p1test <- SGGPpred(s1, xtest), NA)
  expect_true(cor(c(p1test$mean), ytest) > .8)
  expect_true(all(p1test$var > 1e-8))
  expect_true(all(p1test$var < var(ytest)))
  
  # Append points, all three methods
  set.seed(0) # Fails on test, but never in console
  for (sel.method in c("Greedy", "UCB", "TS")) {
    expect_error(s1.app <- SGGPappend(s1, 100, sel.method), NA)
    expect_true(nrow(s1.app$design) > 90)
    expect_true(nrow(s1.app$design) < 100) # Can't get 100 since first block is size 1
    expect_equal(s1.app$design, s1.app$design_unevaluated)
    expect_true(all(s1.app$uo[1,] == 1)) # First block is initial
    expect_true(sum(s1.app$uo[2,]) == d+1) # 2nd block only has one 2
    # Make sure 3rd dim is least explored
    s1.app.colMeans <- colMeans(s1.app$uo[1:s1.app$uoCOUNT,])
    # print(c(sel.method,s1.app.colMeans))
    expect_true(s1.app.colMeans[1]+.1 > s1.app.colMeans[3], info = paste("s1s3",s1.app.colMeans, sel.method, collapse = " "))
    expect_true(s1.app.colMeans[2]+.1 > s1.app.colMeans[3], info = paste("s2s3",s1.app.colMeans, sel.method, collapse = " "))
  }
  set.seed(Sys.time())
  rm(s1, s1.app, s1.app.colMeans)
  
  # Again create with only supp, but use MCMC
  expect_error(s1 <- SGGPcreate(d, 0, Xs=xsup, Ys=ysup,
                                supp_args=list(laplaceapprox=FALSE,
                                               numPostSamples=7)), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  
})



test_that("2. Create, append, predict with only supp, MVout, yes PCA, yes sepOPD", {
  d <- 5
  # f <- function(x){1*(cos(x[3]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + log(x[3]+.2))}
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {cbind(apply(x, 1, f1), apply(x, 1, f2), apply(x, 1, f3))}
    else {c(f1(x), f2(x), f3(x))}
  }
  d_out <- 3
  d_outpca <- 2
  
  nsup <- 40
  xsup <- matrix(runif(nsup*d), nsup, d)
  ysup <- f(xsup)
  ntest <- 500
  xtest <- matrix(runif(ntest*d), ntest, d)
  ytest <- f(xtest)
  
  # Create with only supp
  expect_error(s1 <- SGGPcreate(d, 0, corr="Cauchy",
                                Xs=xsup, Ys=ysup, supp_args=list(use_PCA=TRUE, separateoutputparameterdimensions=TRUE)), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  # Reduced dim from PCA
  expect_equal(ncol(s1$Ys), d_out)
  expect_equal(ncol(s1$ys), d_outpca)
  # theta correct dim
  expect_equal(ncol(s1$thetaMAP), d_outpca)
  
  # Predict with only supp on supp points
  expect_error(p1sup <- SGGPpred(s1, xsup), NA)
  expect_equal(p1sup$mean, ysup, tol=1e-6)
  # expect_true(all(p1sup$var < 1e-8))
  
  # Predict with only supp on test points
  expect_error(p1test <- SGGPpred(s1, xtest), NA)
  expect_true(all(diag(cor(p1test$mean, ytest)) > .8))
  # expect_true(all(p1test$var > 1e-8))
  # expect_true(all(p1test$var < var(ytest)))
  
  # Append points, all three methods
  for (sel.method in c("Greedy", "UCB", "TS")) {
    expect_error(s1.app <- SGGPappend(s1, 100, sel.method), NA)
    expect_true(nrow(s1.app$design) > 90)
    expect_true(nrow(s1.app$design) < 100) # Can't get 100 since first block is size 1
    expect_equal(s1.app$design, s1.app$design_unevaluated)
    expect_true(all(s1.app$uo[1,] == 1)) # First block is initial
    expect_true(sum(s1.app$uo[2,]) == d+1) # 2nd block only has one 2
    # Make sure 3rd dim is least explored
    s1.app.colMeans <- colMeans(s1.app$uo[1:s1.app$uoCOUNT,])
    expect_true(s1.app.colMeans[1]+.1 > s1.app.colMeans[3])
    expect_true(s1.app.colMeans[2]+.1 > s1.app.colMeans[3])
  }
  
})

test_that("3. Create, append, predict with only supp, MVout, no PCA, yes sepOPD", {
  d <- 5
  # f <- function(x){1*(cos(x[3]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + log(x[3]+.2))}
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {cbind(apply(x, 1, f1), apply(x, 1, f2), apply(x, 1, f3))}
    else {c(f1(x), f2(x), f3(x))}
  }
  d_out <- 3
  d_outpca <- 2
  
  nsup <- 40
  xsup <- matrix(runif(nsup*d), nsup, d)
  ysup <- f(xsup)
  ntest <- 500
  xtest <- matrix(runif(ntest*d), ntest, d)
  ytest <- f(xtest)
  
  # Create with only supp
  expect_error(s1 <- SGGPcreate(d, 0, corr="Cauchy",
                                Xs=xsup, Ys=ysup, supp_args=list(use_PCA=FALSE, separateoutputparameterdimensions=TRUE)), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  # No PCA, no reduced dims
  expect_equal(ncol(s1$Ys), d_out)
  expect_equal(ncol(s1$ys), d_out)
  # theta correct dim
  expect_equal(ncol(s1$thetaMAP), d_out)
  
  # Predict with only supp on supp points
  expect_error(p1sup <- SGGPpred(s1, xsup), NA)
  expect_equal(p1sup$mean, ysup, tol=1e-6)
  # expect_true(all(p1sup$var < 1e-8))
  
  # Predict with only supp on test points
  expect_error(p1test <- SGGPpred(s1, xtest), NA)
  expect_true(all(diag(cor(p1test$mean, ytest)) > .8))
  # expect_true(all(p1test$var > 1e-8))
  # expect_true(all(p1test$var < var(ytest)))
  
  # Append points, all three methods
  for (sel.method in c("Greedy", "UCB", "TS")) {
    expect_error(s1.app <- SGGPappend(s1, 100, sel.method), NA)
    expect_true(nrow(s1.app$design) > 90)
    expect_true(nrow(s1.app$design) < 100) # Can't get 100 since first block is size 1
    expect_equal(s1.app$design, s1.app$design_unevaluated)
    expect_true(all(s1.app$uo[1,] == 1)) # First block is initial
    expect_true(sum(s1.app$uo[2,]) == d+1) # 2nd block only has one 2
    # Make sure 3rd dim is least explored
    s1.app.colMeans <- colMeans(s1.app$uo[1:s1.app$uoCOUNT,])
    expect_true(s1.app.colMeans[1]+.1 > s1.app.colMeans[3])
    expect_true(s1.app.colMeans[2]+.1 > s1.app.colMeans[3])
  }
  
})


test_that("4. Create, append, predict with only supp, MVout, yes PCA, no sepOPD", {
  d <- 5
  # f <- function(x){1*(cos(x[3]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + log(x[3]+.2))}
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {cbind(apply(x, 1, f1), apply(x, 1, f2), apply(x, 1, f3))}
    else {c(f1(x), f2(x), f3(x))}
  }
  d_out <- 3
  d_outpca <- 2
  nopd <- d_outpca
  
  nsup <- 40
  xsup <- matrix(runif(nsup*d), nsup, d)
  ysup <- f(xsup)
  ntest <- 500
  xtest <- matrix(runif(ntest*d), ntest, d)
  ytest <- f(xtest)
  
  # Create with only supp
  expect_error(s1 <- SGGPcreate(d, 0, corr="M32",
                                Xs=xsup, Ys=ysup, supp_args=list(use_PCA=TRUE, separateoutputparameterdimensions=FALSE)), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  # Using PCA, reduced dims
  expect_equal(ncol(s1$Ys), d_out)
  expect_equal(ncol(s1$ys), d_outpca)
  # theta correct dim
  expect_is(s1$thetaMAP, "numeric")
  
  # Predict with only supp on supp points
  expect_error(p1sup <- SGGPpred(s1, xsup), NA)
  expect_equal(p1sup$mean, ysup, tol=1e-6)
  # expect_true(all(p1sup$var < 1e-8))
  
  # Predict with only supp on test points
  expect_error(p1test <- SGGPpred(s1, xtest), NA)
  expect_true(all(diag(cor(p1test$mean, ytest)) > .6)) # Usually good, but can go below .8
  # expect_true(all(p1test$var > 1e-8))
  # expect_true(all(p1test$var < var(ytest)))
  
  # Append points, all three methods
  for (sel.method in c("Greedy", "UCB", "TS")) {
    expect_error(s1.app <- SGGPappend(s1, 100, sel.method), NA)
    expect_true(nrow(s1.app$design) > 90)
    expect_true(nrow(s1.app$design) < 100) # Can't get 100 since first block is size 1
    expect_equal(s1.app$design, s1.app$design_unevaluated)
    expect_true(all(s1.app$uo[1,] == 1)) # First block is initial
    expect_true(sum(s1.app$uo[2,]) == d+1) # 2nd block only has one 2
    # Make sure 3rd dim is least explored
    s1.app.colMeans <- colMeans(s1.app$uo[1:s1.app$uoCOUNT,])
    expect_true(s1.app.colMeans[1]+.1 > s1.app.colMeans[3])
    expect_true(s1.app.colMeans[2]+.1 > s1.app.colMeans[3])
  }
  
})



test_that("5. Create, append, predict with only supp, MVout, no PCA, no sepOPD", {
  d <- 5
  # f <- function(x){1*(cos(x[3]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + log(x[3]+.2))}
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {cbind(apply(x, 1, f1), apply(x, 1, f2), apply(x, 1, f3))}
    else {c(f1(x), f2(x), f3(x))}
  }
  d_out <- 3
  d_outpca <- 2
  nopd <- d_outpca
  
  nsup <- 40
  xsup <- matrix(runif(nsup*d), nsup, d)
  ysup <- f(xsup)
  ntest <- 500
  xtest <- matrix(runif(ntest*d), ntest, d)
  ytest <- f(xtest)
  
  # Create with only supp
  expect_error(s1 <- SGGPcreate(d, 0, corr="PowerExp",
                                Xs=xsup, Ys=ysup, supp_args=list(use_PCA=FALSE, separateoutputparameterdimensions=FALSE)), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  # No PCA, no reduced dims
  expect_equal(ncol(s1$Ys), d_out)
  expect_equal(ncol(s1$ys), d_out)
  # theta correct dim
  expect_is(s1$thetaMAP, "numeric")
  
  # Predict with only supp on supp points
  expect_error(p1sup <- SGGPpred(s1, xsup), NA)
  expect_equal(p1sup$mean, ysup, tol=1e-6)
  # expect_true(all(p1sup$var < 1e-8))
  
  # Predict with only supp on test points
  expect_error(p1test <- SGGPpred(s1, xtest), NA)
  expect_true(all(diag(cor(p1test$mean, ytest)) > .6)) # corr can go below .8
  # expect_true(all(p1test$var > 1e-8))
  # expect_true(all(p1test$var < var(ytest)))
  
  # Append points, all three methods
  for (sel.method in c("Greedy", "UCB", "TS")) {
    expect_error(s1.app <- SGGPappend(s1, 100, sel.method), NA)
    expect_true(nrow(s1.app$design) > 90)
    expect_true(nrow(s1.app$design) < 100) # Can't get 100 since first block is size 1
    expect_equal(s1.app$design, s1.app$design_unevaluated)
    expect_true(all(s1.app$uo[1,] == 1)) # First block is initial
    expect_true(sum(s1.app$uo[2,]) == d+1) # 2nd block only has one 2
    # Make sure 3rd dim is least explored
    s1.app.colMeans <- colMeans(s1.app$uo[1:s1.app$uoCOUNT,])
    expect_true(s1.app.colMeans[1]+.3 > s1.app.colMeans[3]) # Hard to get these 100%
    expect_true(s1.app.colMeans[2]+.3 > s1.app.colMeans[3])
  }
  
})
