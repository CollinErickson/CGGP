# files to run tests against other software

# for each test, will define: Ntotal, Nappendsizes (only we can use this now),
#  true function, num of input dim.
#  Num of output dim? Not yet
# They should return predictions as list with mean and var
run_lagp <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    x <- matrix(runif(Ntotal*d), Ntotal, d)
    y <- apply(x, 1, f)
  }
  mod.agp <- laGP::newGPsep(X=x, Z=y, d=laGP::darg(d=list(mle = TRUE, max = 100), X=x)$start,
                            g=laGP::garg(g=list(mle = TRUE), y=y)$start)
  laGP::updateGPsep(mod.agp, x, y)
  pred <- laGP::predGPsep(mod.agp, xtest, lite=T)
  list(mean=pred$mean, var=pred$s2)
}

run_MRFA <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    x <- matrix(runif(Ntotal*d), Ntotal, d)
    y <- apply(x, 1, f)
  }
  mod <- MRFA::MRFA_fit(X=x, Y=y)
  pred <- predict(mod, xtest, lambda = min(mod$lambda))
  list(mean=pred$y_hat)
}


run_SGGP <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  browser()
  # if (missing(x) && missing(y)) {
  #   x <- matrix(runif(Ntotal*d), Ntotal, d)
  #   y <- apply(x, 1, f)
  # }
  for (i in 1:length(Nappend)) {
    ni <- Nappend[i]
    if (i==1) {sg <- SGGPcreate(d=d, batchsize=ni)}
    else {sg <- SGGPappend(sg, batchsize = ni)}
    ynew <- apply(sg$design_unevaluated, 1, f)
    sg <- SGGPfit(sg, Ynew=ynew)
  }
  
  pred <- predict(sg, xtest)
  pred
}

f <- TestFunctions::borehole
d <- 8
ntest <- 1e3
xtest <- matrix(runif(ntest*d), ntest, d)
ytest <- apply(xtest, 1, f)
lagpout <- run_lagp(Ntotal=1e3, f = f, d=d, xtest=xtest)
lagp_rmse <- sqrt(mean((lagpout$mean - ytest)^2))
mrfaout <- run_MRFA(Ntotal=1e3, f = f, d=d, xtest=xtest)
mrfa_rmse <- sqrt(mean((mrfaout$mean - ytest)^2))
sggpout <- run_SGGP(Nappend = rep(200, 5), f = f, d=d, xtest=xtest)
sggp_rmse <- sqrt(mean((sggpout$mean - ytest)^2))

data.frame(lagp_rmse, mrfa_rmse, sggp_rmse)
