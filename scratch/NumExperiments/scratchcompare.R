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


run_GPflow <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    x <- matrix(runif(Ntotal*d), Ntotal, d)
    y <- apply(x, 1, f)
  }
  require("reticulate")
  env1 <- use_condaenv("myenv")
  py$x <- x
  py$y <- as.matrix(y)
  py$xtest <- xtest
  py_run_string("import gpflow")
  py_run_string("import numpy as np")
  py_run_string("import matplotlib")
  
  py_run_string(paste0("k = gpflow.kernels.Matern52( ", d,", lengthscales=0.3)"))
  py_run_string("m = gpflow.models.GPR(x, y, kern=k)")
  py_run_string("m.likelihood.variance = 0.0001")
  py_run_string("gpflow.train.ScipyOptimizer().minimize(m)")
  py_run_string("r.gpflowmean, r.gpflowvar = m.predict_y(xtest)")
  # mod <- MRFA::MRFA_fit(X=x, Y=y)
  # pred <- predict(mod, xtest, lambda = min(mod$lambda))
  list(mean=gpflowmean, var=gpflowvar)
}


run_SGGP <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  # browser()
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
lagp_stats <- valstats(lagpout$mean, lagpout$var, ytest)
mrfaout <- run_MRFA(Ntotal=1e3, f = f, d=d, xtest=xtest)
mrfa_rmse <- sqrt(mean((mrfaout$mean - ytest)^2))
mrfa_stats <- valstats(mrfaout$mean, rep(NA, ntest), ytest)
sggpout <- run_SGGP(Nappend = rep(200, 5), f = f, d=d, xtest=xtest)
sggp_rmse <- sqrt(mean((sggpout$mean - ytest)^2))
sggp_stats <- valstats(sggpout$mean, sggpout$var, ytest)
gpflowout <- run_GPflow(Ntotal=1e3, f = f, d=d, xtest=xtest)
gpflow_rmse <- sqrt(mean((gpflowout$mean - ytest)^2))
gpflow_stats <- valstats(gpflowout$mean, gpflowout$var, ytest)

data.frame(lagp_rmse, mrfa_rmse, sggp_rmse, gpflow_rmse)
rbind(laGP=lagp_stats, MRFA=mrfa_stats, SGGP=sggp_stats, GPflow=gpflow_stats)
