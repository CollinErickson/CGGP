# This will run test against other software, just one output dimension for now.

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
  sdy <- sd(y)
  mny <- mean(y)
  y <- (y-mny) / sdy
  # mod.agp <- laGP::newGPsep(X=x, Z=y, d=laGP::darg(d=list(mle = TRUE, max = 100), X=x)$start,
  #                           g=laGP::garg(g=list(mle = TRUE), y=y)$start)
  # laGP::updateGPsep(mod.agp, x, y)
  pred <- laGP::aGPsep(X=x, Z=y, XX=xtest, method="alc")
  # browser()
  pred$mean <- pred$mean * sdy + mny
  pred$var <- pred$var * sdy^2
  list(mean=pred$mean, var=pred$var)
}

run_MRFA <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    x <- matrix(runif(Ntotal*d), Ntotal, d)
    y <- apply(x, 1, f)
  }
  mod <- MRFA::MRFA_fit(X=x, Y=y, verbose=FALSE)
  pred <- predict(mod, xtest, lambda = min(mod$lambda))
  list(mean=pred$y_hat, var=rep(NaN, nrow(xtest)))
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


run_GPflow.SVGP <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  # Model type can be GPR, VGP, or SVGP
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    x <- matrix(runif(Ntotal*d), Ntotal, d)
    y <- f(x) #apply(x, 1, f)
  }
  require("reticulate")
  env1 <- use_condaenv("myenv")
  py$x <- x
  colMeansY <- colMeans(as.matrix(y))
  py$y <- as.matrix(y) - colMeansY
  py$xtest <- xtest
  Sys.sleep(3) # Try sleep to avoid it not finding x, maybe it just needs more time?
  py_run_string("import gpflow")
  py_run_string("import numpy as np")
  # browser()
  py_run_string(paste0("k = gpflow.kernels.Matern52( ", d,", lengthscales=0.3)"))
  # cat("model type is ", model.type, '\n')
  py$z = x[sample(1:nrow(x), min(nrow(x), 100), replace=FALSE),]
  py_run_string(paste0("m = gpflow.models.SVGP(x, y, Z=z, kern=k, likelihood=gpflow.likelihoods.Gaussian(), minibatch_size=100)"))
  py_run_string(paste0("m.X.set_batch_size(100)"))
  py_run_string(paste0("m.Y.set_batch_size(100)"))
  # py_run_string("m.likelihood.variance = 0.0001")
  # py_run_string("gpflow.train.ScipyOptimizer().minimize(m)")
  
  # for (iii in 1:20) {py_run_string("gpflow.train.ScipyOptimizer().minimize(m)");py_run_string("r.gpflowmean, r.gpflowvar = m.predict_y(xtest)");print(colMeans(gpflowmean))}
  start.time <- Sys.time()
  # for (iii in 1:3) {
  while(as.numeric(Sys.time() - start.time, unit="secs") < 180) {
    py_run_string("gpflow.train.AdamOptimizer().minimize(m)")
    py_run_string("r.gpflowmean, r.gpflowvar = m.predict_y(xtest)")
    cat(colMeans(as.matrix(y)), colMeans(gpflowmean), as.numeric(Sys.time() - start.time, unit="secs"), "\n")
  }
  
  py_run_string("r.gpflowmean, r.gpflowvar = m.predict_y(xtest)")
  list(mean=gpflowmean + colMeansY, var=gpflowvar)
}

run_svm <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {#browser()
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    x <- matrix(runif(Ntotal*d), Ntotal, d)
    y <- apply(x, 1, f)
  }
  require(e1071)
  mod <- svm(x, y)
  pred <- predict(mod, xtest)
  list(mean=pred, var=rep(NaN, nrow(xtest)))
}


run_SGGP <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  require("SGGP")
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  # if (missing(x) && missing(y)) {
  #   x <- matrix(runif(Ntotal*d), Ntotal, d)
  #   y <- apply(x, 1, f)
  # }
  if (any(abs(floor(Nappend)-Nappend)>1e-8)) {stop("Need integers")}
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



# Need a generic function that passes to specific ones
run_one <- function(package, f, d, n) {
  # if (n!= 500) {stop('bad n')}
  f <- eval(parse(text=paste0("TestFunctions::", f)))
  
  ntest <- 1e3
  xtest <- matrix(runif(ntest*d), ncol=d)
  ytest <- apply(xtest, 1, f)
  
  if (package == "SGGP") {
    # out <- run_SGGP(Nappend=c(100,200,300,400,500), f=f, d=d, xtest=xtest)
    out <- run_SGGP(Nappend=floor(n * (1:5)/5), f=f, d=d, xtest=xtest)
  } else if (package == "laGP") {
    out <- run_lagp(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "GPflow") {
    out <- run_GPflow(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "GPflow.SVGP") {
    out <- run_GPflow.SVGP(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "MRFA") {
    out <- run_MRFA(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "svm") {
    out <- run_svm(Ntotal=n, f=f, d=d, xtest=xtest)
  } else {
    stop(paste("Package", package, "not recognized"))
  }
  # browser()
  out <- SGGP::valstats(predmean=out[[1]], predvar=out[[2]],Yval=ytest) #, bydim=FALSE)
  out
}




require("comparer")

excomp <- ffexp$new(
  eval_func = run_one,
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight")[2:5], d=c(3,6,7,8,10)[2:5],
                row.names = c("beam","OTL","piston","borehole","wingweight")[2:5], stringsAsFactors = F),
  package=c("MRFA", "svm", "laGP", "SGGP"), # GPflow.SVGP is slow/bad
  # n_increments=list(n_increments=list(c(100,200,300,400,500))),
  n = c(10000), #, 250, 500, 750, 1000),
  parallel=FALSE,
  parallel_cores = 37,
  # replicate=1:10,
  folder_path= "/home/collin/scratch/SGGP/scratch/ExternalComparison/ComparerRun5" #"./scratch/sggpout"
)

excomp$rungrid
# try because it gave delete error before, but shouldn't need it now
# try(excomp$recover_parallel_temp_save(delete_after = FALSE))
# excomp$save_self()
# excomp$run_all(parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE)
excomp$run_all()

cat("Completed all runs in ExternalComparer2.R\n")

# excomp$save_self()

cat("Saved self\n")

if (F) {
  excomp$plot_run_times()
  plyr::dlply(excomp$outcleandf, "d")
  ggplot(data=excomp$outcleandf, mapping=aes(n, RMSE)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  ggplot(data=excomp$outcleandf, mapping=aes(n, score)) + geom_point() + facet_grid(f ~ package, scales="free_y")
  ggplot(data=excomp$outcleandf, mapping=aes(n, CRPscore)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  ggplot(data=excomp$outcleandf, mapping=aes(n, runtime)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  # saveRDS(excomp, "./scratch/ExternalComparison/ExComp1_completed.rds")
}
