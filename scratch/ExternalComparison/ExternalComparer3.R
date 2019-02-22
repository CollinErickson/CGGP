# Compare ours with and without supp

run_lagp <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed, use_agp=FALSE) {
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    # x <- matrix(runif(Ntotal*d), Ntotal, d)
    x <- lhs::maximinLHS(Ntotal, d)
    y <- apply(x, 1, f)
  }
  sdy <- sd(y)
  mny <- mean(y)
  y <- (y-mny) / sdy
  if (use_agp) {
    pred <- laGP::aGPsep(X=x, Z=y, XX=xtest, method="alc")
    pred$var <- pred$var * sdy^2
  } else {
    mod.agp <- laGP::newGPsep(X=x, Z=y, d=laGP::darg(d=list(mle = TRUE, max = 100), X=x)$start,
                              g=laGP::garg(g=list(mle = TRUE), y=y)$start)
    laGP::updateGPsep(mod.agp, x, y)
    pred <- laGP::predGPsep(mod.agp, xtest, lite=T)
    pred$var <- pred$s2 * sdy^2
  }
  # browser()
  pred$mean <- pred$mean * sdy + mny
  list(mean=pred$mean, var=pred$var, n=nrow(x))
}

run_MRFA <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    x <- lhs::maximinLHS(Ntotal, d)
    y <- apply(x, 1, f)
  }
  mod <- MRFA::MRFA_fit(X=x, Y=y, verbose=FALSE)
  pred <- predict(mod, xtest, lambda = min(mod$lambda))
  list(mean=pred$y_hat, var=rep(NaN, nrow(xtest)), n=nrow(x))
}


run_svm <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {#browser()
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (missing(x) && missing(y)) {
    x <- lhs::maximinLHS(Ntotal, d)
    y <- apply(x, 1, f)
  }
  require(e1071)
  mod <- svm(x, y)
  pred <- predict(mod, xtest)
  list(mean=pred, var=rep(NaN, nrow(xtest)), n=nrow(x))
}



run_SGGP <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed) {
  require("SGGP")
  if (!missing(seed)) {set.seed(seed)}
  # browser()
  if (!missing(Nlhs)) {
    # xsup <- matrix(runif(Nlhs*d), Nlhs, d)
    xsup <- lhs::maximinLHS(Nlhs, d)
    ysup <- apply(xsup, 1, f)
    sg <- SGGPcreate(d=d, batchsize=0, Xs=xsup, Ys=ysup)
  } else {
    xsup <- NULL
    ysup <- NULL
  }
  if (!is.null(Nappend) && length(Nappend) > 0) {
    if (any(abs(floor(Nappend)-Nappend)>1e-8)) {stop("Need integers")}
    for (i in 1:length(Nappend)) {
      ni <- Nappend[i]
      if (exists('sg')) {
        if (!is.null(sg$design)) {ni <- ni - nrow(sg$design)}
        if (!is.null(sg$Xs)) {ni <- ni - nrow(sg$Xs)}
      }
      if (missing(Nlhs) && i==1) { # First time create it
        sg <- SGGPcreate(d, ni)
      } else{
        sg <- SGGPappend(sg, batchsize = ni)
      }
      ynew <- apply(sg$design_unevaluated, 1, f)
      sg <- SGGPfit(sg, Ynew=ynew, Xs=xsup, Ys=ysup)
    }
  }
  
  pred <- predict(sg, xtest)
  list(mean=pred$mean, var=pred$var,
       n=if (!is.null(sg$design)) {nrow(sg$design)} else {0} + if (!is.null(sg$Xs)) {nrow(sg$Xs)} else {0})
}




# Need a generic function that passes to specific ones
run_one <- function(package, f, d, n) {
  # if (n!= 500) {stop('bad n')}
  f <- eval(parse(text=paste0("TestFunctions::", f)))
  
  ntest <- 1e3
  xtest <- matrix(runif(ntest*d), ncol=d)
  ytest <- apply(xtest, 1, f)
  
  if (package == "SGGP") {
    out <- run_SGGP(Nappend=floor(n * (1:5)/5), f=f, d=d, xtest=xtest)
  } else if (package == "SGGPsupp") {
    out <- run_SGGP(Nappend=floor(n*(2:5)/5), Nlhs=floor(.2*n), f=f, d=d, xtest=xtest)
  } else if (package == "SGGPsupponly") {
    out <- run_SGGP(Nappend=c(), Nlhs=n, f=f, d=d, xtest=xtest)
  } else if (package == "laGP") {
    out <- run_lagp(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "aGP") {
    out <- run_lagp(Ntotal=n, f=f, d=d, xtest=xtest, use_agp=TRUE)
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
  outstats <- SGGP::valstats(predmean=out[[1]], predvar=out[[2]],Yval=ytest) #, bydim=FALSE)
  if (out$n > n) {warning(paste("n too big for", package, n, f, d))}
  outstats
}




require("comparer")

packxn <- rbind(
  data.frame(package="SGGP", n=c(100,200,400,800,1600,3200,6400,12800), stringsAsFactors=F),
  data.frame(package="SGGPsupp", n=c(100,200,400,800,1600,3200,6400)),
  data.frame(package="SGGPsupponly", n=c(100,200,400,800)),
  data.frame(package="MRFA", n=c(100,200,400,800,1600,3200,6400,12800)),
  data.frame(package="svm", n=c(100,200,400,800,1600,3200,6400,12800)),
  data.frame(package="laGP", n=c(100,200,400,800)),
  data.frame(package="aGP", n=c(100,200,400,800,1600,3200,6400,12800))
)

funcstouse <- 1:5
excomp <- ffexp$new(
  eval_func = run_one_parallel,
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight")[funcstouse],
                d=c(3,6,7,8,10)[funcstouse],
                row.names = c("beam","OTL","piston","borehole","wingweight")[funcstouse], stringsAsFactors = F),
  # package=c("SGGP", "SGGPsupp", "SGGPsupponly", "MRFA", "svm", "laGP", "aGP"), # GPflow.SVGP is slow/bad
  # n = c(100, 200), #, 250, 500, 750, 1000),
  packxn=packxn,
  parallel=T,
  parallel_cores = 4,
  # replicate=1:10,
  # folder_path= "/home/collin/scratch/SGGP/scratch/ExternalComparison/ComparerRun6" #"./scratch/sggpout"
  folder_path="./scratch/ExternalComparison/ExComp3/"
)

excomp$rungrid
# try because it gave delete error before, but shouldn't need it now
try(excomp$recover_parallel_temp_save(delete_after = FALSE))
excomp$save_self()
excomp$run_all(parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE)
# excomp$run_all()

cat("Completed all runs in ExternalComparer3.R\n")

# excomp$save_self()

cat("Saved self\n")

if (F) {
  excomp$plot_run_times()
  plyr::dlply(excomp$outcleandf, "d")
  require('ggplot2')
  ecdf <- excomp$outcleandf[excomp$completed_runs,]
  ggplot(data=ecdf, mapping=aes(n, RMSE, color=as.factor(n))) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  ggplot(data=ecdf, mapping=aes(n, score, color=as.factor(n))) + geom_point() + facet_grid(f ~ package, scales="free_y")
  ggplot(data=ecdf, mapping=aes(n, CRPscore)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  ggplot(data=ecdf, mapping=aes(n, runtime)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  # saveRDS(excomp, "./scratch/ExternalComparison/ExComp1_completed.rds")
}
