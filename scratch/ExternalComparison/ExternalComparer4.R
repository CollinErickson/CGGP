# Comparison for paper

# decentLHS <- sFFLHD::decentLHS
decentLHS <- function(n, d, ndes, max.time) {
  if (missing(ndes) && missing(max.time)) {
    stop("Must give ndes or max.time to decentLHS")
  }
  start.time <- Sys.time()
  bestdes <- NULL
  bestcrit <- Inf
  i <- 1
  while (TRUE) {
    # Make new LHS
    x <- lhs::randomLHS(n, d)
    crit <- MaxPro::MaxProMeasure(x)
    
    # Check if best yet
    if (crit < bestcrit) {
      bestcrit <- crit
      bestdes <- x
    }
    
    # Check if done
    i <- i+1
    if (!missing(ndes) && (i > ndes)) {
      break
    }
    if (!missing(max.time) && (as.numeric(Sys.time() - start.time, units="secs") > max.time)) {
      break
    }
  }
  # Return bestdes
  bestdes
}

run_lagp <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed, use_agp=FALSE) {
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  sdy <- sd(y)
  mny <- mean(y)
  y <- (y-mny) / sdy
  if (use_agp) {
    stop("Don't use this version")
    fit.time.start <- Sys.time()
    fit.time.end <- Sys.time()
    pred.time.start <- Sys.time()
    pred <- laGP::aGPsep(X=x, Z=y, XX=xtest, method="alc")
    pred$var <- pred$var * sdy^2
  } else {
    fit.time.start <- Sys.time()
    mod.agp <- laGP::newGPsep(X=x, Z=y, d=laGP::darg(d=list(mle = TRUE, max = 100), X=x)$start,
                              g=laGP::garg(g=list(mle = TRUE), y=y)$start)
    laGP::updateGPsep(mod.agp, x, y)
    fit.time.end <- Sys.time()
    
    pred.time.start <- Sys.time()
    pred <- laGP::predGPsep(mod.agp, xtest, lite=T)
    pred$var <- pred$s2 * sdy^2
  }
  # browser()
  pred$mean <- pred$mean * sdy + mny
  pred.time.end <- Sys.time()
  list(mean=pred$mean, var=pred$var, n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}

run_MRFA <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  fit.time.start <- Sys.time()
  mod <- MRFA::MRFA_fit(X=x, Y=y, verbose=FALSE)
  fit.time.end <- Sys.time()
  pred.time.start <- Sys.time()
  pred <- predict(mod, xtest, lambda = min(mod$lambda))
  pred.time.end <- Sys.time()
  list(mean=pred$y_hat, var=rep(NaN, nrow(xtest)), n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}


run_svm <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {#browser()
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  require(e1071)
  fit.time.start <- Sys.time()
  mod <- svm(x, y)
  fit.time.end <- Sys.time()
  pred.time.start <- Sys.time()
  pred <- predict(mod, xtest)
  pred.time.end <- Sys.time()
  list(mean=pred, var=rep(NaN, nrow(xtest)), n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}

run_mlegp <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {#browser()
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  require(mlegp)
  fit.time.start <- Sys.time()
  mod <- mlegp(x, y)
  fit.time.end <- Sys.time()
  pred.time.start <- Sys.time()
  pred <- predict(mod, xtest, se.fit)
  pred.time.end <- Sys.time()
  list(mean=pred$fit, var=pred$se.fit^2, n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}



run_CGGPsupp <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed) {
  require("CGGP")
  if (!missing(seed)) {set.seed(seed)}
  xsup <- lhs::maximinLHS(Nlhs, d)
  ysup <- apply(xsup, 1, f)
  fit.time.start <- Sys.time()
  sg <- CGGPcreate(d=d, batchsize=0, Xs=xsup, Ys=ysup)
  notdone <- TRUE
  while (notdone) {
    print(sg)
    Nalready <- (if (!is.null(sg$design)) {nrow(sg$design)} else {0}) + nrow(sg$Xs)
    ni <- if (Nalready<1000) {200} else if (Nalready<10000) {500} else {2000}
    
    if (ni +Nalready > Ntotal) {
      ni <- Ntotal - Nalready
      notdone <- FALSE
    }
    sg <- CGGPappend(sg, batchsize = ni)
    
    ynew <- apply(sg$design_unevaluated, 1, f)
    sg <- CGGPfit(sg, Ynew=ynew, Xs=xsup, Ys=ysup)
  }
  fit.time.end <- Sys.time()
  print(sg)
  Nevaluated <- nrow(sg$design) + nrow(sg$Xs)
  if (Nevaluated > Ntotal) {stop("CGGP has more points than Ntotal")}
  
  pred.time.start <- Sys.time()
  pred <- predict(sg, xtest)
  pred.time.end <- Sys.time()
  list(mean=pred$mean, var=pred$var,
       n=Nevaluated,
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}




# Need a generic function that passes to specific ones
run_one <- function(package, f, d, npd, replicate) {
  n <- npd * d
  # if (n!= 500) {stop('bad n')}
  f <- eval(parse(text=paste0("TestFunctions::", f)))
  
  ntest <- 1e4
  # xtest <- matrix(runif(ntest*d), ncol=d)
  xtest <- decentLHS(ntest, d, max.time = 10)
  ytest <- apply(xtest, 1, f)
  
  # if (package == "CGGP") {
  #   out <- run_CGGP(Nappend=floor(n * (1:5)/5), f=f, d=d, xtest=xtest)
  if (package == "CGGPsupp") {
    # out <- run_CGGP(Nappend=floor(n*(2:5)/5), Nlhs=floor(.2*n), f=f, d=d, xtest=xtest)
    out <- run_CGGPsupp(Ntotal=n, Nlhs=10*d, f=f, d=d, xtest=xtest)
    # } else if (package == "CGGPsupponly") {
    #   out <- run_CGGP(Nappend=c(), Nlhs=n, f=f, d=d, xtest=xtest)
  } else if (package == "laGP") {
    out <- run_lagp(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "aGP") {
    out <- run_lagp(Ntotal=n, f=f, d=d, xtest=xtest, use_agp=TRUE)
  } else if (package == "MRFA") {
    out <- run_MRFA(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "svm") {
    out <- run_svm(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "mlegp") {
    out <- run_mlegp(Ntotal=n, f=f, d=d, xtest=xtest)
  } else {
    stop(paste("Package", package, "not recognized"))
  }
  # browser()
  outstats <- CGGP::valstats(predmean=out[[1]], predvar=out[[2]],Yval=ytest) #, bydim=FALSE)
  if (out$n > n) {warning(paste("n too big for", package, n, f, d))}
  outstats
}




require("comparer")

excomp <- ffexp$new(
  eval_func = run_one, #run_one_parallel,
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"),
                d=c(3,6,7,8,10),
                row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
  package=c("CGGPsupp", "MRFA", "svm", "aGP"),
  npd=c(10, 30, 100, 300, 1000, 3000, 10000),
  parallel=T,
  parallel_cores = 10,
  replicate=1:5,
  # folder_path= "/home/collin/scratch/CGGP/scratch/ExternalComparison/ExComp4"
  folder_path="./scratch/ExternalComparison/ExComp4/"
)
excomp$run_all(save_output = F, parallel = F, run_order = "random")

excomp$rungrid
# try because it gave delete error before, but shouldn't need it now
try(excomp$recover_parallel_temp_save(delete_after = FALSE))
excomp$save_self()
excomp$run_all(parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE,
               write_start_files=T, write_error_files=T)
# excomp$run_all()

cat("Completed all runs in ExternalComparer4.R\n")

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
