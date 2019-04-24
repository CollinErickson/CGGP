# Run ExternalComparer5 on CGGP and CGGPsupp using ICC-TS


decentLHS <- function(n, d, ndes, max.time) {
  if (missing(ndes) && missing(max.time)) {
    stop("Must give ndes or max.time to decentLHS")
  }
  
  # MaxProMeasure in 10d with 3e4 pts takes 3 min, crashes with 1e5 pts.
  if (n > 1e4) {
    return(lhs::randomLHS(n, d))
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




run_CGGP <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed, selection.method, correlation) {#browser()
  require("CGGP")
  if (!missing(seed)) {set.seed(seed)}
  if (!missing(Nlhs) && Nlhs!=0) {stop("Nlhs given to run_CGGP")}
  fit.time.start <- Sys.time()
  sg <- CGGPcreate(d=d, batchsize=min(Ntotal,200), corr=correlation)
  sg <- CGGPfit(sg, apply(sg$design, 1, f))
  notdone <- (nrow(sg$design) < Ntotal)
  while (notdone) {
    # print(sg)
    Nalready <- nrow(sg$design)
    ni <- if (Nalready<1000) {200} else if (Nalready<10000) {500} else if (Nalready<20000) {2000} else {10000}
    
    if (ni +Nalready > Ntotal) {
      ni <- Ntotal - Nalready
      notdone <- FALSE
    }
    sg <- CGGPappend(sg, batchsize = ni, selectionmethod = selection.method)
    
    if (!is.null(sg$design_unevaluated)) {
      ynew <- apply(sg$design_unevaluated, 1, f)
      sg <- CGGPfit(sg, Ynew=ynew)
    } else {
      # print('Nothing new to evaluate, hopefully !notdone')
    }
  }
  fit.time.end <- Sys.time()
  # print(sg)
  Nevaluated <- if (!is.null(sg$design)) nrow(sg$design) else {0} + nrow(sg$Xs)
  if (Nevaluated > Ntotal) {stop("CGGP has more points than Ntotal")}
  
  pred.time.start <- Sys.time()
  pred <- predict(sg, xtest)
  pred.time.end <- Sys.time()
  list(mean=pred$mean, var=pred$var,
       n=Nevaluated,
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}

run_CGGPsupp <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed, HandlingSuppData, selection.method, correlation) {#browser()
  require("CGGP")
  if (!missing(seed)) {set.seed(seed)}
  xsup <- lhs::maximinLHS(Nlhs, d)
  ysup <- apply(xsup, 1, f)
  fit.time.start <- Sys.time()
  sg <- CGGPcreate(d=d, batchsize=0, Xs=xsup, Ys=ysup, corr=correlation)
  notdone <- (Nlhs < Ntotal)
  while (notdone) {
    # print(sg)
    Nalready <- (if (!is.null(sg$design)) {nrow(sg$design)} else {0}) + nrow(sg$Xs)
    ni <- if (Nalready<1000) {200} else if (Nalready<10000) {500} else if (Nalready<20000) {2000} else {10000}
    
    if (ni +Nalready > Ntotal) {
      ni <- Ntotal - Nalready
      notdone <- FALSE
    }
    sg <- CGGPappend(sg, batchsize = ni, selection.method)
    
    if (!is.null(sg$design_unevaluated)) {
      ynew <- apply(sg$design_unevaluated, 1, f)
      sg <- CGGPfit(sg, Ynew=ynew, Xs=xsup, Ys=ysup, HandlingSuppData = HandlingSuppData)
    } else {
      # print('Nothing new to evaluate, hopefully !notdone')
    }
  }
  fit.time.end <- Sys.time()
  # print(sg)
  Nevaluated <- if (!is.null(sg$design)) nrow(sg$design) else {0} + nrow(sg$Xs)
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
run_one <- function(package, selection.method, correlation, HandlingSuppData,
                    f, d, npd, replicate) {#browser()
  # package <- psch$package
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
    out <- run_CGGPsupp(Ntotal=n, Nlhs=10*d, f=f, d=d, xtest=xtest, HandlingSuppData=HandlingSuppData, selection.method=selection.method, correlation=correlation)
    # } else if (package == "CGGPsupponly") {
    #   out <- run_CGGP(Nappend=c(), Nlhs=n, f=f, d=d, xtest=xtest)
  } else if (package == "CGGPsupponly") {
    out <- run_CGGPsupponly(Ntotal=n, f=f, d=d, xtest=xtest, correlation=correlation)
  } else if (package == "CGGPoneshot") {
    out <- run_CGGPoneshot(Ntotal=n, f=f, d=d, xtest=xtest, correlation=correlation)
  } else if (package == "CGGP") {
    out <- run_CGGP(Ntotal=n, f=f, d=d, xtest=xtest, selection.method=selection.method, correlation=correlation)
  } else if (package == "laGP") {
    out <- run_lagp(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "aGP") {
    out <- run_lagp_bobby(Ntotal=n, f=f, d=d, xtest=xtest, use_agp=TRUE)
  } else if (package == "aGP2") {
    out <- run_lagp_matt(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "MRFA") {
    out <- run_MRFA(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "svm") {
    out <- run_svm(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "mlegp") {
    out <- run_mlegp(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "GPfit") {
    out <- run_GPfit(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "BASS") {
    out <- run_BASS(Ntotal=n, f=f, d=d, xtest=xtest)
  } else {
    stop(paste("Package", package, "not recognized"))
  }
  # browser()
  outstats <- CGGP::valstats(predmean=out[[1]], predvar=out[[2]],Yval=ytest) #, bydim=FALSE)
  if (out$n > n) {warning(paste("n too big for", package, n, f, d))}
  outstats$predtime <- out$pred.time
  outstats$fittime  <- out$fit.time
  # Forgot to add n_used, outstats$n_used <- out$n
  outstats
}







expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

eg1TS <- expand.grid(selection.method = c("TS"),
                   correlation = c("CauchySQ", "Matern32", "PowerExp", "Cauchy", "CauchySQT"), stringsAsFactors=F)
eg2aTS <- expand.grid.df(eg1TS, data.frame(HandlingSuppData=c("Correct"), stringsAsFactors=F), data.frame(package="CGGPsupp", stringsAsFactors=F))
eg2bTS <- expand.grid.df(eg1TS, data.frame(HandlingSuppData="NA", stringsAsFactors=F), data.frame(package=c("CGGP"), stringsAsFactors=F))
# eg2c <- expand.grid(selection.method="NA", correlation = c("CauchySQ", "Matern32", "PowerExp", "Cauchy", "CauchySQT"), HandlingSuppData="NA", package=c("CGGPoneshot", "CGGPsupponly"), stringsAsFactors=F)
# eg2d <- data.frame(selection.method="NA", correlation="NA", HandlingSuppData="NA",
#                    package=c("MRFA", "svm", "aGP", "aGP2", "laGP", "mlegp", "GPfit", "BASS"), stringsAsFactors=F)
eg3TS <- rbind(eg2aTS, eg2bTS)

require("comparer")

excompTS <- ffexp$new(
  eval_func = run_one,
  varlist = c("decentLHS", "run_CGGP", "run_CGGPsupp"),
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"),
                d=c(3,6,7,8,10),
                row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
  psch=eg3TS,
  npd=c(10, 30, 100, 300, 1000, 3000, 10000),
  parallel=if (version$os =="linux-gnu") {TRUE} else {!FALSE},
  parallel_cores = if (version$os =="linux-gnu") {35} else {3},
  replicate=1:10, #:5,
  folder_path= if (version$os =="linux-gnu") {"/home/collin/scratch/SGGP/scratch/ExternalComparison/ExComp5TS/"}
  else {"./scratch/ExternalComparison/ExComp5TS/"}
)
# Remove ones that can't do full size
package.name <- excompTS$arglist$psch$package[excompTS$rungrid$psch]
npd.excompTS <- excompTS$arglist$npd[excompTS$rungrid$npd]
n.excompTS <- npd.excompTS * excompTS$arglist$fd$d[excompTS$rungrid$fd]


table(excompTS$completed_runs)
try(excompTS$recover_parallel_temp_save(delete_after = FALSE))
table(excompTS$completed_runs)
excompTS$save_self()

excompTS$run_all(parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE,
               write_start_files=T, write_error_files=T)


cat("Completed all runs in ExternalComparer5TS.R\n")

excompTS$save_self()

cat("Saved self\n")
