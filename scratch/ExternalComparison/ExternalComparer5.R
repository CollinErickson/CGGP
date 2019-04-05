# Comparison for paper
# 5 is with updated grid size, extended in each dimension
#   also adding BASS into this file, last time had to do it separately

# decentLHS <- sFFLHD::decentLHS
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

run_lagp_bobby <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed, use_agp=FALSE) {
  require(laGP)
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  
  fit.time.start <- Sys.time()
  ## fixing a tiny nugget is very helpful on this problem
  g <- 1/10000000
  ## macro-scale analysis on a random subset of the data
  # browser()
  n <- min(Ntotal, 1000)
  d2 <- darg(list(mle = TRUE, max = 100), x)
  subs <- sample(1:Ntotal, n, replace = FALSE)
  gpsepi <- newGPsep(x[subs, ], y[subs], rep(d2$start, d), g=g, dK=TRUE)
  that <- mleGPsep(gpsepi, param="d", tmin=d2$min, tmax=d2$max, ab=d2$ab, maxit=200)
  # p <- predGPsep(gpsepi, xpred, lite=TRUE)
  # rmse.sub[r] <- sqrt(mean((p$mean - ypred.0)^2))
  deleteGPsep(gpsepi)
  
  ## scale the inputs according to the macro-analysis lengthscales
  scale <- sqrt(that$d)
  xs <- x
  xtests <- xtest
  for(j in 1:ncol(xs)) {
    xs[,j] <- xs[,j] / scale[j]
    xtests[,j] <- xtests[,j] / scale[j]
  }
  
  fit.time.end <- Sys.time()
  pred.time.start <- Sys.time()
  ## laGP analysis on scaled inputs distributed over 8 cores. 
  ## If you don't have 8 cores, change omp.threads to something smaller
  # out <- aGP(xs, y, xtests, d=list(start=1, max=20), g=g, verb=0)
  out.sep <- aGPsep(xs, y, xtests, d=list(start=1, max=20), g=g, verb=0,
                    end= min(50, Ntotal-1))
  pred.time.end <- Sys.time()
  
  # browser()
  list(mean=out.sep$mean, var=out.sep$var, n=nrow(Xs),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}

run_lagp_matt <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  require(laGP)
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  
  fit.time.start <- Sys.time()
  
  Xs <- x
  Ys <- y
  Xp <- xtest
  sdy <- sd(Ys)
  mny <- mean(Ys)
  ys <- (Ys-mny) / sdy
  # yp <- (Yp-mny) / sdy
  if(nrow(Xs)>80){
    gpisep <- newGPsep(Xs[1:80,], ys[1:80],d=rep(1,d), g=10^(-8),dK=TRUE)
  }else{
    gpisep <- newGPsep(Xs, ys,d=rep(1,d), g=10^(-8),dK=TRUE)
  }
  
  gpisep <-mleGPsep(gpisep)
  XslaGP<-Xs
  XplaGP<-Xp
  for(lcv in 1:d){
    XslaGP[,lcv] <- Xs[,lcv]/gpisep$d[lcv]
    XplaGP[,lcv] <- Xp[,lcv]/gpisep$d[lcv]
  }
  if(nrow(Xs)>80){
    gpi <- newGP(XslaGP[1:80,], ys[1:80],d=1, g=10^(-8),dK=TRUE)
  }else{
    gpi <- newGP(XslaGP, ys,d=1, g=10^(-8),dK=TRUE)
  }
  predsGloblaGP <- predGP(gpi, XslaGP)
  predpGloblaGP <- predGP(gpi, XplaGP)
  if(nrow(Xs)>80){
    predlaGP = aGP(XslaGP,ys-predsGloblaGP$mean,XplaGP,verb=0)
    fit.time.end <- Sys.time()
    pred.time.start <- Sys.time()
    
    predturn = list("mean"=mny+sdy*(predpGloblaGP$mean+predlaGP$mean),
                    "var"= sdy^2*predlaGP$var)
  }else{
    fit.time.end <- Sys.time()
    pred.time.start <- Sys.time()
    
    predturn = list("mean"=mny+sdy*(predpGloblaGP$mean),
                    "var"= sdy^2*diag(predpGloblaGP$Sigma))
  }
  pred.time.end <- Sys.time()
  
  # return(predturn)
  list(mean=predturn$mean, var=predturn$var, n=nrow(Xs),
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
  pred <- predict(mod, xtest, se.fit=T)
  pred.time.end <- Sys.time()
  list(mean=pred$fit, var=pred$se.fit^2, n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}

run_GPfit <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {#browser()
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  require(GPfit)
  fit.time.start <- Sys.time()
  mod <- GP_fit(x, y)
  fit.time.end <- Sys.time()
  pred.time.start <- Sys.time()
  pred <- predict(mod, xtest)
  pred.time.end <- Sys.time()
  list(mean=pred$Y_hat, var=pred$MSE, n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}

run_BASS <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  fit.time.start <- Sys.time()
  mod <- BASS::bass(xx=x, y=y)
  fit.time.end <- Sys.time()
  pred.time.start <- Sys.time()
  pred <- predict(mod, xtest)
  pred.time.end <- Sys.time()
  list(mean=colMeans(pred), var=apply(pred, 2, var), n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
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
      print('Nothing new to evaluate, hopefully !notdone')
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
      print('Nothing new to evaluate, hopefully !notdone')
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

run_CGGPsupponly <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed, correlation) {#browser()
  require("CGGP")
  if (!missing(seed)) {set.seed(seed)}
  if (Ntotal > 2000) {stop("CGGPsupponly can't run with more than 2000")}
  xsup <- lhs::maximinLHS(Ntotal, d)
  ysup <- apply(xsup, 1, f)
  fit.time.start <- Sys.time()
  sg <- CGGPcreate(d=d, batchsize=0, Xs=xsup, Ys=ysup, corr=correlation)
  fit.time.end <- Sys.time()
  # print(sg)
  Nevaluated <- nrow(sg$Xs)
  if (Nevaluated > Ntotal) {stop("CGGP has more points than Ntotal")}
  
  pred.time.start <- Sys.time()
  pred <- predict(sg, xtest)
  pred.time.end <- Sys.time()
  list(mean=pred$mean, var=pred$var,
       n=Nevaluated,
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}

run_CGGPoneshot <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed, correlation) {#browser()
  require("CGGP")
  if (!missing(seed)) {set.seed(seed)}
  # xsup <- lhs::maximinLHS(Nlhs, d)
  # ysup <- apply(xsup, 1, f)
  fit.time.start <- Sys.time()
  sg <- CGGPcreate(d=d, batchsize=Ntotal, corr=correlation)
  sg <- CGGPfit(sg, Y=apply(sg$design, 1, f))
  fit.time.end <- Sys.time()
  # print(sg)
  Nevaluated <- nrow(sg$design)
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
  outstats
}


expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

eg1 <- expand.grid(selection.method = c("UCB", "Greedy"),
                   correlation = c("CauchySQ", "Matern32", "PowerExp", "Cauchy", "CauchySQT"), stringsAsFactors=F)
eg2a <- expand.grid.df(eg1, data.frame(HandlingSuppData=c("Ignore", "Correct"), stringsAsFactors=F), data.frame(package="CGGPsupp", stringsAsFactors=F))
eg2b <- expand.grid.df(eg1, data.frame(HandlingSuppData="NA", stringsAsFactors=F), data.frame(package=c("CGGP"), stringsAsFactors=F))
eg2c <- expand.grid(selection.method="NA", correlation = c("CauchySQ", "Matern32", "PowerExp", "Cauchy", "CauchySQT"), HandlingSuppData="NA", package=c("CGGPoneshot", "CGGPsupponly"), stringsAsFactors=F)
eg2d <- data.frame(selection.method="NA", correlation="NA", HandlingSuppData="NA",
                   package=c("MRFA", "svm", "aGP", "aGP2", "laGP", "mlegp", "GPfit", "BASS"), stringsAsFactors=F)
eg3 <- rbind(eg2a, eg2b, eg2c, eg2d)

require("comparer")

excomp <- ffexp$new(
  eval_func = run_one,
  varlist = c("decentLHS", "run_CGGP", "run_CGGPoneshot", "run_CGGPsupp",
              "run_CGGPsupponly", "run_GPfit", "run_lagp", "run_lagp_bobby",
              "run_mlegp", "run_MRFA", "run_svm", "run_BASS", "run_lagp_matt"),
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"),
                d=c(3,6,7,8,10),
                row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
  # package=c("CGGPsupp", "CGPPoneshot", "CGGPsupponly", "CGGP",
  #           "MRFA", "svm", "aGP", "laGP", "mlegp", "GPfit"),
  psch=eg3,
  npd=c(10, 30, 100, 300, 1000, 3000, 10000),
  parallel=if (version$os =="linux-gnu") {TRUE} else {FALSE},
  parallel_cores = if (version$os =="linux-gnu") {34} else {3},
  replicate=1:10, #:5,
  folder_path= if (version$os =="linux-gnu") {"/home/collin/scratch/SGGP/scratch/ExternalComparison/ExComp5/"}
  else {"./scratch/ExternalComparison/ExComp5/"}
)
# Remove ones that can't do full size
package.name <- excomp$arglist$psch$package[excomp$rungrid$psch]
npd.excomp <- excomp$arglist$npd[excomp$rungrid$npd]
n.excomp <- npd.excomp * excomp$arglist$fd$d[excomp$rungrid$fd]
excomp$completed_runs[package.name == "CGGPsupponly" & n.excomp > 1200] <- TRUE
excomp$completed_runs[package.name == "laGP" & n.excomp > 2000] <- TRUE
excomp$completed_runs[package.name == "mlegp" & n.excomp > 500] <- TRUE
excomp$completed_runs[package.name == "GPfit" & n.excomp > 200] <- TRUE
# agp is giving errors
# excomp$completed_runs[package.name == "aGP"] <- TRUE
table(paste(package.name, excomp$completed_runs))
# plyr::ddply(data.frame(package.name, compruns=excomp$completed_runs), "package.name", function(d) {data.frame(done=sum(d$compruns), notdone=sum(!d$compruns))})
rbind(plyr::ddply(data.frame(package.name, compruns=excomp$completed_runs), "package.name", function(d) {data.frame(done=sum(d$compruns), notdone=sum(!d$compruns))}), data.frame(package.name="all", done=sum(excomp$completed_runs), notdone=sum(!excomp$completed_runs)))


# excomp$run_one(1752)
# excomp$run_one(7299)
# excomp$run_all(save_output = T, parallel = F, parallel_temp_save = T, run_order = "random")

# excomp$rungrid
# try because it gave delete error before, but shouldn't need it now
table(excomp$completed_runs)
try(excomp$recover_parallel_temp_save(delete_after = FALSE))
table(excomp$completed_runs)
# plyr::ddply(data.frame(package.name, compruns=excomp$completed_runs), "package.name", function(d) {data.frame(done=sum(d$compruns), notdone=sum(!d$compruns))})
rbind(plyr::ddply(data.frame(package.name, compruns=excomp$completed_runs), "package.name", function(d) {data.frame(done=sum(d$compruns), notdone=sum(!d$compruns))}), data.frame(package.name="all", done=sum(excomp$completed_runs), notdone=sum(!excomp$completed_runs)))
excomp$save_self()
excomp$run_all(parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE,
               write_start_files=T, write_error_files=T)
# Getting errors, run by package.name to see which is causing it
# excomp$parallel_cores <- 10
# excomp$run_all(to_run = which(!excomp$completed_runs & (package.name == "aGP")),
#                parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE,
#                write_start_files=T, write_error_files=T)
# excomp$run_all()

cat("Completed all runs in ExternalComparer4.R\n")

excomp$save_self()

cat("Saved self\n")

if (F) {
  excomp <- readRDS("./scratch/ExternalComparison/ExternalComparer4_mostlydone.rds")
  excompbass <- readRDS("./scratch/ExternalComparison/ExternalComparer4bass_quarterdone.rds")
  excomp$plot_run_times()
  plyr::dlply(excomp$outcleandf, "d")
  require('ggplot2')
  ecdf <- rbind(excomp$outcleandf[excomp$completed_runs & !is.na(excomp$outcleandf$package),],
                excompbass$outcleandf[excompbass$completed_runs & !is.na(excompbass$outcleandf$package),])
  ecdf$n <- ecdf$npd * ecdf$d
  ggplot(data=ecdf, mapping=aes(n, RMSE, color=package)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=ecdf, mapping=aes(n, score, color=package)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_x_log10()
  ggplot(data=ecdf[ecdf$package!="mlegp",], mapping=aes(n, score, color=package)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_x_log10()
  ggplot(data=ecdf, mapping=aes(n, CRPscore)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  ggplot(data=ecdf, mapping=aes(n, runtime)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  # saveRDS(excomp, "./scratch/ExternalComparison/ExComp1_completed.rds")
}
