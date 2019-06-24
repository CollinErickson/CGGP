# External comparisons, version 6
# This (v 6) is supposed to be for the first paper submission
# 5 is with updated grid size, extended in each dimension
# 6 adds BART. Also adding BASS into this file, last time had to do it separately

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
  
  list(mean=out.sep$mean, var=out.sep$var, n=nrow(xs),
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
  predsGloblaGP <- predGP(gpi, XslaGP, lite=T)
  predpGloblaGP <- predGP(gpi, XplaGP, lite=T)
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
                    "var"= sdy^2*(predpGloblaGP$s2)) # used to be diag(.$Sigma), gave error
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


run_svm <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
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

run_mlegp <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
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

run_GPfit <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
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

run_BART <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  fit.time.start <- Sys.time()
  mod <- BART::wbart(x.train=x, y.train=y, x.test=xtest)
  fit.time.end <- Sys.time()
  pred.time.start <- Sys.time()
  # pred <- predict(mod, xtest) # Doing pred inside fit
  pred.time.end <- Sys.time()
  list(mean=mod$yhat.test.mean, var=apply(mod$yhat.test, 2, var), n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}


run_CGGP <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed, selection.method, correlation) {
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
    # ni <- if (Nalready<1000) {200} else if (Nalready<10000) {500} else if (Nalready<20000) {2000} else {10000}
    ni <- if (Nalready<1000) {200} else if (Nalready<4000) {500} else if (Nalready<20000) {2000} else {10000}
    
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

run_CGGPsupp <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed, HandlingSuppData, selection.method, correlation) {
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
    # ni <- if (Nalready<1000) {200} else if (Nalready<10000) {500} else if (Nalready<20000) {2000} else {10000}
    ni <- if (Nalready<1000) {200} else if (Nalready<4000) {500} else if (Nalready<20000) {2000} else {10000}
    
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

run_CGGPsupponly <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed, correlation) {
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

run_CGGPoneshot <- function(Ntotal, Nappend, Nlhs, f, d, x, y, xtest, ytest, seed, correlation) {
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
                    f, d, npd, replicate) {
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
  } else if (package == "BART") {
    out <- run_BART(Ntotal=n, f=f, d=d, xtest=xtest)
  } else {
    stop(paste("Package", package, "not recognized"))
  }
  
  outstats <- CGGP::valstats(predmean=out[[1]], predvar=out[[2]],Yval=ytest) #, bydim=FALSE)
  if (out$n > n) {warning(paste("n too big for", package, n, f, d))}
  outstats$predtime <- out$pred.time
  outstats$fittime  <- out$fit.time
  # Forgot to add n_used, 
  outstats$n_used <- out$n
  outstats
}


expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

eg1 <- expand.grid(selection.method = c("UCB", "TS", "Greedy"),
                   correlation = c("CauchySQ", "Matern32", "PowerExp", "Cauchy", "CauchySQT"), stringsAsFactors=F)
eg2a <- expand.grid.df(eg1, data.frame(HandlingSuppData=c("Ignore", "Correct"), stringsAsFactors=F), data.frame(package="CGGPsupp", stringsAsFactors=F))
eg2b <- expand.grid.df(eg1, data.frame(HandlingSuppData="NA", stringsAsFactors=F), data.frame(package=c("CGGP"), stringsAsFactors=F))
eg2c <- expand.grid(selection.method="NA", correlation = c("CauchySQ", "Matern32", "PowerExp", "Cauchy", "CauchySQT"), HandlingSuppData="NA", package=c("CGGPoneshot", "CGGPsupponly"), stringsAsFactors=F)
eg2d <- data.frame(selection.method="NA", correlation="NA", HandlingSuppData="NA",
                   package=c("MRFA", "svm", "aGP", "aGP2", "laGP", "mlegp", "GPfit", "BASS", "BART"), stringsAsFactors=F)
eg3 <- rbind(eg2a, eg2b, eg2c, eg2d)

require("comparer")

excomp <- ffexp$new(
  eval_func = run_one,
  varlist = c("decentLHS", "run_CGGP", "run_CGGPoneshot", "run_CGGPsupp",
              "run_CGGPsupponly", "run_GPfit", "run_lagp", "run_lagp_bobby",
              "run_mlegp", "run_MRFA", "run_svm", "run_BASS", 
              "run_lagp_matt", "run_BART"),
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"),
                d=c(3,6,7,8,10),
                row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
  # package=c("CGGPsupp", "CGPPoneshot", "CGGPsupponly", "CGGP",
  #           "MRFA", "svm", "aGP", "laGP", "mlegp", "GPfit"),
  psch=eg3,
  npd=c(10, 30, 100, 300, 1000, 3000, 10000),
  parallel=if (version$os =="linux-gnu") {TRUE} else {FALSE},
  parallel_cores = if (version$os =="linux-gnu") {10} else {3},
  replicate=1:10, #:5,
  folder_path= if (version$os =="linux-gnu") {"/home/collin/scratch/SGGP/scratch/ExternalComparison/ExComp6/"}
  else {"./scratch/ExternalComparison/ExComp6/"}
)
# Remove ones that can't do full size
package.name <- excomp$arglist$psch$package[excomp$rungrid$psch]
npd.excomp <- excomp$arglist$npd[excomp$rungrid$npd]
n.excomp <- npd.excomp * excomp$arglist$fd$d[excomp$rungrid$fd]
excomp$completed_runs[package.name == "CGGPsupponly" & n.excomp > 1200] <- TRUE
excomp$completed_runs[package.name == "laGP" & n.excomp > 3000] <- TRUE
excomp$completed_runs[package.name == "mlegp" & n.excomp > 800] <- TRUE
excomp$completed_runs[package.name == "GPfit" & n.excomp > 400] <- TRUE
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
# excomp$run_all(parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE,
#                write_start_files=T, write_error_files=T)
# Getting errors, run by package.name to see which is causing it
# excomp$parallel_cores <- 10
if (F) {
  excomp$run_all(to_run = which(!excomp$completed_runs & (package.name == "BART")),
                 parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE,
                 write_start_files=T, write_error_files=T)
  # excomp$run_all()
}
cat("Completed all runs in ExternalComparer6.R\n")

excomp$save_self()

cat("Saved self\n")

# =============================================
# Load results and make plots of results ======
# =============================================
# For best plots, go to very bottom
if (F) {
  excomp <- readRDS("C:/Users/cbe117/Documents/GitHub/CGGP/scratch/ExternalComparison/ExComp6_almostall.rds")
  excomp$plot_run_times()
  # plyr::dlply(excomp$outcleandf, "d")
  require('ggplot2');require('dplyr');require('magrittr');
  table(excomp$completed_runs)
  excomp$rungrid2()[!excomp$completed_runs,]
  ecdf <- excomp$outcleandf[excomp$completed_runs & !is.na(excomp$outcleandf$package),]
  ecdf$n <- ecdf$npd * ecdf$d
  ggplot(data=ecdf, mapping=aes(n, RMSE, color=package)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=ecdf, mapping=aes(n, score, color=package)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_x_log10()
  ggplot(data=ecdf[ecdf$package!="mlegp",], mapping=aes(n, score, color=package)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_x_log10()
  ggplot(data=ecdf, mapping=aes(n, CRPscore)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  ggplot(data=ecdf, mapping=aes(n, runtime)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=ecdf, mapping=aes(n, fittime)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=ecdf, mapping=aes(n, predtime)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  # saveRDS(excomp, "./scratch/ExternalComparison/ExComp1_completed.rds")
  ggplot(data=ecdf %>% filter(package %in% c("CGGP","CGGPsupp", "CGGPoneshot")), mapping=aes(n, RMSE, color=correlation)) + geom_point() + facet_grid(f ~ interaction(package,correlation), scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=ecdf %>% filter(package %in% c("CGGP","CGGPsupp")), mapping=aes(n, RMSE, color=correlation)) + geom_point() + facet_grid(f ~ interaction(package,correlation), scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=ecdf, mapping=aes(n, RMSE, color=package)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  
  
  ##########################
  #### Plots for thesis ####
  ##########################
  
  #### Internal comparisons
  # First internal, CGGP vsCGGPsupp, compare correlations (only 3) and sel.method
  ggplot(data=ecdf %>% filter(package %in% c("CGGP", "CGGPsupp")), mapping=aes(n, RMSE, color=correlation, shape=correlation)) + geom_point(size=4) + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  # Only wingweight and OTL, supp, and Cauchy/M32/PowerExp
  ggplot(data=ecdf %>% filter(package=="CGGPsupp", correlation !="Cauchy", correlation!="CauchySQ", f=="wingweight" | f=="OTL_Circuit"), mapping=aes(n, RMSE, color=selection.method, shape=selection.method)) + geom_point(size=4) + facet_grid(f ~ correlation, scales="free_y") + scale_y_log10() + scale_x_log10()
  # Now compare CGGP vs CGGPsupp
  ggplot(data=ecdf %>% filter(package %in% c("CGGP","CGGPsupp"), correlation=="PowerExp", f=="wingweight" | f=="OTL_Circuit"), mapping=aes(n, RMSE, color=selection.method, shape=selection.method)) + geom_point(size=4) + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=ecdf %>% filter(package %in% c("CGGP","CGGPsupp"), correlation=="PowerExp", f=="wingweight" | f=="OTL_Circuit"), mapping=aes(n, score, color=selection.method, shape=selection.method)) + geom_point(size=4) + facet_grid(f ~ package, scales="free_y") + scale_x_log10()
  # Making plots better
  # First show that power exp is only reliable corr func
  incompdf1 <- ecdf %>% filter(package %in% c("CGGP","CGGPsupp"), correlation !="Cauchy", correlation!="CauchySQ", HandlingSuppData!="Ignore", f %in% c("wingweight","OTL_Circuit","borehole", "piston"), selection.method=="UCB") %>% rename(ICC=selection.method); if (nrow(incompdf1)!=4*3*10*2*7) {stop("incompdf1 wrong")}
  # incompdf1$ICC[incompdf1$ICC=="Greedy"] <- "ICC-MAP"; incompdf1$ICC[incompdf1$ICC=="UCB"] <- "ICC-UCB"
  incompdf1$package[incompdf1$package=="CGGP"] <- "No";incompdf1$package[incompdf1$package=="CGGPsupp"] <- "Yes"
  incompdf1$f[incompdf1$f=="borehole"] <- "Borehole";incompdf1$f[incompdf1$f=="OTL_Circuit"] <- "OTL circuit";incompdf1$f[incompdf1$f=="piston"] <- "Piston";incompdf1$f[incompdf1$f=="wingweight"] <- "Wing weight";
  
  pi1 <- ggplot(data=incompdf1, mapping=aes(n, RMSE, color=package, shape=package)) + labs(color="Supp. data?", shape="Supp. data?") + geom_point(size=4) + facet_grid(f ~ correlation, scales="free_y") + scale_y_log10() + scale_x_log10(); pi1
  pi2 <- ggplot(data=incompdf1, mapping=aes(n, CRPscore, color=package, shape=package)) + labs(color="Supp. data?", shape="Supp. data?") + geom_point(size=4) + facet_grid(f ~ correlation, scales="free_y") + scale_y_log10() + scale_x_log10()+ylab("CRP score"); pi2
  pi3 <- ggplot(data=incompdf1, mapping=aes(n, runtime, color=package, shape=package)) + labs(color="Supp. data?", shape="Supp. data?") + geom_point(size=4) + facet_grid(f ~ correlation, scales="free_y") + scale_y_log10() + scale_x_log10()+ylab("Run time (sec)"); pi3
  
  
  
  # SAVE IMAGES, set SAVEPLOT to FALSE to not save images
  SAVEPLOT <- T
  maybe_save <- function(filepath, p,folderpath="./scratch/thesis/", device='eps', width=4, height=4) {
    if (SAVEPLOT) {
      ggsave(paste0(folderpath, "/", filepath, ".", device), p, device=device, width=width, height=height, units="in")
    } else {p}
  }
  maybe_save("InternalCompRMSE_corr", device="eps", width=8, height=8, pi1)
  maybe_save("InternalCompCRPscore_corr", device="eps", width=8, height=8, pi2)
  maybe_save("InternalCompRuntime_corr", device="eps", width=8, height=8, pi3)
  
  
  # ===========================
  #### External comparisons
  # ===========================
  # First
  ggplot(data=ecdf %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS", "BART"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f %in% c("wingweight","OTL_Circuit","borehole")), mapping=aes(n, RMSE, color=package, shape=package)) + geom_point(size=4) + 
    facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18,19))
  ggplot(data=ecdf %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f %in% c("wingweight","OTL_Circuit","borehole")), mapping=aes(n, score, color=package, shape=package)) + geom_point(size=4) + 
    facet_grid(f ~ package, scales="free_y") + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18))
  ggplot(data=ecdf %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f %in% c("wingweight","OTL_Circuit","borehole")), mapping=aes(n, runtime, color=package, shape=package)) + geom_point(size=4) + 
    facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18))
  ecdf %>% group_by(package, correlation, selection.method)
  # shared plots
  ggplot(data=ecdf %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA")), mapping=aes(n, RMSE, color=package, shape=package)) + geom_point(size=4) + 
    facet_wrap(. ~ f) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18))
  ggplot(data=ecdf %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA")), mapping=aes(n, runtime, color=package, shape=package)) + geom_point(size=4) + 
    facet_wrap(. ~ f) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18))
  # Get mean lines for other packages
  ecdf_mean <- plyr::ddply(ecdf, c("package", "f", "n", "selection.method", "correlation", "HandlingSuppData"),
                           function(df) {data.frame(package=df$package[1],f=df$f[1], d=df$d[1], n=df$n[1],
                                                    selection.method=df$selection.method[1], correlation=df$correlation[1],
                                                    RMSEmeanlog=exp(mean(log(df$RMSE))),
                                                    CRPscoremeanlog=exp(mean(log(df$CRPscore))),
                                                    runtimemeanlog=exp(mean(log(df$runtime))), stringsAsFactors = F)})
  ggplot(data=ecdf_mean %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA")), mapping=aes(n, RMSEmeanlog, color=package, shape=package)) + geom_point(size=4) + 
    facet_wrap(. ~ f) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18))
  ggplot(data=ecdf_mean %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f  %in% c("borehole", "OTL_Circuit", "piston", "wingweight")), mapping=aes(n, RMSEmeanlog, color=package, shape=package)) + geom_point(size=4) + 
    facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18)) 
  ggplot() + geom_line(data=ecdf_mean %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f  %in% c("borehole", "OTL_Circuit", "piston", "wingweight")), mapping=aes(n, RMSEmeanlog, color=package), size=1) + 
    facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18))
  ggplot() + geom_line(data=ecdf_mean %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f  %in% c("borehole", "OTL_Circuit", "piston", "wingweight")), mapping=aes(n, RMSEmeanlog, color=package), size=1) + 
    facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(15,16,17,18))
  ecdf_mean_toplot <- ecdf_mean %>% filter(package %in% c("CGGPsupp","MRFA","aGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f  %in% c("borehole", "OTL_Circuit", "piston", "wingweight"))
  ecdf_mean_toplot$linesize <- 1 + 3*(ecdf_mean_toplot$package=="CGGPsupp")
  ggplot() +geom_line(data=ecdf_mean_toplot, mapping=aes(n, RMSEmeanlog, color=package, linetype=package), size=ecdf_mean_toplot$linesize) + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(17,15,16,18)) + labs(color="gn", linetype="gn", shape="gn")
  # Using smoother for mean
  # Put all packages on same plot, only 4 plots, 1 for each function, facet wrapped
  ecdf2 <- ecdf; ecdf2$f[ecdf2$f=="piston"] <- "Piston"; ecdf2$f[ecdf2$f=="OTL_Circuit"] <- "OTL circuit"; ecdf2$f[ecdf2$f=="wingweight"] <- "Wing weight"; ecdf2$f[ecdf2$f=="borehole"] <- "Borehole"
  ecdf2$package[ecdf$package=="aGP"] <- "laGP"; ecdf2 <- ecdf2 %>% filter(package!="CGGP");ecdf2$package[ecdf2$package=="CGGPsupp"] <- "CGGP"
  colnames(ecdf2)[colnames(ecdf2)=="package"] <- "Package"
  colnames(ecdf_mean_toplot)[colnames(ecdf_mean_toplot)=="package"] <- "Package"
  # Loess smoothing
  p1 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","MRFA","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, RMSE, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + stat_smooth(se=F, alpha=.6, geom='line') +scale_linetype_manual(values=c(2,1,3,4), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package"); p1
  p2 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, CRPscore, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + stat_smooth(se=F, alpha=.6, geom='line') +scale_linetype_manual(values=c(2,1,3,4), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package"); p2
  p3 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","MRFA","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, runtime, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + stat_smooth(se=F, alpha=.6, geom='line') +scale_linetype_manual(values=c(2,1,3,4), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package"); p3
  # Mean smoothing
  p1 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","MRFA","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, RMSE, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + scale_linetype_manual(values=c(2,1,3,4), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package")+ stat_summary(fun.y=mean, geom="line", alpha=.6); p1
  p2 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, CRPscore, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + scale_linetype_manual(values=c(2,1,3,4), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package") + stat_summary(fun.y=mean, geom="line", alpha=.6); p2
  p3 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","MRFA","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, runtime, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) +scale_linetype_manual(values=c(2,1,3,4), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package") + stat_summary(fun.y=mean, geom="line", alpha=.6); p3
  # Mean smoothing, but all normal lines (no dots/dashes). Use true means, not the log scale. For run times, only use means, no dots.
  p1 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","MRFA","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, RMSE, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + scale_linetype_manual(values=c(1,1,1,1), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package")+ stat_summary(fun.y=mean, geom="line", alpha=.6); p1
  p2 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, CRPscore, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + scale_linetype_manual(values=c(1,1,1,1), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package") + stat_summary(fun.y=mean, geom="line", alpha=.6); p2
  p3 <- ggplot(data=ecdf2 %>% filter(Package %in% c("CGGP","MRFA","laGP","BASS"), selection.method %in% c("NA","UCB"),correlation %in% c("PowerExp","NA"), f!='beambending') %>% 
                 mutate(pointsize=1+(Package=="CGGP")), mapping=aes(n, runtime, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package') + facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) +scale_linetype_manual(values=c(1,1,1,1), name='Package') + scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package") + stat_summary(fun.y=mean, geom="line", alpha=.6); p3
  
  maybe_save("ExternalCompRMSE_2", device="png", width=8, height=8, p1)
  maybe_save("ExternalCompCRPscore_2", device="png", width=8, height=8, p2)
  maybe_save("ExternalCompruntime_2", device="png", width=8, height=8, p3)
  
  # -------------------------------
  # Going to redo, make it cleaner
  # -------------------------------
  # Filter out only the 4 we want, on the four functions
  ecdf3 <- ecdf
  ecdf3 <- ecdf3 %>% filter(package %in% c("CGGPsupp", "MRFA", "aGP", "BASS"), selection.method %in% c("NA", "UCB"), correlation  %in% c("PowerExp", "NA"), HandlingSuppData %in% c("NA", "Correct"), f!="beambending")
  ecdf3$package[ecdf3$package=="CGGPsupp"] <- "CGGP"; ecdf3$package[ecdf3$package=="aGP"] <- "laGP"
  colnames(ecdf3)[colnames(ecdf3)=="package"] <- "Package"
  ecdf3$f[ecdf3$f=="piston"] <- "Piston"; ecdf3$f[ecdf3$f=="OTL_Circuit"] <- "OTL circuit"; ecdf3$f[ecdf3$f=="wingweight"] <- "Wing weight"; ecdf3$f[ecdf3$f=="borehole"] <- "Borehole"
  if (nrow(ecdf3) != 10*4*4*7 - 16) {stop("ecdf3 is wrong")} # Missing 1 BASS and 15 MRFA
  # Mean smoothing, but all normal lines (no dots/dashes). Use true means, not the log scale. For run times, only use means, no dots.
  p1 <- ggplot(data=ecdf3,
               mapping=aes(n, RMSE, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + scale_linetype_manual(values=c(1,1,1,1), name='Package') + 
    scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package") + 
    stat_summary(fun.y=mean, geom="line", alpha=.6) +
    geom_point(mapping=aes(size=Package)) + scale_size_manual(values=c(2,4,2,2), name='Package'); p1
  p2 <- ggplot(data=ecdf3 %>% filter(Package!= "MRFA"), 
               mapping=aes(n, CRPscore, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) + scale_linetype_manual(values=c(1,1,1,1), name='Package') + 
    scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package") + 
    stat_summary(fun.y=mean, geom="line", alpha=.6) + ylab("CRP score") +
    geom_point(mapping=aes(size=Package)) + scale_size_manual(values=c(2,4,2,2), name='Package'); p2
  p3 <- ggplot(data=ecdf3, 
               mapping=aes(n, runtime, color=Package, shape=Package, group=Package, linetype=Package, size=Package)) + 
    facet_wrap(. ~ f, scales="free") + scale_y_log10() + scale_x_log10() + 
    scale_shape_manual(values=c(15,17,16,18)) +scale_linetype_manual(values=c(1,1,1,1), name='Package') + 
    scale_color_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"), name="Package") + 
    stat_summary(fun.y=mean, geom="line", alpha=.6) + ylab("Run time (sec)") + 
    geom_point(mapping=aes(size=Package))+ scale_size_manual(values=c(2,4,2,2), name='Package'); p3
  # Save these
  maybe_save("ExternalCompRMSE_2", device="png", width=8, height=8, p1)
  maybe_save("ExternalCompCRPscore_2", device="png", width=8, height=8, p2)
  maybe_save("ExternalCompruntime_2", device="png", width=8, height=8, p3)
  
}
