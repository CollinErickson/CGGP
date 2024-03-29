SGGP_internal_predwithonlysupp <- function(SGGP, xp, fullBayesian=FALSE, theta=NULL, outdims=NULL) {
  if (!inherits(SGGP, "SGGP")) {
    stop("First argument to SGGP must be an SGGP object")
  }
  # Require that you run SGGPfit first
  if (is.null(SGGP$supplemented)) {
    stop("You must run SGGPfit on SGGP object before using SGGPpredict")
  }
  if (!SGGP$supplemented) {
    stop("Must be supplemented to use SGGPpred_supponly")
  }
  # We could check for design_unevaluated, maybe give warning?
  
  # Full Bayesian
  if (fullBayesian) {
    if (!is.null(theta)) {stop("Don't give in theta for fullBayesian")}
    preds <- lapply(1:SGGP$numPostSamples, 
                    function(i) {SGGP_internal_predwithonlysupp(SGGP, xp, theta=SGGP$thetaPostSamples[,i])})
    means <- sapply(preds, function(x) {x$mean})
    mn <- as.matrix(apply(means, 1, mean))
    vars <- sapply(preds, function(x) {x$var})
    # This is for normal mixture, need t mixture?
    vr <- apply(vars, 1, mean) + apply(means^2, 1, mean) - apply(means, 1, mean)^2
    GP <- list(mean=mn, var=vr)
    # p <- SGGPpred(xp, SGGP)
    # stripchart(data.frame(t(vars)))
    # stripchart(data.frame(t(p$var)), add=T, col=2, pch=4)
    # print(cbind(p$var, vr, vr / p$var))
    return(GP)
  }
  
  # If theta is given (for full Bayesian prediction), need to recalculate pw
  if (!is.null(theta) && length(theta)!=length(SGGP$thetaMAP)) {stop("Theta is wrong length")}
  if (!is.null(theta) && theta==SGGP$thetaMAP) {
    # If you give in theta=thetaMAP, set it to NULL to avoid recalculating.
    theta <- NULL
  }
  if (is.null(theta)) {
    thetaMAP <- SGGP$thetaMAP
    recalculate_pw <- FALSE
  } else {
    thetaMAP <- theta
    rm(theta)
    # pw <- SGGP_internal_calcpw(SGGP, SGGP$y, theta=thetaMAP)
    recalculate_pw <- TRUE
  }
  
  
  separateoutputparameterdimensions <- is.matrix(SGGP$thetaMAP)
  # nnn is numberofoutputparameterdimensions
  nnn <- if (separateoutputparameterdimensions) {
    ncol(SGGP$ys)
  } else {
    1
  }
  if (nnn > 1) {
    # meanall <- matrix(NaN, nrow(xp), ncol=nnn)
    meanall2 <- matrix(0, nrow(xp), ncol=ncol(SGGP$Ys))
    # varall <- matrix(NaN, nrow(xp), ncol=ncol(SGGP$Ys))
    tempvarall <- matrix(0, nrow(xp), ncol=ncol(SGGP$Ys))
  }
  
  if (!is.null(outdims) && (nnn==1 || any(abs(SGGP$M-diag(1,nrow(SGGP$M),ncol(SGGP$M))) > 1e-10))) {
    stop("outdims can only be given when multiple outputs, no PCA, and separate correlation parameters")
  }
  opd_values <- if (is.null(outdims)) {1:nnn} else {outdims}
  for (opdlcv in opd_values) {# 1:nnn) {
    thetaMAP.thisloop <- if (nnn==1) thetaMAP else thetaMAP[, opdlcv]
    if (!recalculate_pw) { # use already calculated
      # pw.thisloop <- if (nnn==1) SGGP$pw else SGGP$pw[,opdlcv]
      # sigma2MAP.thisloop <- if (nnn==1) SGGP$sigma2MAP else SGGP$sigma2MAP[opdlcv]
      sigma2MAP.thisloop <- SGGP$sigma2MAP
    } else { # recalculate pw and sigma2MAP
      stop("not imp 3209842")
      y.thisloop <- if (nnn==1) SGGP$y else SGGP$y[,opdlcv]
      pw.thisloop <- SGGP_internal_calcpw(SGGP, y.thisloop, theta=thetaMAP.thisloop)
      sigma2MAP.thisloop <- SGGP_internal_calcsigma2anddsigma2(SGGP=SGGP, y=y.thisloop,
                                                               theta=thetaMAP.thisloop,
                                                               return_lS=FALSE)$sigma2
      # if (length(sigma2MAP.thisloop) != 1) {stop("sigma2MAP not scalar #923583")}
      # sigma2MAP.thisloop <- sigma2MAP.thisloop[1,1] # Convert 1x1 matrix to scalar
      # It can be vector when there is multiple output
      sigma2MAP.thisloop <- as.vector(sigma2MAP.thisloop)
      rm(y.thisloop)
    }
    mu.thisloop <- if (nnn==1) SGGP$mu else SGGP$mu[opdlcv] # Not used for PCA, added back at end
    
    # Cp is sigma(x_0) in paper, correlation vector between design points and xp
    # Cp = matrix(1,dim(xp)[1],SGGP$ss)
    # for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    #   V = SGGP$CorrMat(xp[,dimlcv], SGGP$xb, thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    #   Cp = Cp*V[,SGGP$designindex[,dimlcv]]
    # }
    # MSE_v = array(0, c(SGGP$d, SGGP$maxlevel + 1,dim(xp)[1])) # Add 1 to maxlevel so it doesn't go outside of array size
    # for (dimlcv in 1:SGGP$d) {
    #   MSE_v[dimlcv, 1,] = 1
    # }
    # for (dimlcv in 1:SGGP$d) {
    #   for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
    #     MSE_v[dimlcv, levellcv+1,] = SGGP_internal_MSEpredcalc(xp[,dimlcv],
    #                                                            SGGP$xb[1:SGGP$sizest[levellcv]],
    #                                                            thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
    #                                                            CorrMat=SGGP$CorrMat)
    #     MSE_v[dimlcv, levellcv+1,] = pmin(MSE_v[dimlcv, levellcv+1,], MSE_v[dimlcv, levellcv,])
    #   }
    # }
    # 
    # ME_t = prod(MSE_v[,1,],1)
    # for (blocklcv in 1:SGGP$uoCOUNT) {
    #   ME_v = rep(1,dim(xp)[1])
    #   for (dimlcv in 1:SGGP$d) {
    #     levelnow = SGGP$uo[blocklcv,dimlcv]
    #     ME_v = ME_v*(MSE_v[dimlcv,1,]-MSE_v[dimlcv,levelnow+1,])
    #   }
    #   ME_t = ME_t-SGGP$w[blocklcv]*ME_v
    # }
    
    
    if (!SGGP$supplemented) {  
      stop("Can't be here 0293582")
      # # Return list with mean and var predictions
      # if(is.vector(pw.thisloop)){
      #   if (nnn == 1) {
      #     mean = (mu.thisloop+Cp%*%pw.thisloop)
      #     var=sigma2MAP.thisloop*ME_t
      #   }
      #   
      #   # With sepparout and PCA (or not), do this
      #   if (nnn > 1) {
      #     meanall2 <- meanall2 + outer(c(Cp%*%pw.thisloop), SGGP$M[opdlcv, ])
      #     
      #     # This variance calculation was wrong when using PCA with separate theta for each output dim.
      #     # var <- (as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%diag(sigma2MAP.thisloop)%*%(SGGP$M))))[,opdlcv]
      #     # print("DEFINITELY WRONG! Need to do something with var and M here. Did something, maybe this is right, maybe not?")
      #     
      #     # This should be correct variance. Needs to be tested better.
      #     # Pick out the current dimension, set other values to zero
      #     tempM <- SGGP$M
      #     tempM[-opdlcv,] <- 0
      #     tempsigma2.thisloop <- sigma2MAP.thisloop
      #     tempsigma2.thisloop[-opdlcv] <- 0
      #     tempvar <- (as.vector(ME_t)%*%t(diag(t(tempM)%*%diag(tempsigma2.thisloop)%*%(tempM))))
      #     # print((as.vector(ME_t)%*%t(diag(t(tempM)%*%diag(tempsigma2.thisloop)%*%(tempM)))) %>% c %>% summary)
      #   }
      # }else{ # y was a matrix, so PCA
      #   if(length(sigma2MAP.thisloop)==1){
      #     stop("When is it a matrix but sigma2MAP a scalar???")
      #     # mean = ( matrix(rep(mu.thisloop,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+
      #     #            (Cp%*%pw.thisloop)%*%(SGGP$M))
      #     # var=as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%(sigma2MAP.thisloop)%*%(SGGP$M)))
      #     
      #   }else{
      #     mean = ( matrix(rep(mu.thisloop,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+
      #                (Cp%*%pw.thisloop)%*%(SGGP$M))
      #     var=as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%diag(sigma2MAP.thisloop)%*%(SGGP$M)))
      #   }
      # }
    } else { # SGGP$supplemented is TRUE
      if (!recalculate_pw) {
        # pw_uppad.thisloop <- if (nnn==1) SGGP$pw_uppad else SGGP$pw_uppad[,opdlcv]
        supppw.thisloop <- if (nnn==1) SGGP$supppw else SGGP$supppw[,opdlcv]
        Sti.thisloop <- if (nnn==1) SGGP$Sti else SGGP$Sti[,,opdlcv]
      } else {
        stop("Give theta in not implemented in SGGPpred. Need to fix sigma2MAP here too!")
      }
      
      
      Cps = matrix(1,dim(xp)[1],dim(SGGP$Xs)[1])
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(xp[,dimlcv], SGGP$Xs[,dimlcv], thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Cps = Cps*V
      }
      
      # yhatp = Cp%*%pw_uppad.thisloop + Cps%*%supppw.thisloop
      yhatp = Cps%*%supppw.thisloop
      
      # MSE_ps = list(matrix(0,dim(xp)[1],dim(SGGP$Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1)) 
      # for (dimlcv in 1:SGGP$d) {
      #   for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      #     MSE_ps[[(dimlcv)*SGGP$maxlevel+levellcv]] =(
      #       -SGGP_internal_postvarmatcalc(xp[,dimlcv],
      #                                     SGGP$Xs[,dimlcv],
      #                                     SGGP$xb[1:SGGP$sizest[levellcv]],
      #                                     thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
      #                                     CorrMat=SGGP$CorrMat))
      #   }
      # }
      
      # for (blocklcv in 1:SGGP$uoCOUNT) {
      #   ME_ps = matrix(1,nrow=dim(xp)[1],ncol=dim(SGGP$Xs)[1])
      #   for (dimlcv in 1:SGGP$d) {
      #     levelnow = SGGP$uo[blocklcv,dimlcv]
      #     ME_ps = ME_ps*MSE_ps[[(dimlcv)*SGGP$maxlevel+levelnow]]
      #   }
      #   Cps = Cps-SGGP$w[blocklcv]*(ME_ps)
      # }
      # ME_adj = rowSums((Cps%*%Sti.thisloop)*Cps)
      
      
      # ME_t = ME_t-ME_adj
      # Return list with mean and var predictions
      if(is.vector(supppw.thisloop)){
        if (nnn == 1) {
          mean = (SGGP$mu + yhatp)
          # print("You didn't fix var here")
          # browser()
          var=SGGP$sigma2MAP[1]* (1-diag(Cps %*% SGGP$Sti %*% t(Cps)))  # ME_t
          # could return cov mat here
        }
        
        # With sepparout and PCA (or not), do this
        if (nnn > 1) {
          warning("This won't work...")
          meanall2 <- meanall2 + outer(c(yhatp), SGGP$M[opdlcv,])
          leftvar <- if (is.null(SGGP$leftover_variance)) {0} else {SGGP$leftover_variance}
          # var <- (as.vector(ME_t)%*%t(leftvar+diag(t(SGGP$M)%*%diag(SGGP$sigma2MAP)%*%(SGGP$M))))[,opdlcv]
          # print("Need to do something with var and M here. Did something, maybe this is right, maybe not?")
          # Same fix as above
          tempM <- SGGP$M
          tempM[-opdlcv,] <- 0
          tempsigma2.thisloop <- sigma2MAP.thisloop
          tempsigma2.thisloop[-opdlcv] <- 0
          # tempvar <- (as.vector(ME_t)%*%t(leftvar+diag(t(tempM)%*%diag(tempsigma2.thisloop)%*%(tempM))))
          tempvar <- NaN
        }
        
      }else{ # supppw is matrix, so predicting multiple columns at once
        warning("This won't work either 2350729")
        if(length(SGGP$sigma2MAP)==1){
          stop("Does this ever happen? #952570")
          # mean = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ yhatp%*%(SGGP$M))
          # var=as.vector(ME_t)%*%t(SGGP$leftover_variance+diag(t(SGGP$M)%*%(SGGP$sigma2MAP)%*%(SGGP$M)))
        }else{
          mean = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ yhatp%*%(SGGP$M))
          # var=as.vector(ME_t)%*%t(SGGP$leftover_variance+diag(t(SGGP$M)%*%diag(SGGP$sigma2MAP)%*%(SGGP$M)))
          var <- NaN
        }
      }
      
      
    }
    # rm(Cp,ME_t, MSE_v, V) # Just to make sure nothing is carrying through
    # if (nnn > 1) {meanall[,opdlcv] <- mean}
    # if (nnn > 1) {varall[,opdlcv] <- var}
    if (nnn > 1) {tempvarall <- tempvarall + tempvar}
  }
  # If PCA values were calculated separately, need to do transformation on both before mu is added, then add mu back
  # if (nnn > 1) {meanall <- sweep(sweep(meanall,2,SGGP$mu) %*% SGGP$M,2,SGGP$mu, `+`)}
  if (nnn > 1) {meanall2 <- sweep(meanall2, 2, SGGP$mu, `+`)}
  
  if (nnn > 1) {
    GP <- list(mean=meanall2, var=tempvarall)
  } else {
    GP <- list(mean=mean, var=var)
  }
  # Check for negative variances, set them all to be a tiny number
  GP$var <- pmax(GP$var, .Machine$double.eps)
  return(GP)
}