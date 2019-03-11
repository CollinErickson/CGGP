#' Predict with SGGP object
#' 
#' Predict using SG with y values at xp?
#' Shouldn't y values already be stored in SG?
#'
#' @param xp x value to predict at
#' @param SGGP SG object
#' @param theta Leave as NULL unless you want to use a value other than thetaMAP.
#' Much slower.
#' @param outdims If multiple outputs fit without PCA and with separate
#' parameters, you can predict just for certain dimensions to speed it up.
#' Will leave other columns in the output, but they will be wrong.
#'
#' @return Predicted mean values
#' @export
#' @family SGGP core functions
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' SG <- SGGPfit(SG, Y=y)
#' SGGPpred(SG, matrix(c(.1,.1,.1),1,3))
#' cbind(SGGPpred(SG, SG$design)$mean, y) # Should be near equal
SGGPpred <- function(SGGP, xp, theta=NULL, outdims=NULL) {
  if (!inherits(SGGP, "SGGP")) {
    stop("First argument to SGGP must be an SGGP object")
  }
  # Require that you run SGGPfit first
  if (is.null(SGGP$supplemented)) {
    stop("You must run SGGPfit on SGGP object before using SGGPpredict")
  }
  if (SGGP$supplemented && (is.null(SGGP[["Y"]]) || length(SGGP$Y)==0)) {
    return(SGGP_internal_predwithonlysupp(SGGP=SGGP, xp=xp, theta=theta, outdims=outdims))
  }
  # We could check for design_unevaluated, maybe give warning?
  
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
  # nopd is numberofoutputparameterdimensions
  nopd <- if (separateoutputparameterdimensions) {
    ncol(SGGP$y)
  } else {
    1
  }
  if (nopd > 1) {
    # meanall <- matrix(NaN, nrow(xp), ncol=nopd)
    meanall2 <- matrix(0, nrow(xp), ncol=ncol(SGGP$Y))
    # varall <- matrix(NaN, nrow(xp), ncol=ncol(SGGP$Y))
    tempvarall <- matrix(0, nrow(xp), ncol=ncol(SGGP$Y))
  }
  
  if (!is.null(outdims) && nopd==1) {
    stop("outdims can only be given when multiple outputs and separate correlation parameters")
  }
  opd_values <- if (is.null(outdims)) {1:nopd} else {outdims}
  for (opdlcv in opd_values) {# 1:nopd) {
    thetaMAP.thisloop <- if (nopd==1) thetaMAP else thetaMAP[, opdlcv]
    if (!recalculate_pw) { # use already calculated
      pw.thisloop <- if (nopd==1) SGGP$pw else SGGP$pw[,opdlcv]
      # sigma2MAP.thisloop <- if (nopd==1) SGGP$sigma2MAP else SGGP$sigma2MAP[opdlcv]
      sigma2MAP.thisloop <- SGGP$sigma2MAP # Need this, not above, bc of PCA?
      cholS.thisloop <- if (nopd==1) SGGP$cholSs else SGGP$cholSs[[opdlcv]]
    } else { # recalculate pw and sigma2MAP
      y.thisloop <- if (nopd==1) SGGP$y else SGGP$y[,opdlcv]
      
      lik_stuff <- SGGP_internal_faststuff1(SGGP=SGGP, y=y.thisloop,theta=thetaMAP.thisloop)
      cholS.thisloop = lik_stuff$cholS
      sigma2MAP.thisloop <-  lik_stuff$sigma2
      pw.thisloop = lik_stuff$pw
      # if (length(sigma2MAP.thisloop) != 1) {stop("sigma2MAP not scalar #923583")}
      # sigma2MAP.thisloop <- sigma2MAP.thisloop[1,1] # Convert 1x1 matrix to scalar
      # It can be vector when there is multiple output
      sigma2MAP.thisloop <- as.vector(sigma2MAP.thisloop)
      rm(y.thisloop)
    }
    mu.thisloop <- if (nopd==1) SGGP$mu else SGGP$mu[opdlcv] # Not used for PCA, added back at end
    
    
    # Cp is sigma(x_0) in paper, correlation vector between design points and xp
    Cp = matrix(0,dim(xp)[1],SGGP$ss)
    GGGG = list(matrix(1,dim(xp)[1],length(SGGP$xb)),SGGP$d)
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      V = SGGP$CorrMat(xp[,dimlcv], SGGP$xb[1:SGGP$sizest[max(SGGP$uo[,dimlcv])]], thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],returnlogs=TRUE)
      GGGG[[dimlcv]] = exp(V)
      Cp = Cp+V[,SGGP$designindex[,dimlcv]]
    }
    Cp = exp(Cp)
    
    
    ME_t = matrix(1,dim(xp)[1],1)
    MSE_v = list(matrix(0,dim(xp)[1],2),(SGGP$d+1)*(SGGP$maxlevel+1)) 
    Q  = max(SGGP$uo[1:SGGP$uoCOUNT,])
    for (dimlcv in 1:SGGP$d) {
      for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
        Q  = max(SGGP$uo[1:SGGP$uoCOUNT,])
        gg = (dimlcv-1)*Q
        INDSN = 1:SGGP$sizest[levellcv]
        INDSN = INDSN[sort(SGGP$xb[1:SGGP$sizest[levellcv]],index.return = TRUE)$ix]
        MSE_v[[(dimlcv)*SGGP$maxlevel+levellcv]] = (SGGP_internal_postvarmatcalcfaster(GGGG[[dimlcv]],
                                                                                       c(),
                                                                                       as.matrix(cholS.thisloop[[gg+levellcv]]),
                                                                                       c(),
                                                                                       INDSN,
                                                                                       SGGP$numpara,
                                                                                       returndiag=TRUE))
      }
    }
    
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv]) > 0.5){
        ME_s = matrix(1,nrow=dim(xp)[1],1)
        for (dimlcv in 1:SGGP$d) {
          levelnow = SGGP$uo[blocklcv,dimlcv]
          ME_s = ME_s*MSE_v[[(dimlcv)*SGGP$maxlevel+levelnow]]
        }
        ME_t = ME_t-SGGP$w[blocklcv]*ME_s
      }
    }
    
    
    if (!SGGP$supplemented) {  
      
      # Return list with mean and var predictions
      if(is.vector(pw.thisloop)){
        if (nopd == 1) {
          mean = (mu.thisloop+Cp%*%pw.thisloop)
          var=sigma2MAP.thisloop*ME_t
        }
        
        # With sepparout and PCA (or not), do this
        if (nopd > 1) {
          # meanall2 <- meanall2 + outer(c(Cp%*%pw.thisloop), SGGP$M[opdlcv, ])
          meanall2[,opdlcv] <- c(Cp%*%pw.thisloop)
          
          # This should be correct variance. Needs to be tested better.
          # Pick out the current dimension, set other values to zero
          tempM <- diag(nopd) #SGGP$M
          tempM[-opdlcv,] <- 0
          tempsigma2.thisloop <- sigma2MAP.thisloop
          tempsigma2.thisloop[-opdlcv] <- 0
          tempvar <- (as.vector(ME_t)%*%t(diag(t(tempM)%*%diag(tempsigma2.thisloop)%*%(tempM))))
        }
      }else{ # y was a matrix
        if(length(sigma2MAP.thisloop)==1){
          stop("When is it a matrix but sigma2MAP a scalar???")
        }else{
          mean = (matrix(rep(mu.thisloop,each=dim(xp)[1]), ncol=ncol(SGGP$Y), byrow=FALSE) +
                     (Cp%*%pw.thisloop))
          var=as.vector(ME_t)%*%t(sigma2MAP.thisloop)
        }
      }
    } else { # SGGP$supplemented is TRUE
      if (!recalculate_pw) {
        pw_uppad.thisloop <- if (nopd==1) SGGP$pw_uppad else SGGP$pw_uppad[,opdlcv]
        supppw.thisloop <- if (nopd==1) SGGP$supppw else SGGP$supppw[,opdlcv]
        Sti.thisloop <- if (nopd==1) SGGP$Sti else SGGP$Sti[,,opdlcv]
      } else {
        stop("Give theta in not implemented in SGGPpred. Need to fix sigma2MAP here too!")
      }
      
      Cps = matrix(0,dim(xp)[1],dim(SGGP$Xs)[1])
      GGGG2 = list(matrix(0,nrow=dim(xp)[1],ncol=dim(SGGP$Xs)[1]),SGGP$d)
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(xp[,dimlcv], 
                         SGGP$Xs[,dimlcv], 
                         thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
                         returnlogs=TRUE)
        Cps = Cps+V
        
        V = SGGP$CorrMat(SGGP$Xs[,dimlcv], 
                         SGGP$xb[1:SGGP$sizest[max(SGGP$uo[,dimlcv])]], 
                         thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
                         returnlogs=TRUE)
        GGGG2[[dimlcv]] = exp(V)
      }
      Cps = exp(Cps)
      
      yhatp = Cp%*%pw_uppad.thisloop + Cps%*%supppw.thisloop
      
      
      MSE_ps = matrix(NaN,nrow=dim(SGGP$Xs)[1]*dim(xp)[1],ncol=(SGGP$d)*(SGGP$maxlevel))
      Q  = max(SGGP$uo[1:SGGP$uoCOUNT,])
      for (dimlcv in 1:SGGP$d) {
        gg = (dimlcv-1)*Q
        for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
          INDSN = 1:SGGP$sizest[levellcv]
          INDSN = INDSN[sort(SGGP$xb[1:SGGP$sizest[levellcv]],index.return = TRUE)$ix]
          REEALL= SGGP_internal_postvarmatcalcfasterasym6(GGGG[[dimlcv]],
                                                                                              GGGG2[[dimlcv]],
                                                                                              as.matrix(cholS.thisloop[[gg+levellcv]]),
                                                                                              INDSN)
          MSE_ps[,(dimlcv-1)*SGGP$maxlevel+levellcv]  =  as.vector(REEALL)
        }
      }
      
      
      Cps2 = as.vector(Cps)
      rcpp_fastmatclcr(SGGP$uo[1:SGGP$uoCOUNT,], SGGP$w[1:SGGP$uoCOUNT], MSE_ps, Cps2,SGGP$maxlevel)
      Cps = matrix(Cps2,ncol=dim(SGGP$Xs)[1] , byrow = FALSE)
      
      
      ME_adj = rowSums((Cps%*%Sti.thisloop)*Cps)
      
      
      ME_t = ME_t-ME_adj
      # Return list with mean and var predictions
      if(is.vector(pw.thisloop)){
        if (nopd == 1) {
          mean = (SGGP$mu+ yhatp)
          if (length(SGGP$sigma2MAP)>1) {warning("If this happens, you should fix var here")}
          # Should this be sigma2MAP.thisloop?
          var=SGGP$sigma2MAP[1]*ME_t
        }
        
        # With sepparout and PCA (or not), do this
        if (nopd > 1) {
          # meanall2 <- meanall2 + outer(c(yhatp), SGGP$M[opdlcv,])
          meanall2[,opdlcv] <- yhatp
          leftvar <- if (is.null(SGGP$leftover_variance)) {0} else {SGGP$leftover_variance}
          tempM <- diag(nopd) #SGGP$M
          tempM[-opdlcv,] <- 0
          tempsigma2.thisloop <- sigma2MAP.thisloop
          tempsigma2.thisloop[-opdlcv] <- 0
          tempvar <- (as.vector(ME_t)%*%t(leftvar+diag(t(tempM)%*%diag(tempsigma2.thisloop)%*%(tempM))))
        }
        
      }else{
        if(length(SGGP$sigma2MAP)==1){
          stop("This should never happen #952570")
        }else{
          leftvar <- if (is.null(SGGP$leftover_variance)) {0} else {SGGP$leftover_variance}
          mean = matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=ncol(SGGP$Y), byrow=FALSE)+ yhatp
          var=as.vector(ME_t)%*%t(leftvar+SGGP$sigma2MAP)
        }
      }
      
      
    }
    rm(Cp,ME_t, MSE_v, V) # Just to make sure nothing is carrying through
    if (nopd > 1) {tempvarall <- tempvarall + tempvar}
  }
  # If PCA values were calculated separately, need to do transformation on both before mu is added, then add mu back
  # if (nopd > 1) {meanall <- sweep(sweep(meanall,2,SGGP$mu) %*% SGGP$M,2,SGGP$mu, `+`)}
  if (nopd > 1) {meanall2 <- sweep(meanall2, 2, SGGP$mu, `+`)}
  
  if (nopd > 1) {
    GP <- list(mean=meanall2, var=tempvarall)
  } else {
    GP <- list(mean=mean, var=var)
  }
  # Check for negative variances, set them all to be a tiny number
  GP$var <- pmax(GP$var, .Machine$double.eps)
  return(GP)
}

#' Calculate MSE prediction along a single dimension
#'
#' @param xp Points at which to calculate MSE
#' @param xl Levels along dimension, vector???
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#'
#' @return MSE predictions
#' @export
#'
#' @examples
#' SGGP_internal_MSEpredcalc(c(.4,.52), c(0,.25,.5,.75,1), theta=c(.1,.2),
#'              CorrMat=SGGP_internal_CorrMatCauchySQ)
SGGP_internal_MSEpredcalc <- function(xp,xl,theta,CorrMat) {
  S = CorrMat(xl, xl, theta)
  n = length(xl)
  cholS = chol(S)
  
  Cp = CorrMat(xp, xl, theta)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  
  MSE_val = 1 - rowSums(t(CiCp)*((Cp)))
  return(MSE_val)
}



#' Calculate posterior variance, faster version
#'
#' @param GMat1 Matrix 1
#' @param GMat2 Matrix 2
#' @param cholS Cholesky factorization of S
#' @param INDSN Indices, maybe
#'
#' @return Variance posterior
#' @export
#'
#' @examples
#' SGGP_internal_postvarmatcalc(c(.4,.52), c(0,.25,.5,.75,1),
#'              xo=c(.11), theta=c(.1,.2,.3),
#'              CorrMat=SGGP_internal_CorrMatCauchySQT)
SGGP_internal_postvarmatcalcfasterasym6 <- function(GMat1,GMat2,cholS,INDSN) {
  
  CoinvC1o = backsolve(cholS,backsolve(cholS,t(GMat1[,INDSN]), transpose = TRUE))
  Sigma_mat =  (t(CoinvC1o)%*%(t(GMat2[,INDSN])))
  
  return(Sigma_mat )
  
}
