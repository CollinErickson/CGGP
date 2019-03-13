#' Predict with only supplemental data
#'
#' @param CGGP CGGP object with supplemental data
#' @param xp Matrix of points to predict at
#' @param theta Correlation parameters
#' @param outdims Output dimensions to predict for
#'
#' @return Predictions
# @export
#' @noRd
#'
#' @examples
#' d <- 3
#' n <- 30
#' Xs <- matrix(runif(d*n), n, d)
#' Ys <- apply(Xs, 1, function(x){x[1]/(x[3]+.2)+exp(x[3])*cos(x[2]^2)})
#' cg <- CGGPcreate(d, Xs=Xs, Ys=Ys, batchsize=0)
#' xp <- matrix(runif(d*10), ncol=d)
#' predict(cg, xp)
CGGP_internal_predwithonlysupp <- function(CGGP, xp, theta=NULL, outdims=NULL) {
  if (!inherits(CGGP, "CGGP")) {
    stop("First argument to CGGP must be an CGGP object")
  }
  # Require that you run CGGPfit first
  if (is.null(CGGP$supplemented)) {
    stop("You must run CGGPfit on CGGP object before using CGGPpredict")
  }
  if (!CGGP$supplemented) {
    stop("Must be supplemented to use CGGPpred_supponly")
  }
  # We could check for design_unevaluated, maybe give warning?
  
  # If theta is given, need to recalculate pw
  if (!is.null(theta) && length(theta)!=length(CGGP$thetaMAP)) {stop("Theta is wrong length")}
  if (!is.null(theta) && theta==CGGP$thetaMAP) {
    # If you give in theta=thetaMAP, set it to NULL to avoid recalculating.
    theta <- NULL
  }
  if (is.null(theta)) {
    thetaMAP <- CGGP$thetaMAP
    recalculate_pw <- FALSE
  } else {
    thetaMAP <- theta
    rm(theta)
    # pw <- CGGP_internal_calcpw(CGGP, CGGP$y, theta=thetaMAP)
    recalculate_pw <- TRUE
  }
  
  
  separateoutputparameterdimensions <- is.matrix(CGGP$thetaMAP)
  # nnn is numberofoutputparameterdimensions
  nnn <- if (separateoutputparameterdimensions) {
    ncol(CGGP$ys)
  } else {
    1
  }
  if (nnn > 1) {
    # meanall <- matrix(NaN, nrow(xp), ncol=nnn)
    meanall2 <- matrix(0, nrow(xp), ncol=ncol(CGGP$Ys))
    # varall <- matrix(NaN, nrow(xp), ncol=ncol(CGGP$Ys))
    tempvarall <- matrix(0, nrow(xp), ncol=ncol(CGGP$Ys))
  }
  
  if (!is.null(outdims) && (nnn==1 || any(abs(CGGP$M-diag(1,nrow(CGGP$M),ncol(CGGP$M))) > 1e-10))) {
    stop("outdims can only be given when multiple outputs, no PCA, and separate correlation parameters")
  }
  opd_values <- if (is.null(outdims)) {1:nnn} else {outdims}
  for (opdlcv in opd_values) {# 1:nnn) {
    thetaMAP.thisloop <- if (nnn==1) thetaMAP else thetaMAP[, opdlcv]
    if (!recalculate_pw) { # use already calculated
      # pw.thisloop <- if (nnn==1) CGGP$pw else CGGP$pw[,opdlcv]
      # sigma2MAP.thisloop <- if (nnn==1) CGGP$sigma2MAP else CGGP$sigma2MAP[opdlcv]
      sigma2MAP.thisloop <- CGGP$sigma2MAP
    } else { # recalculate pw and sigma2MAP
      stop("not imp 3209842")
      y.thisloop <- if (nnn==1) CGGP$y else CGGP$y[,opdlcv]
      pw.thisloop <- CGGP_internal_calcpw(CGGP, y.thisloop, theta=thetaMAP.thisloop)
      sigma2MAP.thisloop <- CGGP_internal_calcsigma2anddsigma2(CGGP=CGGP, y=y.thisloop,
                                                               theta=thetaMAP.thisloop,
                                                               return_lS=FALSE)$sigma2
      # if (length(sigma2MAP.thisloop) != 1) {stop("sigma2MAP not scalar #923583")}
      # sigma2MAP.thisloop <- sigma2MAP.thisloop[1,1] # Convert 1x1 matrix to scalar
      # It can be vector when there is multiple output
      sigma2MAP.thisloop <- as.vector(sigma2MAP.thisloop)
      rm(y.thisloop)
    }
    mu.thisloop <- if (nnn==1) CGGP$mu else CGGP$mu[opdlcv] # Not used for PCA, added back at end
    
    # Cp is sigma(x_0) in paper, correlation vector between design points and xp
    # Cp = matrix(1,dim(xp)[1],CGGP$ss)
    # for (dimlcv in 1:CGGP$d) { # Loop over dimensions
    #   V = CGGP$CorrMat(xp[,dimlcv], CGGP$xb, thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
    #   Cp = Cp*V[,CGGP$designindex[,dimlcv]]
    # }
    # MSE_v = array(0, c(CGGP$d, CGGP$maxlevel + 1,dim(xp)[1])) # Add 1 to maxlevel so it doesn't go outside of array size
    # for (dimlcv in 1:CGGP$d) {
    #   MSE_v[dimlcv, 1,] = 1
    # }
    # for (dimlcv in 1:CGGP$d) {
    #   for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
    #     MSE_v[dimlcv, levellcv+1,] = CGGP_internal_MSEpredcalc(xp[,dimlcv],
    #                                                            CGGP$xb[1:CGGP$sizest[levellcv]],
    #                                                            thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
    #                                                            CorrMat=CGGP$CorrMat)
    #     MSE_v[dimlcv, levellcv+1,] = pmin(MSE_v[dimlcv, levellcv+1,], MSE_v[dimlcv, levellcv,])
    #   }
    # }
    # 
    # ME_t = prod(MSE_v[,1,],1)
    # for (blocklcv in 1:CGGP$uoCOUNT) {
    #   ME_v = rep(1,dim(xp)[1])
    #   for (dimlcv in 1:CGGP$d) {
    #     levelnow = CGGP$uo[blocklcv,dimlcv]
    #     ME_v = ME_v*(MSE_v[dimlcv,1,]-MSE_v[dimlcv,levelnow+1,])
    #   }
    #   ME_t = ME_t-CGGP$w[blocklcv]*ME_v
    # }
    
    
    if (!CGGP$supplemented) {  
      stop("Error in predwithonlysupp, there is no supp data #0293582")
    } else { # CGGP$supplemented is TRUE
      if (!recalculate_pw) {
        supppw.thisloop <- if (nnn==1) CGGP$supppw else CGGP$supppw[,opdlcv]
        Sti.thisloop <- if (nnn==1) CGGP$Sti else CGGP$Sti[,,opdlcv]
      } else {
        stop("Give theta in not implemented in CGGPpred. Need to fix sigma2MAP here too!")
      }
      
      
      Cps = matrix(1,dim(xp)[1],dim(CGGP$Xs)[1])
      for (dimlcv in 1:CGGP$d) { # Loop over dimensions
        V = CGGP$CorrMat(xp[,dimlcv], CGGP$Xs[,dimlcv], thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
        Cps = Cps*V
      }
      
      yhatp = Cps%*%supppw.thisloop

      # Return list with mean and var predictions
      if(is.vector(supppw.thisloop)){
        if (nnn == 1) {
          mean = (CGGP$mu + yhatp)
          var=CGGP$sigma2MAP[1]* (1-diag(Cps %*% CGGP$Sti %*% t(Cps)))  # ME_t
          # could return cov mat here
        }
        
        # With sepparout and PCA (or not), do this
        if (nnn > 1) {
          warning("This won't work...")
          meanall2 <- meanall2 + outer(c(yhatp), CGGP$M[opdlcv,])
          leftvar <- if (is.null(CGGP$leftover_variance)) {0} else {CGGP$leftover_variance}
          tempM <- CGGP$M
          tempM[-opdlcv,] <- 0
          tempsigma2.thisloop <- sigma2MAP.thisloop
          tempsigma2.thisloop[-opdlcv] <- 0
          # tempvar <- (as.vector(ME_t)%*%t(leftvar+diag(t(tempM)%*%diag(tempsigma2.thisloop)%*%(tempM))))
          tempvar <- NaN
        }
        
      }else{ # supppw is matrix, so predicting multiple columns at once
        warning("This won't work either 2350729")
        if(length(CGGP$sigma2MAP)==1){
          stop("Does this ever happen? #952570")
          # mean = ( matrix(rep(CGGP$mu,each=dim(xp)[1]), ncol=dim(CGGP$M)[2], byrow=FALSE)+ yhatp%*%(CGGP$M))
          # var=as.vector(ME_t)%*%t(CGGP$leftover_variance+diag(t(CGGP$M)%*%(CGGP$sigma2MAP)%*%(CGGP$M)))
        }else{
          mean = ( matrix(rep(CGGP$mu,each=dim(xp)[1]), ncol=dim(CGGP$M)[2], byrow=FALSE)+ yhatp%*%(CGGP$M))
          # var=as.vector(ME_t)%*%t(CGGP$leftover_variance+diag(t(CGGP$M)%*%diag(CGGP$sigma2MAP)%*%(CGGP$M)))
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
  # if (nnn > 1) {meanall <- sweep(sweep(meanall,2,CGGP$mu) %*% CGGP$M,2,CGGP$mu, `+`)}
  if (nnn > 1) {meanall2 <- sweep(meanall2, 2, CGGP$mu, `+`)}
  
  if (nnn > 1) {
    GP <- list(mean=meanall2, var=tempvarall)
  } else {
    GP <- list(mean=mean, var=var)
  }
  # Check for negative variances, set them all to be a tiny number
  GP$var <- pmax(GP$var, .Machine$double.eps)
  return(GP)
}
