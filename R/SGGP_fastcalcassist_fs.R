#' Calculate predictive weights for SGGP
#' 
#' Predictive weights are Sigma^{-1}*y in standard GP.
#' This calculation is much faster since we don't need to
#' solve the full system of equations.
#'
#' @param SGGP SGGP object
#' @param y Measured values for SGGP$design
#' @param theta Correlation parameters
#' @param return_lS Should lS be returned?
#'
#' @return Vector with predictive weights
#' @export
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' SGGP_internal_calcpw(SG=SG, y=y, theta=SG$thetaMAP)
SGGP_internal_calcpw <- function(SGGP, y, theta, return_lS=FALSE) {
  Q  = max(SGGP$uo[1:SGGP$uoCOUNT,]) # Max value of all blocks
  # Now going to store choleskys instead of inverses for stability
  #CiS = list(matrix(1,1,1),Q*SGGP$d) # A list of matrices, Q for each dimension
  
  cholS = list(matrix(1,1,1),Q*SGGP$d) # To store choleskys
  lS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$d) # Save log determinant of matrices
  # Loop over each dimension
  for (dimlcv in 1:SGGP$d) {
    # Loop over each possible needed correlation matrix
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      Xbrn = SGGP$xb[1:SGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = FALSE)
      Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
      solvetry <- try({
        cS = chol(S)
        cholS[[(dimlcv-1)*Q+levellcv]]= cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      })
      if (inherits(solvetry, "try-error")) {return(Inf)}
      lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
    }
  }
  
  if(!is.matrix(y)){
    pw = rep(0, length(y)) # Predictive weight for each measured point
    # Loop over blocks selected
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
      B = y[IS]
      rcpp_kronDBS(unlist(cholS[gg+SGGP$uo[blocklcv,]]), B, SGGP$gridsizest[blocklcv,])
      pw[IS] = pw[IS]+SGGP$w[blocklcv] * B
    }
    if (return_lS) {
      return(list(pw=pw, lS=lS))
    }else{
      return(pw)
    }
  }else{
    numout = dim(y)[2]
    pw = matrix(0,nrow=dim(y)[1],ncol=numout) # Predictive weight for each measured point
    # Loop over blocks selected
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
      VVV1 = unlist(cholS[gg+SGGP$uo[blocklcv,]]);
      VVV2 = SGGP$gridsizest[blocklcv,];
      for(outdimlcv in 1:numout){
        B = y[IS,outdimlcv]
        rcpp_kronDBS(VVV1, B, VVV2)
        pw[IS,outdimlcv] = pw[IS,outdimlcv]+SGGP$w[blocklcv] * B
      }
    }
    if (return_lS) {
      return(list(pw=pw, lS=lS))
    }else{
      return(pw)
    }
  }
}

#' Calculate derivative of pw
#'
#' @inheritParams SGGP_internal_calcpw
#' @param return_lS Should lS and dlS be returned?
#'
#' @return derivative matrix of pw with respect to logtheta
#' @export
#' @import Rcpp
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' SGGP_internal_calcpwanddpw(SG=SG, y=y, theta=SG$thetaMAP)
SGGP_internal_calcpwanddpw <- function(SGGP, y, theta, return_lS=FALSE) {
  Q  = max(SGGP$uo[1:SGGP$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*SGGP$d) # To store choleskys
  dMatdtheta = list(matrix(1,1,1),Q*SGGP$d)
  
  if(return_lS){
    lS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$d) # Save log determinant of matrices
    dlS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$numpara*SGGP$d)
  }
  
  # Loop over each dimension
  for (dimlcv in 1:SGGP$d) {
    # Loop over depth of each dim
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      Xbrn = SGGP$xb[1:SGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      
      cS = chol(S)
      
      cholS[[(dimlcv-1)*Q+levellcv]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      dMatdtheta[[(dimlcv-1)*Q+levellcv]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
      for(paralcv in 1:SGGP$numpara){
        dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv] = t(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv])
      }
      if(return_lS){
        lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
        for(paralcv in 1:SGGP$numpara){
          dlS[levellcv, SGGP$numpara*(dimlcv-1)+paralcv] = -sum(diag(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv]))
        }
      }
    }
  }
  
  pw = rep(0, length(y)) # predictive weights
  dpw = matrix(0, nrow = SGGP$numpara*SGGP$d, ncol = length(y)) # derivative of predictive weights
  gg = (1:SGGP$d-1)*Q
  for (blocklcv in 1:SGGP$uoCOUNT) {
    IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
    B = SGGP$w[blocklcv]*y[IS]
    dB = rcpp_gkronDBS(unlist(cholS[gg+SGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]]), B, SGGP$gridsizest[blocklcv,])
    dpw[,IS] = dpw[,IS] +dB
    pw[IS] = pw[IS] + B
  }
  dpw =t(dpw)
  out <- list(pw=pw,
              dpw=dpw)
  if (return_lS) {
    out$lS <- lS
    out$dlS <- dlS
  }
  
  out
}



SGGP_internal_calcsigma2 <- function(SGGP, y, theta, return_lS=FALSE) {
  Q  = max(SGGP$uo[1:SGGP$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*SGGP$d) # To store choleskys
  if(return_lS){
    lS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$d) # Save log determinant of matrices
  }
  
  # Loop over each dimension
  for (dimlcv in 1:SGGP$d) {
    # Loop over depth of each dim
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      Xbrn = SGGP$xb[1:SGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = FALSE)
      S = Sstuff
      cS = chol(S)
      
      cholS[[(dimlcv-1)*Q+levellcv]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      if(return_lS){
        lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
      }
    }
  }
  if(is.matrix(y)){
    numout = dim(y)[2]
    sigma2 = rep(0,numout) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
      VVV1=unlist(cholS[gg+SGGP$uo[blocklcv,]])
      VVV3=SGGP$gridsizest[blocklcv,]
      for(outdimlcv in 1:numout){
        B0 = y[IS,outdimlcv]
        B = (SGGP$w[blocklcv]/dim(y)[1])*B0
        rcpp_kronDBS(VVV1,B,VVV3)
        sigma2[outdimlcv] = sigma2[outdimlcv]  + t(B0)%*%B
      }
    }
    out <- list(sigma2=sigma2)
    if (return_lS) {
      out$lS <- lS
    }
  }else{
    sigma2 = 0 # Predictive weight for each measured point
    dsigma2 = rep(0,nrow=SGGP$d) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
      B0 = y[IS]
      B = (SGGP$w[blocklcv]/length(y))*B0
      rcpp_kronDBS(unlist(cholS[gg+SGGP$uo[blocklcv,]]),B, SGGP$gridsizest[blocklcv,])
      sigma2 = sigma2 + t(B0)%*%B
    }
    out <- list(sigma2=sigma2)
    if (return_lS) {
      out$lS <- lS
    }
  }
  return(out)
}

SGGP_internal_calcsigma2anddsigma2 <- function(SGGP, y, theta, return_lS=FALSE) {
  Q  = max(SGGP$uo[1:SGGP$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*SGGP$d) # To store choleskys
  dMatdtheta = list(matrix(1,1,1),Q*SGGP$d)
  
  if(return_lS){
    lS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$d) # Save log determinant of matrices
    dlS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$numpara*SGGP$d) 
  }
  
  # Loop over each dimension
  for (dimlcv in 1:SGGP$d) {
    # Loop over depth of each dim
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      Xbrn = SGGP$xb[1:SGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      
      cS = chol(S)
      
      cholS[[(dimlcv-1)*Q+levellcv]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      dMatdtheta[[(dimlcv-1)*Q+levellcv]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
      for(paralcv in 1:SGGP$numpara){
        dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv] = t(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv])
      }
      if(return_lS){
        lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
        for(paralcv in 1:SGGP$numpara){
          dlS[levellcv, SGGP$numpara*(dimlcv-1)+paralcv] = -sum(diag(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv]))
        }
      }
    }
  }
  
  if(is.matrix(y)){
    numout = dim(y)[2]
    sigma2 = rep(0,numout) # Predictive weight for each measured point
    dsigma2 = matrix(0,nrow=SGGP$numpara*SGGP$d,ncol=numout) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
      VVV1=unlist(cholS[gg+SGGP$uo[blocklcv,]])
      VVV2=unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]])
      VVV3=SGGP$gridsizest[blocklcv,]
      for(outdimlcv in 1:numout){
        B0 = y[IS,outdimlcv]
        B = (SGGP$w[blocklcv]/dim(y)[1])*B0
        dB = rcpp_gkronDBS(VVV1,VVV2,B,VVV3)
        dsigma2[,outdimlcv] = dsigma2[,outdimlcv] + as.vector(dB%*%B0)
        sigma2[outdimlcv] = sigma2[outdimlcv]  + sum(B0*B)
      }
    }
    out <- list(sigma2=sigma2,
                dsigma2=dsigma2)
    if (return_lS) {
      out$lS <- lS
      out$dlS <- dlS
    }
  }else{
    sigma2 = 0 # Predictive weight for each measured point
    dsigma2 = rep(0,nrow=SGGP$d) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
      B0 = y[IS]
      B = (SGGP$w[blocklcv]/length(y))*B0
      dB = rcpp_gkronDBS(unlist(cholS[gg+SGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]]), B, SGGP$gridsizest[blocklcv,])
      
      dsigma2 = dsigma2 +t(B0)%*%t(dB)
      sigma2 = sigma2 + t(B0)%*%B
    }
    out <- list(sigma2=sigma2,
                dsigma2=dsigma2)
    if (return_lS) {
      out$lS <- lS
      out$dlS <- dlS
    }
  }
  out
}

