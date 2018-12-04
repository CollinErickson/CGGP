#' Calculate predictive weights for SGGP
#' 
#' Predictive weights are Sigma^{-1}*y in standard GP.
#' This calculation is much faster since we don't need to
#' solve the full system of equations.
#'
#' @param SG SGGP object
#' @param y Measured values for SG$design
#' @param logtheta Log of correlation parameters
#' @param return_lS Should lS be returned?
#'
#' @return Vector with predictive weights
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' calculate_pw_C(SG=SG, y=y, logtheta=c(-.1,.1,.3))
calculate_pw_C <- function(SG, y, theta, return_lS=FALSE) {
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max value of all blocks
  # Now going to store choleskys instead of inverses for stability
  #CiS = list(matrix(1,1,1),Q*SG$d) # A list of matrices, Q for each dimension
  
  cholS = list(matrix(1,1,1),Q*SG$d) # To store choleskys
  lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Save log determinant of matrices
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over each possible needed correlation matrix
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      Sstuff = SG$CorrMat(Xbrn, Xbrn , theta[(lcv2-1)*SG$numpara+1:SG$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
      solvetry <- try({
        cS = chol(S)
        cholS[[(lcv2-1)*Q+lcv1]]= cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      })
      if (inherits(solvetry, "try-error")) {return(Inf)}
      lS[lcv1, lcv2] = 2*sum(log(diag(cS)))
    }
  }
  
  if(!is.matrix(y)){
    pw = rep(0, length(y)) # Predictive weight for each measured point
    # Loop over blocks selected
    gg = (1:SG$d-1)*Q
    for (lcv1 in 1:SG$uoCOUNT) {
      IS = SG$dit[lcv1, 1:SG$gridsizet[lcv1]];
      B = y[IS]
      rcpp_kronDBS(unlist(cholS[gg+SG$uo[lcv1,]]), B, SG$gridsizest[lcv1,])
      pw[IS] = pw[IS]+SG$w[lcv1] * B
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
    gg = (1:SG$d-1)*Q
    for (lcv1 in 1:SG$uoCOUNT) {
      IS = SG$dit[lcv1, 1:SG$gridsizet[lcv1]];
      VVV1 = unlist(cholS[gg+SG$uo[lcv1,]]);
      VVV2 = SG$gridsizest[lcv1,];
      for(lcv2 in 1:numout){
        B = y[IS,lcv2]
        rcpp_kronDBS(VVV1, B, VVV2)
        pw[IS,lcv2] = pw[IS,lcv2]+SG$w[lcv1] * B
      }
    }
    if (return_lS) {
      return(list(pw=pw, lS=lS))
    }else{
      return(pw)
    }
  }
}
#' 
#' Calculate derivative of pw
#'
#' @inheritParams calculate_pw
#' @param return_dlS Should dlS be returned?
#'
#' @return derivative matrix of pw with respect to logtheta
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' calculate_pw_and_dpw_C(SG=SG, y=y, logtheta=c(-.1,.1,.3))
calculate_pw_and_dpw_C <- function(SG, y, theta, return_lS=FALSE) {
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*SG$d) # To store choleskys
  dMatdtheta = list(matrix(1,1,1),Q*SG$d)
  
  if(return_lS){
    lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Save log determinant of matrices
    dlS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$numpara*SG$d)
  }
  
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over depth of each dim
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = SG$CorrMat(Xbrn, Xbrn , theta[(lcv2-1)*SG$numpara+1:SG$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      
      cS = chol(S)
      
      cholS[[(lcv2-1)*Q+lcv1]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      dMatdtheta[[(lcv2-1)*Q+lcv1]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
      for(lcv3 in 1:SG$numpara){
        dMatdtheta[[(lcv2-1)*Q+lcv1]][1:nv,nv*(lcv3-1)+1:nv] = t(dMatdtheta[[(lcv2-1)*Q+lcv1]][1:nv,nv*(lcv3-1)+1:nv])
      }
      if(return_lS){
        lS[lcv1, lcv2] = 2*sum(log(diag(cS)))
        for(lcv3 in 1:SG$numpara){
          dlS[lcv1, SG$numpara*(lcv2-1)+lcv3] = -sum(diag(dMatdtheta[[(lcv2-1)*Q+lcv1]][1:nv,nv*(lcv3-1)+1:nv]))
        }
      }
    }
  }
  
  pw = rep(0, length(y)) # predictive weights
  dpw = matrix(0, nrow = SG$numpara*SG$d, ncol = length(y)) # derivative of predictive weights
  gg = (1:SG$d-1)*Q
  for (lcv1 in 1:SG$uoCOUNT) {
    IS = SG$dit[lcv1, 1:SG$gridsizet[lcv1]];
    B = SG$w[lcv1]*y[IS]
    dB = rcpp_gkronDBS(unlist(cholS[gg+SG$uo[lcv1,]]),unlist(dMatdtheta[gg+SG$uo[lcv1,]]), B, SG$gridsizest[lcv1,])
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



calculate_sigma2_and_dsigma2_C <- function(SG, y, theta, return_lS=FALSE) {
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*SG$d) # To store choleskys
  dMatdtheta = list(matrix(1,1,1),Q*SG$d)
  
  if(return_lS){
    lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Save log determinant of matrices
    dlS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$numpara*SG$d) 
  }
  
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over depth of each dim
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = SG$CorrMat(Xbrn, Xbrn , theta[(lcv2-1)*SG$numpara+1:SG$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      
      cS = chol(S)
      
      cholS[[(lcv2-1)*Q+lcv1]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      dMatdtheta[[(lcv2-1)*Q+lcv1]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
      for(lcv3 in 1:SG$numpara){
        dMatdtheta[[(lcv2-1)*Q+lcv1]][1:nv,nv*(lcv3-1)+1:nv] = t(dMatdtheta[[(lcv2-1)*Q+lcv1]][1:nv,nv*(lcv3-1)+1:nv])
      }
      if(return_lS){
        lS[lcv1, lcv2] = 2*sum(log(diag(cS)))
        for(lcv3 in 1:SG$numpara){
          dlS[lcv1, SG$numpara*(lcv2-1)+lcv3] = -sum(diag(dMatdtheta[[(lcv2-1)*Q+lcv1]][1:nv,nv*(lcv3-1)+1:nv]))
        }
      }
    }
  }
  
  if(is.matrix(y)){
    numout = dim(y)[2]
    sigma2 = rep(0,numout) # Predictive weight for each measured point
    dsigma2 = matrix(0,nrow=SG$numpara*SG$d,ncol=numout) # Predictive weight for each measured point
    gg = (1:SG$d-1)*Q
    for (lcv1 in 1:SG$uoCOUNT) {
      IS = SG$dit[lcv1, 1:SG$gridsizet[lcv1]];
      VVV1=unlist(cholS[gg+SG$uo[lcv1,]])
      VVV2=unlist(dMatdtheta[gg+SG$uo[lcv1,]])
      VVV3=SG$gridsizest[lcv1,]
      for(lcv2 in 1:numout){
        B0 = y[IS,lcv2]
        B = SG$w[lcv1]*B0
        dB = rcpp_gkronDBS(VVV1,VVV2,B,VVV3)
        dsigma2[,lcv2] = dsigma2[,lcv2] + as.vector(t(B0)%*%t(dB))/length(y)
        sigma2[lcv2] = sigma2[lcv2]  + t(B0)%*%B/length(y)
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
    dsigma2 = rep(0,nrow=SG$d) # Predictive weight for each measured point
    gg = (1:SG$d-1)*Q
    for (lcv1 in 1:SG$uoCOUNT) {
      IS = SG$dit[lcv1, 1:SG$gridsizet[lcv1]];
      B0 = y[IS]
      B = SG$w[lcv1]*B0
      dB = rcpp_gkronDBS(unlist(cholS[gg+SG$uo[lcv1,]]),unlist(dMatdtheta[gg+SG$uo[lcv1,]]), B, SG$gridsizest[lcv1,])
      
      dsigma2 = dsigma2 +t(B0)%*%t(dB)/length(y)
      
      sigma2 = sigma2 + t(B0)%*%B/length(y)
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

