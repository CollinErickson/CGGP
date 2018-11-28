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
calculate_pw_C <- function(SG, y, logtheta, logsmooth, return_lS=FALSE) {
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max value of all blocks
  # Now going to store choleskys instead of inverses for stability
  #CiS = list(matrix(1,1,1),Q*SG$d) # A list of matrices, Q for each dimension
  CCS = list(matrix(1,1,1),Q*SG$d)
  lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Save log determinant of matrices
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over each possible needed correlation matrix
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]] # xb are the possible points
      Xbrn = Xbrn[order(Xbrn)] # Sort them low to high, is this necessary? Probably just needs to be consistent.
      S = SG$CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2], logsmooth = )
      diag(S) = diag(S) + SG$nugget
      # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
      solvetry <- try({
        cS = chol(S)
        CCS[[(lcv2-1)*Q+lcv1]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      })
      if (inherits(solvetry, "try-error")) {return(Inf)}
      lS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
    }
  }
  
  pw = rep(0, length(y)) # Predictive weight for each measured point
  # Loop over blocks selected
  for (lcv1 in 1:SG$uoCOUNT) {
    B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    rcpp_kronDBS(unlist(CCS[((1:SG$d-1)*Q+SG$uo[lcv1,1:SG$d])]), B, SG$gridsizest[lcv1,])
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * B
  }
  if (return_lS) {
    return(list(pw=pw, lS=lS))
  }
  pw
}

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
calculate_pw_and_dpw_C <- function(SG, y, logLS,logSM, return_lS=FALSE, return_dlS=FALSE) {
  warning("calculate_pw_and_dpw_C is not accurate yet!")
  
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*SG$d) # To store choleskys
  dSdLS = list(matrix(1,1,1),Q*SG$d)
  dSdSM = list(matrix(1,1,1),Q*SG$d)
  logdetS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store log det S
  logdetSLS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store deriv of log det S
  logdetSSM = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store deriv of log det S
  
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over depth of each dim
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      S = SG$CorrMat(Xbrn, Xbrn , logLS=logLS[lcv2], logSM=logSM[lcv2])
      diag(S) = diag(S) + SG$nugget
      cS = chol(S)
      
      cholS[[(lcv2-1)*Q+lcv1]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      dcholSdLS[[(lcv2-1)*Q+lcv1]] = -SG$dCorrMatdLS(Xbrn, Xbrn , logLS=logLS[lcv2], logSM=logSM[lcv2])
      dcholSdSM[[(lcv2-1)*Q+lcv1]] = -SG$dCorrMatdSM(Xbrn, Xbrn , logLS=logLS[lcv2], logSM=logSM[lcv2])
      
      logdetS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
      V = backsolve(CCS[[(lcv2-1)*Q+lcv1]],backsolve(CCS[[(lcv2-1)*Q+lcv1]],dS,transpose=TRUE));
      dlS[lcv1, lcv2] = sum(diag(V))
    }
  }
  
  pw = rep(0, length(y)) # predictive weights
  
  dpw = matrix(0, nrow = length(y), ncol = SG$d) # derivative of predictive weights
  
  for (lcv1 in 1:SG$uoCOUNT) {
    B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    dB = rcpp_gkronDBS(unlist(CCS[((1:SG$d-1)*Q+SG$uo[lcv1,1:SG$d])]), unlist(dCS[((1:SG$d-1)*Q+SG$uo[lcv1,1:SG$d])]), B, SG$gridsizest[lcv1,])
    dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],] = dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],] +
        SG$w[lcv1] * t(dB)
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * B
  }
  
  out <- list(pw=pw,
              dpw=dpw)
  if (return_lS) {
    out$lS <- lS
  }
  if (return_dlS) {
    out$dlS <- dlS
  }
  out
}
