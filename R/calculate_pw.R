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
#' calculate_pw(SG=SG, y=y, logtheta=c(-.1,.1,.3))
calculate_pw <- function(SG, y, logtheta, return_lS=FALSE) {
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
      S = SG$CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      diag(S) = diag(S) + SG$nugget
      # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
      solvetry <- try({
        #CiS[[(lcv2-1)*Q+lcv1]] = solve(S)
        CCS[[(lcv2-1)*Q+lcv1]] = chol(S)
      })
      if (inherits(solvetry, "try-error")) {return(Inf)}
      #lS[lcv1, lcv2] = sum(log(eigen(S)$values))
      lS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
    }
  }
  
  pw = rep(0, length(y)) # Predictive weight for each measured point
  # Loop over blocks selected
  for (lcv1 in 1:SG$uoCOUNT) {
    B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    for (e in SG$d:1) {
      if(SG$gridsizest[lcv1,e] > 1.5){
        B <- matrix(as.vector(B),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
        B <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B, transpose = TRUE))
        #B <-  solve(CS[[((e-1)*Q+SG$uo[lcv1,e])]],B)
        B <- t(B)
      }
      else{
        B = as.vector(B)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
      }
    }
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
#'
#' @return derivative matrix of pw with respect to logtheta
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' calculate_pw_and_dpw(SG=SG, y=y, logtheta=c(-.1,.1,.3))
calculate_pw_and_dpw <- function(SG, y, logtheta, return_lS=FALSE) {
  
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max level of all blocks
  # Now storing choleskys instead of inverses
  CCS = list(matrix(1,1,1),Q*SG$d) # To store choleskys
  dCS = list(matrix(1,1,1),Q*SG$d)
  #CiS = list(matrix(1,1,1),Q*SG$d) # To store correlation matrices
  #dCiS = list(matrix(1,1,1),Q*SG$d) # To store derivatives of corr mats
  lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store log det S
  dlS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store deriv of log det S
  dlS2 = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # ???? added for chols
  
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over depth of each dim
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      S = SG$CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      diag(S) = diag(S) + SG$nugget
      dS = SG$dCorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      
      CCS[[(lcv2-1)*Q+lcv1]] = chol(S)#-CiS[[(lcv2-1)*Q+lcv1]]  %*% 
      dCS[[(lcv2-1)*Q+lcv1]] = dS
      lS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
      V = solve(S,dS);
      dlS[lcv1, lcv2] = sum(diag(V))
    }
  }
  
  pw = rep(0, length(y)) # predictive weights
  
  dpw = matrix(0, nrow = length(y), ncol = SG$d) # derivative of predictive weights
  
  for (lcv1 in 1:SG$uoCOUNT) {
    
    B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    for (e in SG$d:1) {
      if(SG$gridsizest[lcv1,e] > 1.5){
        B <- matrix(as.vector(B),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
        B <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B, transpose = TRUE))
        B <- t(B)
      }
      else{
        B = as.vector(B)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
      }
    }
    
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * B
    
    
    B3 = B
    for (e in  SG$d:1) {
      if(SG$gridsizest[lcv1,e] > 1.5){
        B3 <- matrix(as.vector(B3),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
      }
      
      B2 = B3
      
      if(SG$gridsizest[lcv1,e] > 1.5){
        B3 <- t(B3)
      } else{
        B3 = as.vector(B3)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
      }
      
      
      if(SG$gridsizest[lcv1,e] > 1.5){
        B2 <- -dCS[[((e-1)*Q+SG$uo[lcv1,e])]]%*%B2
        B2 <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B2, transpose = TRUE))
        B2 = t(B2)
      }else{
        B2 =-as.vector(dCS[[((e-1)*Q+SG$uo[lcv1,e])]])*as.vector(B2)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
      }
      
      if(e>1.5){
        for (e2 in  (e-1):1) {
          if(SG$gridsizest[lcv1,e2] > 1.5){
            B2 <- matrix(as.vector(B2),SG$gridsizest[lcv1,e2],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e2])
            B2 <- t(B2)
          }
          else{
            B2= as.vector(B2)
          }
        }
      }
      
      dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],e] = dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],e] +
        SG$w[lcv1] * B2
    }
    
  }
  
  list(pw=pw,
       dpw=dpw)
}
