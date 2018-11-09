

#' Calculate likelihood
#'
#' @param logtheta Log of correlation parameters
#' @param SG SGGP object
#' @param y Measured values of SG$design
#' @param ... Don't use, just forces theta to be named
#'
#' @return Likelihood
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' lik(c(.1,.1,.1), SG=SG, y=y)
lik <- function(logtheta, ..., SG, y) {
  
  # Return Inf if theta is too large. Why????
  if (max(logtheta) >= (4 - 10 ^ (-6))) {
    return(Inf)
  } else{
    
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
    
    # We think pw is Sigma^{-1} * y
    pw = rep(0, length(y)) # For each point
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
    sigma_hat = t(y) %*% pw / length(y)
    
    # Log determinant, keep a sum from smaller matrices
    lDet = 0
    
    # Calculate log det. See page 1586 of paper.
    # Loop over evaluated blocks
    for (lcv1 in 1:SG$uoCOUNT) {
      # Loop over dimensions
      for (lcv2 in 1:SG$d) {
        levelnow = SG$uo[lcv1, lcv2]
        # Add to log det when multiple points. It is zero when single point.
        if (levelnow > 1.5) {
          lDet = lDet + (lS[levelnow, lcv2] - lS[levelnow - 1, lcv2]) * (SG$gridsize[lcv1]) /
            (SG$gridsizes[lcv1, lcv2])
        }
      }
    }
    
    # Where does sum(theta^2) come from? Looks like regularization? Or from coordinate transformation
    # This next line is really wrong? The paranthese closes off the return before including the lDet.
    #ogthetasqrt3 <- log(exp(logtheta)*sqrt(3))
    return(log(c(sigma_hat))+sum(logtheta^2)/length(y) + 1 / length(y) * lDet )
  }
  
}


#' Gradient of likelihood. Is it log likelihood?
#'
#' @param logtheta Log of correlation parameters
#' @param SG SGGP object
#' @param y SG$design measured values
#' @param ... Don't use, just forces theta to be named
#'
#' @return Vector for gradient of likelihood w.r.t. x (theta)
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' glik(c(.1,.1,.1), SG=SG, y=y)
glik <- function(logtheta, ..., SG, y) {
  #theta = x
  
  
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
  
  pw = rep(0, length(y)) # ???
  
  dpw = matrix(0, nrow = length(y), ncol = SG$d) # ???
  
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
  sigma_hat = t(y) %*% pw / length(y)
  
  dsigma_hat = t(y) %*% dpw / length(y)
  
  lDet = 0 # Not needed for glik, only for lik
  
  dlDet = rep(0, SG$d) # Only needed for glik, not lik
  
  for (lcv1 in 1:SG$uoCOUNT) {
    for (lcv2 in 1:SG$d) {
      levelnow = SG$uo[lcv1, lcv2]
      if (levelnow > 1.5) {
        #lDet = lDet + (lS[levelnow, lcv2] - lS[levelnow - 1, lcv2]) * (SG$gridsize[lcv1]) /
        #  (SG$gridsizes[lcv1, lcv2]) # CBE added this, not needed for ddL, can use to return fn and gr at same time
        dlDet[lcv2] = dlDet[lcv2] + (dlS[levelnow, lcv2] - dlS[levelnow - 1, lcv2]) * (SG$gridsize[lcv1]) /
          (SG$gridsizes[lcv1, lcv2])
      }
    }
  }
  
  #logthetasqrt3 <- log(exp(logtheta)*sqrt(3))
  ddL = dsigma_hat / sigma_hat[1] + 2 / length(y) *logtheta +  dlDet / length(y) 
  
  return(ddL)
}



#' Calculate theta MLE given data
#'
#' @param SG Sparse grid objects
#' @param y Output values calculated at SG$design
#' @param logtheta0 Initial logtheta
#' @param tol Relative tolerance for optimization. Can't use absolute tolerance
#' since lik can be less than zero.
#' @param ... Don't use, just forces theta to be named
#'
#' @return theta MLE
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' logthetaMLE(SG=SG, y=y)
logthetaMLE <- function(SG, y,..., logtheta0 = rep(0,SG$d),tol=1e-4) {
  x2 = optim(
    logtheta0,
    fn = lik,
    gr = glik,
    lower = rep(-2, SG$d),
    upper = rep(3.9, SG$d),
    y = y - mean(y),
    SG = SG,
    method = "L-BFGS-B", #"BFGS",
    hessian = FALSE,
    control = list()#reltol=1e-4)#abstol = tol)
    # Is minimizing, default option of optim.
  )
  #return(pmin(2,x2$par)) # CBE adding this
  return(x2$par)
}


#' Likelihood for multivariate SGGP
#'
#' @param logtheta log of correlation parameters
#' @param SG SGGP object
#' @param yMV Matrix with output, number of columns is number of outputs
#'
#' @return Likelihood of logtheta
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y1 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' y2 <- apply(SG$design, 1, function(x){x[1]^1.3+.4*sin(6*x[2])+rnorm(1,0,.01)})
#' y <- cbind(y1, y2)
#' likMV(logtheta=c(.1,.2,.3), SG=SG, yMV=y)
likMV <- function(logtheta, SG, yMV) {
  p = dim(yMV)[2]
  
  logLikMV = 0
  for(i in 1:p){
    logLikMV = logLikMV+lik(logtheta, SG=SG, y=as.vector(yMV[,i])-mean(yMV[,i])) 
  }
  return(logLikMV)
}


#' Gradient of likelihood for multivariate SGGP
#'
#' @param logtheta log of correlation parameters
#' @param SG SGGP object
#' @param yMV Matrix with output, number of columns is number of outputs
#'
#' @return Vector, Gradient of likelihood of logtheta
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y1 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' y2 <- apply(SG$design, 1, function(x){x[1]^1.3+.4*sin(6*x[2])+rnorm(1,0,.01)})
#' y <- cbind(y1, y2)
#' glikMV(logtheta=c(.1,.2,.3), SG=SG, yMV=y)
glikMV <- function(logtheta, SG, yMV) {
  p = dim(yMV)[2]
  glogLikMV = rep(0,length(logtheta)) #p)
  for(i in 1:p){
    glogLikMV = glogLikMV+glik(logtheta, SG=SG, y=as.vector(yMV[,i])-mean(yMV[,i]))[1,] 
  }
  return(glogLikMV)
}

#' Find MLE of logtheta for multivariate output
#'
#' @param SG SG object
#' @param yMV Output matrix
#' @param ... Don't use
#' @param logtheta0 Initial values of logtheta for optimization
#' @param tol Optimization tolerance
#'
#' @return Vector, logtheta MLE
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y1 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' y2 <- apply(SG$design, 1, function(x){x[1]^1.3+.4*sin(6*x[2])+rnorm(1,0,.01)})
#' y <- cbind(y1, y2)
#' logthetaMLE(SG=SG, y=y1)
#' logthetaMLE(SG=SG, y=y2)
#' logthetaMLEMV(SG=SG, yMV=y)
logthetaMLEMV <- function(SG, yMV, ..., logtheta0 = rep(0,SG$d),tol=1e-4) {
  x2 = optim(
    logtheta0,
    fn = likMV,
    gr = glikMV,
    lower = rep(-2, SG$d),
    upper = rep(3.9, SG$d),
    SG = SG,
    yMV = yMV,
    method = "L-BFGS-B", #"BFGS",
    hessian = FALSE,
    control = list()#reltol=1e-4)#abstol = tol)
    # Is minimizing, default option of optim.
  )
  
  return(x2$par)
}
