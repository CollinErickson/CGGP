

#' Calculate likelihood
#'
#' @param x Theta. Terrible variable name.
#' @param SG SGGP object
#' @param y Measured values of SG$design
#'
#' @return Likelihood
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' lik(c(.1,.1,.1), SG=SG, y=y)
lik <- function(logtheta, ..., SG, y) {
  #theta = x
  
  
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max value of all blocks
  CiS = list(matrix(1,1,1),Q*SG$d) # A list of matrices, Q for each dimension
  lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Save log determinant of matrices
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over each possible needed correlation matrix
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]] # xb are the possible points
      Xbrn = Xbrn[order(Xbrn)] # Sort them low to high, is this necessary? Probably just needs to be consistent.
      S = CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      CiS[[(lcv2-1)*Q+lcv1]] = solve(S)
      lS[lcv1, lcv2] = sum(log(eigen(S)$values))
    }
  }
  
  # Return Inf if theta is too large. Why????
  if (max(logtheta) >= (4 - 10 ^ (-6))) {
    return(Inf)
  } else{
    # We think pw is Sigma^{-1} * y
    pw = rep(0, length(y)) # For each point
    # Loop over blocks selected
    for (lcv1 in 1:SG$uoCOUNT) {
      narrowd = which(SG$uo[lcv1,] > 1.5) # Indices where block levels are > 1, THIS LINE IS REDUNDANT!!!!
      Ci = 1
      if (lcv1 == 1) {
        narrowd = 1 # First block has no levels above 1
      } else{
        narrowd = which(SG$uo[lcv1,] > 1.5) # Indices where block levels are > 1
      }
      # Expand correlation matrix on dimensions that are > 1
      for (e in narrowd) {
        Ci = kronecker(Ci, CiS[[((e-1)*Q+SG$uo[lcv1,e])]])
      }
      # Calculate something??? what is dit?
      pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
        SG$w[lcv1] * Ci %*% y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    }
    # This is not sigma hat?
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
    warning('should this be 3*theta???')
    logthetasqrt3 <- log(exp(logtheta)*sqrt(3))
    return(log(sigma_hat)+sum(logthetasqrt3^2)/length(y) + 1 / length(y) * lDet )
  }
  
}


#' Gradient of likelihood. Is it log likelihood?
#'
#' @param x Theta on normal scale
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
  CiS = list(matrix(1,1,1),Q*SG$d) # To store correlation matrices
  dCiS = list(matrix(1,1,1),Q*SG$d) # To store derivatives of corr mats
  lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store log det S
  dlS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store deriv of log det S
  
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over depth of each dim
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      S = CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      dS = dCorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      CiS[[(lcv2-1)*Q+lcv1]] = solve(S)
      dCiS[[(lcv2-1)*Q+lcv1]] = -CiS[[(lcv2-1)*Q+lcv1]]  %*% dS %*% CiS[[(lcv2-1)*Q+lcv1]] 
      lS[lcv1, lcv2] = sum(log(eigen(S)$values))
      dlS[lcv1, lcv2] = sum(eigen(CiS[[(lcv2-1)*Q+lcv1]] %*% dS)$values)
    }
  }
  
  pw = rep(0, length(y)) # ???
  
  dpw = matrix(0, nrow = length(y), ncol = SG$d) # ???
  for (lcv1 in 1:SG$uoCOUNT) {
    if (lcv1 == 1) {
      narrowd = 1
    } else{
      narrowd = which(SG$uo[lcv1,] > 1.5)
    }
    Ci = 1
    for (e in narrowd) {
      Ci = kronecker(Ci, CiS[[((e-1)*Q+SG$uo[lcv1,e])]])
    }
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * Ci %*% y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    
    for (e in narrowd) {
      Ci = 1
    for (e2 in narrowd) {
      if (e == e2) {
        Ci = kronecker(Ci,dCiS[[((e2-1)*Q+SG$uo[lcv1,e2])]])
      } else {
        Ci = kronecker(Ci, CiS[[((e2-1)*Q+SG$uo[lcv1,e2])]])
      }
    }
      
      dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]], e] = dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]], e] +
        SG$w[lcv1] * Ci %*% y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    
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
  warning("   this one also be theta * sqrt(3)??")
  logthetasqrt3 <- log(exp(logtheta)*sqrt(3))
 ddL = dsigma_hat / sigma_hat[1] + 2 / length(y) *logthetasqrt3 +  dlDet / length(y) 
 
 return(ddL)
}



#' Calculate theta MLE given data
#'
#' @param SG Sparse grid objects
#' @param y Output values calculated at SG$design
#' @param theta0 Initial theta. Is this on a log scale???
#' @param tol Tolerance for optimization
#'
#' @return theta MLE
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' thetaMLE(SG=SG, y=y)
logthetaMLE <- function(SG, y,..., logtheta0 = rep(0,SG$d),tol=1e-1) {
  x2 = optim(
    logtheta0,
    fn = lik,
    gr = glik,
    y = y - mean(y),
    SG = SG,
    method = "BFGS",
    hessian = FALSE,
    control = list(abstol = tol)
    # Is minimizing, default option of optim.
  )
  return(x2$par)
}
