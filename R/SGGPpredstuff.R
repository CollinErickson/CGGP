#' ????????????
#'
#' @param xp Points at which to calculate MSE
#' @param xl Levels along dimension, vector???
#' @param theta Correlation parameters
#' @param logtheta Log of correlation parameters
#' @param nugget Nugget to add to diagonal of correlation matrix
#' @param ... Don't use, just forces theta to be named
#'
#' @return MSE predictions
#' @export
#'
#' @examples
#' MSEpred_calc(c(.4,.52), c(0,.25,.5,.75,1), theta=.1, nugget=1e-5)
MSEpred_calc <- function(xp,xl, ..., logtheta, theta, nugget) {
  if (missing(theta)) {theta <- exp(logtheta)}
  S = CorrMat(xl, xl, theta=theta)
  S = S + nugget*diag(nrow(S))
  #t = exp(theta)
  n = length(xl)
  Ci = solve(S)
  
  Cp = CorrMat(xp, xl, theta=theta)
  
  MSE_val = 1-rowSums((Cp%*%Ci)*((Cp)))
  return(MSE_val)
  
}



#' Predict
#' 
#' Predict using SG with y values at xp?
#' Shouldn't y values already be stored in SG?
#'
#' @param xp x value to predict at
#' @param SG SG object
#' @param y Observations for SG
#' @param theta Correlation parameters
#' @param logtheta Log of correlation parameters
#' @param ... Don't use, just forces theta to be named
#'
#' @return Predicted mean values
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' SGGPpred(matrix(c(.1,.1,.1),1,3), SG=SG, y=y, theta=c(.1,.1,.1))
#' cbind(SGGPpred(SG$design, SG=SG, y=y, theta=c(.1,.1,.1))$mean, y) # Should be near equal
SGGPpred <- function(xp,SG, y, ..., logtheta, theta) {
  if (missing(theta)) {theta <- exp(logtheta)}
  # Center outputs
  my = mean(y)
  y = y-my
  
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max level of evaluated blocks
  CiS = list(matrix(1,1,1),Q*SG$d) # Store correlation matrices
  # Loop over dimensions and possible levels
  for (lcv2 in 1:SG$d) {
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]] # Get x's
      Xbrn = Xbrn[order(Xbrn)] # Sort them
      S = CorrMat(Xbrn, Xbrn , theta=theta[lcv2]) # Calculate corr mat
      S = S + SG$nugget*diag(nrow(S))
      CiS[[(lcv2-1)*Q+lcv1]] = solve(S) # Store inversion
      #print(eigen(S)$val)
      #print(det(S))
    }
  }
  
  
  pw = rep(0, length(y)) # ????????????
  # Loop over blocks
  for (lcv1 in 1:SG$uoCOUNT) {
    #narrowd = which(SG$uo[lcv1,] > 1.5) # Dimensions past first, THIS LINE IS REDUNDANT, CBE removed it
    Ci = 1
    if (lcv1 == 1) { # First block is all 1's
      narrowd = 1
    } else{ # All others have at least one
      narrowd = which(SG$uo[lcv1,] > 1.5)
    }
    for (e in narrowd) { # Multiply out to get full matrix
      Ci = kronecker(Ci, CiS[[((e-1)*Q+SG$uo[lcv1,e])]])
    }
    # Update predicted value
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * Ci %*% y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
  }
  sigma_hat = t(y) %*% pw / length(y)
  
  # Cp is sigma(x_0) in paper, correlation vector between design points and xp
  Cp = matrix(1,dim(xp)[1],SG$ss)
  for (e in 1:SG$d) { # Loop over dimensions
    Cp = Cp*CorrMat(xp[,e], SG$design[,e], theta=theta[e]) # Multiply correlation from each dimension
  }
  
  #    Cact = matrix(1,SG$ss,SG$ss)
  #    for (e in 1:SG$d) {
  #      Cact = Cact*CorrMat(SG$design[,e], SG$design[,e], theta[e])
  #    }
  
  
  MSE_v = array(0, c(SG$d, 8,dim(xp)[1]))
  for (lcv1 in 1:SG$d) {
    for (lcv2 in 1:8) {
      MSE_v[lcv1, lcv2,] = abs(MSEpred_calc(xp[,lcv1],SG$xb[1:SG$sizest[lcv2]], theta=theta[lcv1], nugget=SG$nugget))
      if (lcv2 > 1.5) {
        MSE_v[lcv1, lcv2,] = pmin(MSE_v[lcv1, lcv2,], MSE_v[lcv1, lcv2 - 1,])
      }
    }
  }
  
  
  ME_t = rep(1,dim(xp)[1])
  for (lcv1 in 1:SG$uoCOUNT) {
    ME_v = rep(1,dim(xp)[1])
    for (e in 1:SG$d) {
      levelnow = SG$uo[lcv1,e]
      if(levelnow > 1.5){
        ME_v = ME_v*(1-MSE_v[e,levelnow,])
      } else {
        ME_v = ME_v*(1-MSE_v[e,levelnow,])
      }
    }
    ME_t = ME_t-SG$w[lcv1]*ME_v
  }
  
  # Return list with mean and var predictions
  GP = list("mean" = (my+Cp %*% pw), "var"=sigma_hat[1]*ME_t)
  
  return(GP)
}
