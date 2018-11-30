#' ????????????
#'
#' @param xp Points at which to calculate MSE
#' @param xl Levels along dimension, vector???
#' @param theta Correlation parameters
#' @param logtheta Log of correlation parameters
#' @param nugget Nugget to add to diagonal of correlation matrix
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#' @param diag_corrMat Function that gives diagonal of correlation matrix
#' for vector of 1D points.
#' @param ... Don't use, just forces theta to be named
#'
#' @return MSE predictions
#' @export
#'
#' @examples
#' MSEpred_calc(c(.4,.52), c(0,.25,.5,.75,1), theta=.1, nugget=1e-5,
#'              CorrMat=CorrMatMatern32,
#'              diag_corrMat=diag_corrMatMatern32)
MSEpred_calc <- function(xp,xl, theta, CorrMat) {
  S = CorrMat(xl, xl, theta)
  n = length(xl)
  cholS = chol(S)
  Cp = CorrMat(xp, xl, theta)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  
  MSE_val = 1 - rowSums(t(CiCp)*((Cp)))
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
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' SGGPpred(matrix(c(.1,.1,.1),1,3), SG=SG, y=y, theta=c(.1,.1,.1))
#' cbind(SGGPpred(SG$design, SG=SG, y=y, theta=c(.1,.1,.1))$mean, y) # Should be near equal
SGGPpred <- function(xp,SG) {
  # Center outputs
  y = SG$y
  
  my = mean(y)
  y = y-my
  pw <- calculate_pw_C(SG,y,SG$theta) # Already subtracted mean(y) from y
  print(length(y))
  print(length(pw))
  
  sigma_hat = t(y) %*%pw / length(y)
  # Cp is sigma(x_0) in paper, correlation vector between design points and xp
  Cp = matrix(1,dim(xp)[1],SG$ss)
  for (e in 1:SG$d) { # Loop over dimensions
    #Cp = Cp*CorrMat(xp[,e], SG$design[,e], theta=theta[e]) # Multiply correlation from each dimension
    V = SG$CorrMat(xp[,e], SG$xb, SG$theta[(e-1)*SG$numpara+1:SG$numpara])
    Cp = Cp*V[,SG$designindex[,e]]
  }
  MSE_v = array(0, c(SG$d, SG$maxgridsize,dim(xp)[1]))
  for (lcv1 in 1:SG$d) {
    MSE_v[lcv1, 1,] = 1
  }
  for (lcv1 in 1:SG$d) {
    for (lcv2 in 1:max(SG$uo[1:SG$uoCOUNT,lcv1])) {
      MSE_v[lcv1, lcv2+1,] = MSEpred_calc(xp[,lcv1],SG$xb[1:SG$sizest[lcv2]],SG$theta[(lcv1-1)*SG$numpara+1:SG$numpara],CorrMat=SG$CorrMat)
      MSE_v[lcv1, lcv2+1,] = pmin(MSE_v[lcv1, lcv2+1,], MSE_v[lcv1, lcv2,])
    }
  }
  
  ME_t = prod(MSE_v[,1,],1)
  for (lcv1 in 1:SG$uoCOUNT) {
    ME_v = rep(1,dim(xp)[1])
    for (e in 1:SG$d) {
      levelnow = SG$uo[lcv1,e]
      ME_v = ME_v*(MSE_v[e,1,]-MSE_v[e,levelnow+1,])
    }
    ME_t = ME_t-SG$w[lcv1]*ME_v
  }
  
  
  # Return list with mean and var predictions
  GP = list("mean" = (my+Cp %*%pw), "var"=sigma_hat[1]*ME_t)
  
  return(GP)
}


#' Multivariate SGGP prediction
#' 
#' Predict output for an SGGP with multivariate output.
#'
#' @param yMV Output matrix, each row is for a design point
#' @inheritParams SGGPpred
#'
#' @return Matrix of predictions
#' @export
#'
#' @examples
#' 
#' SG <- SGcreate(d=3, batchsize=100)
#' y1 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' y2 <- apply(SG$design, 1, function(x){x[1]^1.3+.4*sin(6*x[2])+rnorm(1,0,.01)})
#' y <- cbind(y1, y2)
#' tx <- cbind(SGGPpredMV(SG$design, SG=SG, yMV=y, theta=c(.1,.1,.1))$mean, y)
#' tx # Columns 1 and 3, and 2 and 4 should be near equal
#' cor(tx)
SGGPpredMV <- function(xp,SG, yMV, ..., logtheta=SG$logtheta, theta) {
  
  p = dim(yMV)[2]
  
  meanMV = matrix(1,dim(xp)[1],p)
  varMV = matrix(1,dim(xp)[1],p)
  for(i in 1:p){
    GP1d = SGGPpred(xp,SG,as.vector(yMV[,i]),logtheta=logtheta, theta=theta)
    
    meanMV[,i] = GP1d$mean
    varMV[,i] = GP1d$var
  }
  
  # Return list with mean and var predictions
  GPMV = list("mean" = meanMV, "var"= varMV)
  
  return(GPMV)
}
