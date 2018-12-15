#' Predict
#' 
#' Predict using SG with y values at xp?
#' Shouldn't y values already be stored in SG?
#'
#' @param xp x value to predict at
#' @param SGGP SG object
# ' @param y Observations for SG
# ' @param theta Correlation parameters
# ' @param logtheta Log of correlation parameters
# ' @param ... Don't use, just forces theta to be named
#'
#' @return Predicted mean values
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SGGP$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' SGGPpred(matrix(c(.1,.1,.1),1,3), SG=SG, y=y, theta=c(.1,.1,.1))
#' cbind(SGGPpred(SGGP$design, SG=SG, y=y, theta=c(.1,.1,.1))$mean, y) # Should be near equal
SGGPpred <- function(xp,SGGP) {
  # Require that you run SGGPfit
  
  # Cp is sigma(x_0) in paper, correlation vector between design points and xp
  Cp = matrix(1,dim(xp)[1],SGGP$ss)
  for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    V = SGGP$CorrMat(xp[,dimlcv], SGGP$xb, SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    Cp = Cp*V[,SGGP$designindex[,dimlcv]]
  }
  MSE_v = array(0, c(SGGP$d, SGGP$maxgridsize,dim(xp)[1]))
  for (dimlcv in 1:SGGP$d) {
    MSE_v[dimlcv, 1,] = 1
  }
  for (dimlcv in 1:SGGP$d) {
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      MSE_v[dimlcv, levellcv+1,] = SGGP_internal_MSEpredcalc(xp[,dimlcv],SGGP$xb[1:SGGP$sizest[levellcv]],SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],CorrMat=SGGP$CorrMat)
      MSE_v[dimlcv, levellcv+1,] = pmin(MSE_v[dimlcv, levellcv+1,], MSE_v[dimlcv, levellcv,])
    }
  }
  
  ME_t = prod(MSE_v[,1,],1)
  for (blocklcv in 1:SGGP$uoCOUNT) {
    ME_v = rep(1,dim(xp)[1])
    for (dimlcv in 1:SGGP$d) {
      levelnow = SGGP$uo[blocklcv,dimlcv]
      ME_v = ME_v*(MSE_v[dimlcv,1,]-MSE_v[dimlcv,levelnow+1,])
    }
    ME_t = ME_t-SGGP$w[blocklcv]*ME_v
  }
  
  
  # Return list with mean and var predictions
  if(is.vector(SGGP$pw)){
    GP = list("mean" = (SGGP$mu+Cp%*%SGGP$pw), "var"=SGGP$sigma2MAP[1]*ME_t)
  }else{
    if(length(SGGP$sigma2MAP)==1){
      GP = list("mean" = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+(Cp%*%SGGP$pw)%*%(SGGP$M)), "var"=as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%(SGGP$sigma2MAP)%*%(SGGP$M))))
      
    }else{
    GP = list("mean" = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+(Cp%*%SGGP$pw)%*%(SGGP$M)), "var"=as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%diag(SGGP$sigma2MAP)%*%(SGGP$M))))
    }
  }
  
  return(GP)
}

#' Calculate MSE over one dimension
#'
#' @param xp Points at which to calculate MSE
#' @param xl Levels along dimension, vector???
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#' for vector of 1D points.
#'
#' @return MSE predictions
#' @export
#'
#' @examples
#' MSEpred_calc(c(.4,.52), c(0,.25,.5,.75,1), theta=.1, nugget=1e-5,
#'              CorrMat=CorrMatMatern32,
#'              diag_corrMat=diag_corrMatMatern32)
SGGP_internal_MSEpredcalc <- function(xp,xl,theta,CorrMat) {
  S = CorrMat(xl, xl, theta)
  n = length(xl)
  cholS = chol(S)
  
  Cp = CorrMat(xp, xl, theta)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  
  MSE_val = 1 - rowSums(t(CiCp)*((Cp)))
  return(MSE_val)
}