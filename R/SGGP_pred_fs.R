#' Predict with SGGP object
#' 
#' Predict using SG with y values at xp?
#' Shouldn't y values already be stored in SG?
#'
#' @param xp x value to predict at
#' @param SGGP SG object
#'
#' @return Predicted mean values
#' @export
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' SG <- SGGPfit(SG, Y=y)
#' SGGPpred(matrix(c(.1,.1,.1),1,3), SGGP=SG)
#' cbind(SGGPpred(SG$design, SG=SG)$mean, y) # Should be near equal
SGGPpred <- function(xp,SGGP) {
  # Require that you run SGGPfit first
  if (is.null(SGGP$supplemented)) {
    stop("You must run SGGPfit on SGGP object before using SGGPpredict")
  }
  # We could check for design_unevaluated, maybe give warning?
  
  # Cp is sigma(x_0) in paper, correlation vector between design points and xp
  Cp = matrix(1,dim(xp)[1],SGGP$ss)
  for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    V = SGGP$CorrMat(xp[,dimlcv], SGGP$xb, SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    Cp = Cp*V[,SGGP$designindex[,dimlcv]]
  }
  MSE_v = array(0, c(SGGP$d, SGGP$maxlevel + 1,dim(xp)[1])) # Add 1 to maxlevel so it doesn't go outside of array size
  for (dimlcv in 1:SGGP$d) {
    MSE_v[dimlcv, 1,] = 1
  }
  for (dimlcv in 1:SGGP$d) {
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      MSE_v[dimlcv, levellcv+1,] = SGGP_internal_MSEpredcalc(xp[,dimlcv],
                                                             SGGP$xb[1:SGGP$sizest[levellcv]],
                                                             SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
                                                             CorrMat=SGGP$CorrMat)
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
  }
  

  if (!SGGP$supplemented) {  
    
    # Return list with mean and var predictions
    if(is.vector(SGGP$pw)){
      GP = list("mean" = (SGGP$mu+Cp%*%SGGP$pw), "var"=SGGP$sigma2MAP[1]*ME_t)
    }else{ # y was a matrix, so PCA
      browser()
      if(length(SGGP$sigma2MAP)==1){
        GP = list("mean" = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+(Cp%*%SGGP$pw)%*%(SGGP$M)),
                  "var"=as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%(SGGP$sigma2MAP)%*%(SGGP$M))))
        
      }else{
        GP = list("mean" = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+(Cp%*%SGGP$pw)%*%(SGGP$M)),
                  "var"=as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%diag(SGGP$sigma2MAP)%*%(SGGP$M))))
      }
    }
  } else {
    stop("Fix pred for sggp$supp")
    Cps = matrix(1,dim(xp)[1],dim(SGGP$Xs)[1])
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      V = SGGP$CorrMat(xp[,dimlcv], SGGP$Xs[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
      Cps = Cps*V
    }
    
    yhatp = Cp%*%SGGP$pw_uppad + Cps%*%SGGP$supppw
    
    MSE_ps = list(matrix(0,dim(xp)[1],dim(SGGP$Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1)) 
    for (dimlcv in 1:SGGP$d) {
      for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
        MSE_ps[[(dimlcv)*SGGP$maxlevel+levellcv]] =(-SGGP_internal_postvarmatcalc(xp[,dimlcv],
                                                                                  SGGP$Xs[,dimlcv],
                                                                                  SGGP$xb[1:SGGP$sizest[levellcv]],
                                                                                  SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
                                                                                  CorrMat=SGGP$CorrMat))
      }
    }
    
    for (blocklcv in 1:SGGP$uoCOUNT) {
      ME_ps = matrix(1,nrow=dim(xp)[1],ncol=dim(SGGP$Xs)[1])
      for (dimlcv in 1:SGGP$d) {
        levelnow = SGGP$uo[blocklcv,dimlcv]
        ME_ps = ME_ps*MSE_ps[[(dimlcv)*SGGP$maxlevel+levelnow]]
      }
      Cps = Cps-SGGP$w[blocklcv]*(ME_ps)
    }
    ME_adj = rowSums((Cps%*%SGGP$Sti)*Cps)
    
    
    ME_t = ME_t-ME_adj
    # Return list with mean and var predictions
    if(is.vector(SGGP$pw)){
      GP = list("mean" = (SGGP$mu+ yhatp), "var"=SGGP$sigma2MAP[1]*ME_t)
    }else{
      if(length(SGGP$sigma2MAP)==1){
        GP = list("mean" = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ yhatp%*%(SGGP$M)),
                  "var"=as.vector(ME_t)%*%t(SGGP$leftover_variance+diag(t(SGGP$M)%*%(SGGP$sigma2MAP)%*%(SGGP$M))))
        
      }else{
        GP = list("mean" = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ yhatp%*%(SGGP$M)),
                  "var"=as.vector(ME_t)%*%t(SGGP$leftover_variance+diag(t(SGGP$M)%*%diag(SGGP$sigma2MAP)%*%(SGGP$M))))
      }
    }
    
    
  }
  
  return(GP)
}

#' Calculate MSE prediction along a single dimension
#'
#' @param xp Points at which to calculate MSE
#' @param xl Levels along dimension, vector???
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#'
#' @return MSE predictions
#' @export
#'
#' @examples
#' SGGP_internal_MSEpredcalc(c(.4,.52), c(0,.25,.5,.75,1), theta=c(.1,.2,.3),
#'              CorrMat=SGGP_internal_CorrMatCauchySQ)
SGGP_internal_MSEpredcalc <- function(xp,xl,theta,CorrMat) {
  S = CorrMat(xl, xl, theta)
  n = length(xl)
  cholS = chol(S)
  
  Cp = CorrMat(xp, xl, theta)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  
  MSE_val = 1 - rowSums(t(CiCp)*((Cp)))
  return(MSE_val)
}