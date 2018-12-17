#' Calculate negative posterior
#'
#' @param theta Correlation parameters
#' @param SG SGGP object
#' @param y Measured values of SGGP$design
#' @param ... Don't use, just forces theta to be named
#'
#' @return Likelihood
#' @export
#' @useDynLib SGGP
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SGGP$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' lik(c(.1,.1,.1), SG=SG, y=y)
SGGP_internal_predscore <- function(thetanew,SGGP,xp,Yp) {
  # Return Inf if theta is too large. Why????
  if (max(thetanew) >= 0.9999 || min(thetanew) <= -0.9999) {
    return(Inf)
  } else{
    calc_sigma2 <- SGGP_internal_calcsigma2(SGGP=SGGP, y=SGGP$y, thetanew)
    sigma2new = calc_sigma2$sigma2
    pwnew <- SGGP_internal_calcpw(SGGP=SGGP, SGGP$y,thetanew)
    
    # Cp is sigma(x_0) in paper, correlation vector between design points and xp
    Cp = matrix(1,dim(xp)[1],SGGP$ss)
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      V = SGGP$CorrMat(xp[,dimlcv], SGGP$xb, thetanew[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
      Cp = Cp*V[,SGGP$designindex[,dimlcv]]
    }
    MSE_v = array(0, c(SGGP$d, SGGP$maxgridsize,dim(xp)[1]))
    for (dimlcv in 1:SGGP$d) {
      MSE_v[dimlcv, 1,] = 1
    }
    for (dimlcv in 1:SGGP$d) {
      for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
        MSE_v[dimlcv, levellcv+1,] = SGGP_internal_MSEpredcalc(xp[,dimlcv],SGGP$xb[1:SGGP$sizest[levellcv]],thetanew[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],CorrMat=SGGP$CorrMat)
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
      Yhatp = (SGGP$mu+Cp%*%SGGP$pw)
      scale_variance =sigma2new*ME_t
    }else{
      if(length(SGGP$sigma2MAP)==1){
        Yhatp = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+(Cp%*%pwnew)%*%(SGGP$M))
        scale_variance =as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%(sigma2new)%*%(SGGP$M)))
        
      }else{
        Yhatp = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+(Cp%*%pwnew)%*%(SGGP$M))
        scale_variance =as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%diag(sigma2new)%*%(SGGP$M)))
      }
    }
    
    pred_score =  mean((Yhatp-Yp)^2/scale_variance+log(scale_variance))
    pred_score =  mean((Yhatp-Yp)^2)
    return(pred_score)
  }
}

#' Calculate negative posterior
#'
#' @param theta Correlation parameters
#' @param SG SGGP object
#' @param y Measured values of SGGP$design
#' @param ... Don't use, just forces theta to be named
#'
#' @return Likelihood
#' @export
#' @useDynLib SGGP
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SGGP$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' lik(c(.1,.1,.1), SG=SG, y=y)
SGGPpredsupplement <- function(xp,SGGP,Xs,Ys) {
  # Cp is sigma(x_0) in paper, correlation vector between design points and xp
  Cp = matrix(1,dim(xp)[1],SGGP$ss)
  for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    V = SGGP$CorrMat(Xp[,dimlcv], SGGP$xb, SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    Cp = Cp*V[,SGGP$designindex[,dimlcv]]
  }
  Cs = matrix(1,dim(Xs)[1],SGGP$ss)
  for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    V = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb, SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    Cs = Cs*V[,SGGP$designindex[,dimlcv]]
  }
  
  
  MSE_v = array(0, c(SGGP$d, SGGP$maxgridsize,dim(Xs)[1],dim(Xs)[1]))
  for (dimlcv in 1:SGGP$d) {
    MSE_v[dimlcv, 1,,] = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
  }
  for (dimlcv in 1:SGGP$d) {
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      MSE_v[dimlcv, levellcv+1,,] =SGGP_internal_postvarmatcalc(Xs[,dimlcv],Xs[,dimlcv],
                                                                SGGP$xb[1:SGGP$sizest[levellcv]],SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],CorrMat=SGGP$CorrMat)
    }
  }
  
  Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1])
  for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    Sigma_t = Sigma_t*V
  }
  
  for (blocklcv in 1:SGGP$uoCOUNT) {
    ME_v = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
    for (dimlcv in 1:SGGP$d) {
      levelnow = SGGP$uo[blocklcv,dimlcv]
      ME_v = ME_v*(MSE_v[dimlcv,levelnow,,]-MSE_v[dimlcv,levelnow+1,,])
    }
    Sigma_t = Sigma_t-ME_v
  }
  
  
  MSE_v = array(0, c(SGGP$d, SGGP$maxgridsize,dim(Xs)[1],dim(Xp)[1]))
  for (dimlcv in 1:SGGP$d) {
    MSE_v[dimlcv, 1,,] = SGGP$CorrMat(Xs[,dimlcv], Xp[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
  }
  for (dimlcv in 1:SGGP$d) {
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      MSE_v[dimlcv, levellcv+1,,] =SGGP_internal_postvarmatcalc(Xs[,dimlcv],Xp[,dimlcv],
                                                                SGGP$xb[1:SGGP$sizest[levellcv]],SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],CorrMat=SGGP$CorrMat)
    }
  }
  
  Sigma_pt = matrix(1,dim(Xs)[1],dim(Xp)[1])
  for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    V = SGGP$CorrMat(Xs[,dimlcv], Xp[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    Sigma_pt = Sigma_pt*V
  }
  C_pt = Sigma_pt
  
  for (blocklcv in 1:SGGP$uoCOUNT) {
    ME_v = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xp)[1])
    for (dimlcv in 1:SGGP$d) {
      levelnow = SGGP$uo[blocklcv,dimlcv]
      ME_v = ME_v*(MSE_v[dimlcv,levelnow,,]-MSE_v[dimlcv,levelnow+1,,])
    }
    Sigma_pt = Sigma_pt-ME_v
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
  
  
  
  Ys_centered = (Ys - matrix(rep(SGGP$mu,each=dim(Ys)[1]), ncol=dim(Ys)[2], byrow=FALSE))%*%diag(1/SGGP$st)
  ys = Ys_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
  
  yhatp = Cp%*%SGGP$pw
  yhats = Cs%*%SGGP$pw
  
  
  
  #here we go 
  yhatp = yhatp - Cp%*%pw_adj + t(C_pt)%*%solve(Sigma_t,ys-yhats)
  
  
  # Return list with mean and var predictions
  if(is.vector(SGGP$pw)){
    Yhatp = (SGGP$mu+yhatp)
    scale_variance =(SGGP$sigma2MAP)*ME_t
  }else{
    if(length(SGGP$sigma2MAP)==1){
      Yhatp = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+yhatp%*%(SGGP$M))
      scale_variance =as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%(SGGP$sigma2MAP)%*%(SGGP$M)))
      
    }else{
      Yhatp = ( matrix(rep(SGGP$mu,each=dim(xp)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+yhatp%*%(SGGP$M))
      scale_variance =as.vector(ME_t)%*%t(diag(t(SGGP$M)%*%diag(SGGP$sigma2MAP)%*%(SGGP$M)))
    }
  }
  if(is.vector(SGGP$pw)){
    GP = list("mean" = Yhatp, "var"=scale_variance)
  }else{
    if(length(SGGP$sigma2MAP)==1){
      GP =  list("mean" = Yhatp, "var"=scale_variance)
      
    }else{
      GP =  list("mean" = Yhatp, "var"=scale_variance)
    }
  }
  return(GP)
}

#' Calculate theta MLE given data
#'
#' @param SG Sparse grid objects
#' @param y Output values calculated at SGGP$design
#' @param theta0 Initial theta
#' @param tol Relative tolerance for optimization. Can't use absolute tolerance
#' since lik can be less than zero.
#' @param ... Don't use, just forces theta to be named
#' @param return_optim If TRUE, return output from optim().
#' @param lower Lower bound for parameter optimization
#' @param upper Upper bound for parameter optimization
#' @param method Optimization method, must be "L-BFGS-B" when using lower and upper
#' @param use_splitfngr Should give exact same results but with a slight speed up
#' If FALSE return updated SG.
#'
#' @return theta MLE
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SGGP$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' thetaMLE(SG=SG, y=y)
SGGPsupplement <- function(SGGP,Xs,Ys) {
  
  
  SGGP$supplemented = TRUE
  SGGP$Xs = Xs
  SGGP$Ys = Ys
  SGGP$ys = ys
  
  Ys_centered = (Ys - matrix(rep(SGGP$mu,each=dim(Ys)[1]), ncol=dim(Ys)[2], byrow=FALSE))%*%diag(1/SGGP$st)
  ys = Ys_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
  
  
  Cs = matrix(1,dim(Xs)[1],SGGP$ss)
  for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    V = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb, SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    Cs = Cs*V[,SGGP$designindex[,dimlcv]]
  }
  
  MSE_v = array(0, c(SGGP$d, SGGP$maxlevel,dim(Xs)[1],dim(Xs)[1]))
  for (dimlcv in 1:SGGP$d) {
    MSE_v[dimlcv, 1,,] = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
  }
  
  Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1])
  for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
    Sigma_t = Sigma_t*V
  }
  
  MSE_v = array(0, c(SGGP$d, SGGP$maxgridsize,dim(Xs)[1],dim(Xs)[1]))
  #for (dimlcv in 1:SGGP$d) {
  #  MSE_v[dimlcv, 1,,] = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
  #}
  for (dimlcv in 1:SGGP$d) {
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      MSE_v[dimlcv, levellcv+1,,] =SGGP_internal_postvarmatcalc(Xs[,dimlcv],Xs[,dimlcv],
                                                                SGGP$xb[1:SGGP$sizest[levellcv]],SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],CorrMat=SGGP$CorrMat)
    }
  }
  for (blocklcv in 1:SGGP$uoCOUNT) {
    ME_v = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
    for (dimlcv in 1:SGGP$d) {
      levelnow = SGGP$uo[blocklcv,dimlcv]
      ME_v = ME_v*(MSE_v[dimlcv,levelnow,,]-MSE_v[dimlcv,levelnow+1,,])
    }
    Sigma_t = Sigma_t-ME_v
  }
  
  yhats = Cs%*%SGGP$pw
  Sti_resid = solve(Sigma_t,ys-yhats)
  SGGP$Sti = solve(Sigma_t)
  SGGP$sigma2MAP = (SGGP$sigma2MAP*dim(SGGP$design)[1]+colSums((ys-yhats)*Sti_resid))/(dim(SGGP$design)[1]+dim(Xs)[1])
  
  pw_adj_y = t(Cs)%*%Sti_resid
  pw_adj <- SGGP_internal_calcpw(SGGP=SGGP, y=pw_adj_y, theta=SGGP$thetaMAP)
  
  SGGP$pw_uppadj = SGGP$pw-pw_adj
  SGGP$supppw = Sti_resid
  
  return(SGGP)
}


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
SGGP_internal_postvarmatcalc <- function(x1,x2,xo,theta,CorrMat) {
  S = CorrMat(xo, xo, theta)
  n = length(xo)
  cholS = chol(S)
  
  C1o = CorrMat(x1, xo, theta)
  CoinvC1o = backsolve(cholS,backsolve(cholS,t(C1o), transpose = TRUE))
  C2o = CorrMat(x2, xo, theta)
  C12 = CorrMat(x1, x2, theta)
  Sigma_mat = - t(CoinvC1o)%*%t(C2o)  
  return(Sigma_mat)
}
