#' Calculate likelihood
#'
#' @param theta Log of correlation parameters
#' @param SG SGGP object
#' @param y Measured values of SG$design
#' @param ... Don't use, just forces theta to be named
#'
#' @return Likelihood
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' lik(c(.1,.1,.1), SG=SG, y=y)
lik <- function(theta,SG, y) {
  
  # Return Inf if theta is too large. Why????
  if (max(theta) >= (4 - 10 ^ (-6))) {
    return(Inf)
  } else{
    
    calc_pw <- calculate_pw_C(SG=SG, y=y, theta=theta, return_lS=TRUE)
    pw <- calc_pw$pw
    lS <- calc_pw$lS
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
    return(log(c(sigma_hat))-10*sum(log(4-theta)+log(theta+4))/length(y)+1/ length(y) * lDet )
  }
  
}

#' Gradient of likelihood. Is it log likelihood?
#'
#' @param theta Log of correlation parameters
#' @param SG SGGP object
#' @param y SG$design measured values
#' @param return_lik If yes, it returns a list with lik and glik
#' @param ... Don't use, just forces theta to be named
#'
#' @return Vector for gradient of likelihood w.r.t. x (theta)
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' glik(c(.1,.1,.1), SG=SG, y=y)
glik <- function(theta, SG, y, return_lik=FALSE) {
  
  #print(theta)
  calc_pw_dpw <- calculate_pw_and_dpw_C(SG=SG, y=y, theta=theta, return_lS=TRUE)
  pw <- calc_pw_dpw$pw
  dpw <- calc_pw_dpw$dpw
  lS <- calc_pw_dpw$lS
  dlS <- calc_pw_dpw$dlS
  
  
  sigma_hat = t(y) %*% pw / length(y)
  dsigma_hat = c(t(y) %*% dpw) / length(y)
  
  lDet = 0 # Not needed for glik, only for lik
  
  dlDet = rep(0, SG$numpara*SG$d) # Only needed for glik, not lik
  
  for (lcv1 in 1:SG$uoCOUNT) {
    nv = SG$gridsize[lcv1]/SG$gridsizes[lcv1,]
    uonow = SG$uo[lcv1,]
    for (lcv2 in which(uonow>1.5)) {
      if (return_lik) {
        lDet = lDet + (lS[uonow[lcv2], lcv2] - lS[uonow[lcv2] - 1, lcv2])*nv[lcv2]
      }
      IS = (lcv2-1)*SG$numpara+1:SG$numpara
      dlDet[IS] = dlDet[IS] + (dlS[uonow[lcv2], IS] - dlS[uonow[lcv2]-1, IS])*nv[lcv2]
    }
  }
  
  
  ddL = dsigma_hat / sigma_hat[1]+10/(4-theta)/length(y)-10/(theta+4)/length(y)+ dlDet / length(y) 
  
  if (return_lik) {
    return(list(lik=log(c(sigma_hat))-10*sum(log(4-theta)+log(theta+4))/length(y)+ 1 / length(y) * lDet ,
                glik=ddL))
  }
  return(ddL)
}

#' Calculate theta MLE given data
#'
#' @param SG Sparse grid objects
#' @param y Output values calculated at SG$design
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
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' thetaMLE(SG=SG, y=y)
thetaMLE <- function(SG, y, ..., 
                     lower = rep(-3.9, SG$d),
                     upper = rep(3.9, SG$d),
                     method = "L-BFGS-B", 
                     theta0 = rep(0,SG$numpara*SG$d),tol=1e-4, return_optim=FALSE,
                     use_splitfngr=FALSE) {
  if (use_splitfngr) {
    opt.out = splitfngr::optim_share(
      par=theta0,
      fngr = function(par) glik(theta=par, SG=SG, y=y-mean(y), return_lik=T),
      lower = lower, #rep(-2, SG$d),
      upper = upper, #rep(3.9, SG$d),
      method = method, #"L-BFGS-B", #"BFGS", Only L-BFGS-B can use upper/lower
      hessian = TRUE,
      control = list()#reltol=1e-4)#abstol = tol)
    )
  } else {
    opt.out = optim(
      theta0,
      fn = lik,
      gr = glik,
      lower = lower, #rep(-2, SG$d),
      upper = upper, #rep(3.9, SG$d),
      y = y - mean(y),
      SG = SG,
      method = method, #"L-BFGS-B", #"BFGS", Only L-BFGS-B can use upper/lower
      hessian = TRUE,
      control = list()#reltol=1e-4)#abstol = tol)
    )
  }
  if (return_optim) {
    return(opt.out)
  }
  
  # Set new theta
  SG$theta <- opt.out$par
  SG$y = y
  
  # Save pw with SG
  SG$pw <- calculate_pw_C(SG=SG, y=y-mean(y), theta=SG$theta)
  
  SG
}


#' Likelihood for multivariate SGGP
#'
#' @param theta log of correlation parameters
#' @param SG SGGP object
#' @param yMV Matrix with output, number of columns is number of outputs
#'
#' @return Likelihood of theta
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y1 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' y2 <- apply(SG$design, 1, function(x){x[1]^1.3+.4*sin(6*x[2])+rnorm(1,0,.01)})
#' y <- cbind(y1, y2)
#' likMV(theta=c(.1,.2,.3), SG=SG, yMV=y)
likMV <- function(theta, SG, yMV) {
  p = dim(yMV)[2]
  
  logLikMV = 0
  for(i in 1:p){
    logLikMV = logLikMV+lik(theta, SG=SG, y=as.vector(yMV[,i])-mean(yMV[,i])) 
  }
  return(logLikMV)
}


#' Gradient of likelihood for multivariate SGGP
#'
#' @param theta log of correlation parameters
#' @param SG SGGP object
#' @param yMV Matrix with output, number of columns is number of outputs
#'
#' @return Vector, Gradient of likelihood of theta
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y1 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' y2 <- apply(SG$design, 1, function(x){x[1]^1.3+.4*sin(6*x[2])+rnorm(1,0,.01)})
#' y <- cbind(y1, y2)
#' glikMV(theta=c(.1,.2,.3), SG=SG, yMV=y)
glikMV <- function(theta, SG, yMV) {
  p = dim(yMV)[2]
  glogLikMV = rep(0,length(theta)) #p)
  for(i in 1:p){
    glogLikMV = glogLikMV+glik(theta, SG=SG, y=as.vector(yMV[,i])-mean(yMV[,i])) 
  }
  return(glogLikMV)
}

#' Find MLE of theta for multivariate output
#'
#' @param SG SG object
#' @param yMV Output matrix
#' @param ... Don't use
#' @param theta0 Initial values of theta for optimization
#' @param tol Optimization tolerance
#' @param return_optim If TRUE, return output from optim().
#' If FALSE return updated SG.
#'
#' @return Vector, theta MLE
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y1 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' y2 <- apply(SG$design, 1, function(x){x[1]^1.3+.4*sin(6*x[2])+rnorm(1,0,.01)})
#' y <- cbind(y1, y2)
#' thetaMLE(SG=SG, y=y1)
#' thetaMLE(SG=SG, y=y2)
#' thetaMLEMV(SG=SG, yMV=y)
thetaMLEMV <- function(SG, yMV, ..., theta0 = rep(0,SG$numpara*SG$d),tol=1e-4, return_optim=FALSE) {
  opt.out = optim(
    theta0,
    fn = likMV,
    gr = glikMV,
    lower = rep(-2, SG$d),
    upper = rep(3.9, SG$d),
    SG = SG,
    yMV = yMV,
    method = "L-BFGS-B", #"BFGS",
    hessian = FALSE,
    control = list()#reltol=1e-4)#abstol = tol)
  )
  
  if (return_optim) {
    return(opt.out)
  }
  
  # Set new theta
  SG$theta <- opt.out$par
  
  # Save pw with SG
  # Not as simple with multiple y outputs
  # SG$pw <- calculate_pw(SG=SG, y=yMV-mean(yMV), theta=SG$theta)
  
  SG
}
