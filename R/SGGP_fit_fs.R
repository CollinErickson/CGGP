#' Calculate negative posterior
#'
#' @param theta Correlation parameters
#' @param SG SGGP object
#' @param y Measured values of SGGP$design
#' @param ... Don't use, just forces theta to be named
#'
#' @return Likelihood
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SGGP$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' lik(c(.1,.1,.1), SG=SG, y=y)
SGGP_internal_neglogpost <- function(theta,SGGP,y) {
  # Return Inf if theta is too large. Why????
  if (max(theta) >= 0.9999 || min(theta) <= -0.9999) {
    return(Inf)
  } else{
    calc_sigma2 <- SGGP_internal_calcsigma2(SGGP=SGGP, y=y, theta=theta, return_lS=TRUE)
    lS <- calc_sigma2$lS
    sigma2_hat = calc_sigma2$sigma2
    # Log determinant, keep a sum from smaller matrices
    lDet = 0
    # Calculate log det. See page 1586 of paper.
    # Loop over evaluated blocks
    for (blocklcv in 1:SGGP$uoCOUNT) {
      # Loop over dimensions
      for (dimlcv in 1:SGGP$d) {
        levelnow = SGGP$uo[blocklcv, dimlcv]
        # Add to log det when multiple points. It is zero when single point.
        if (levelnow > 1.5) {
          lDet = lDet + (lS[levelnow, dimlcv] - lS[levelnow - 1, dimlcv]) * (SGGP$gridsize[blocklcv]) /
            (SGGP$gridsizes[blocklcv, dimlcv])
        }
      }
    }
    
    
    if(!is.matrix(y)){
      neglogpost = 1/2*(length(y)*log(sigma2_hat[1])-0.501*sum(log(1-theta)+log(theta+1))+lDet)
    }else{
      neglogpost = 1/2*(dim(y)[1]*sum(log(c(sigma2_hat)))-0.501*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
    }
    
    return(neglogpost)
  }
  
}

#' Gradient of likelihood. Is it log likelihood?
#'
#' @param theta Log of correlation parameters
#' @param SG SGGP object
#' @param y SGGP$design measured values
#' @param return_lik If yes, it returns a list with lik and glik
#' @param ... Don't use, just forces theta to be named
#'
#' @return Vector for gradient of likelihood w.r.t. x (theta)
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SGGP$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' glik(c(.1,.1,.1), SG=SG, y=y)
SGGP_internal_gneglogpost <- function(theta, SGGP , y, return_lik=FALSE) {
  
  #print(theta)
  sigma2anddsigma2 <- SGGP_internal_calcsigma2anddsigma2(SGGP=SGGP, y=y, theta=theta, return_lS=TRUE)
  
  lS <- sigma2anddsigma2$lS
  dlS <-sigma2anddsigma2$dlS
  
  
  sigma2_hat = sigma2anddsigma2$sigma2
  dsigma2_hat = sigma2anddsigma2$dsigma2
  
  lDet = 0 # Not needed for glik, only for lik
  
  dlDet = rep(0, SGGP$numpara*SGGP$d) # Only needed for glik, not lik
  
  for (blocklcv in 1:SGGP$uoCOUNT) {
    nv = SGGP$gridsize[blocklcv]/SGGP$gridsizes[blocklcv,]
    uonow = SGGP$uo[blocklcv,]
    for (dimlcv in which(uonow>1.5)) {
      if (return_lik) {
        lDet = lDet + (lS[uonow[dimlcv], dimlcv] - lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
      }
      IS = (dimlcv-1)*SGGP$numpara+1:SGGP$numpara
      dlDet[IS] = dlDet[IS] + (dlS[uonow[dimlcv], IS] - dlS[uonow[dimlcv]-1, IS])*nv[dimlcv]
    }
  }
  
  
  if(!is.matrix(y)){
    neglogpost = 1/2*(length(y)*log(sigma2_hat[1])-0.501*sum(log(1-theta)+log(theta+1))+lDet)
    gneglogpost = 0.501*(1/(1-theta)-1/(theta+1))+ dlDet+ length(y)*dsigma2_hat / sigma2_hat[1]
    gneglogpost =  gneglogpost/2
  }else{
    neglogpost = 1/2*(dim(y)[1]*sum(log(c(sigma2_hat)))-0.501*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
    gneglogpost = 0.501*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
    for(i in 1:dim(y)[2]){
      gneglogpost = gneglogpost + dim(y)[1]*dsigma2_hat[,i] / sigma2_hat[i]
    }
    gneglogpost =  gneglogpost/2
  }
  
  if(return_lik){ nreturn(list(neglogpost=neglogpost,gneglogpost=gneglogpost))}else{return(gneglogpost)}
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
SGGPfit<- function(SGGP, Y, ..., 
                   lower = rep(-0.999, SGGP$d),
                   upper = rep(0.999, SGGP$d),
                   method = "L-BFGS-B", 
                   theta0 = rep(0,SGGP$numpara*SGGP$d),tol=1e-4,laplaceapprox = TRUE) {
  
  #first do the pre-processing
  #for cleanness: Y is always the user input, y is after transformation
  if(!is.matrix(Y)){
    SGGP$mu = mean(Y)
    y = Y-SGGP$mu
  }else{
    SGGP$mu = colMeans(Y)
    SGGP$st = (colMeans(Y^2)- colMeans(Y)^2)^(1/6) #somewhat arbitrary power, but seems to work. 1/2 is standard
    Y_centered = (Y - matrix(rep(SGGP$mu,each=dim(Y)[1]), ncol=dim(Y)[2], byrow=FALSE))%*%diag(1/SGGP$st)
    SigV = 1/dim(Y)[1]*t(Y_centered)%*%Y_centered
    Eigen_result =eigen(SigV+10^(-12)*diag(length(SGGP$mu)))
    percent_explained = cumsum(sqrt(Eigen_result$values))/sum(sqrt(Eigen_result$values))
    num_PC = max(min(which(percent_explained>0.9999)),1)
    y = Y_centered%*%(Eigen_result$vectors[,1:num_PC])
    SGGP$M = t(Eigen_result$vectors[,1:num_PC])%*%diag(SGGP$st)
  }
  
  opt.out = optim(
    theta0,
    fn = SGGP_internal_neglogpost,
    gr = SGGP_internal_gneglogpost,
    lower = lower, 
    upper = upper,
    y = y,
    SGGP = SGGP,
    method = method, #"L-BFGS-B", #"BFGS", Only L-BFGS-B can use upper/lower
    hessian = TRUE,
    control = list()#reltol=1e-4)#abstol = tol)
  )
  
  # Set new theta
  SGGP$thetaMAP <- opt.out$par
  SGGP$sigma2MAP <- SGGP_internal_calcsigma2anddsigma2(SGGP=SGGP, y=y, theta=SGGP$thetaMAP, return_lS=TRUE)$sigma2
  SGGP$pw <- SGGP_internal_calcpw(SGGP=SGGP, y, theta=SGGP$thetaMAP)
  totnumpara = length(SGGP$thetaMAP)
  
  H = matrix(0,nrow=totnumpara,ncol=totnumpara)
  PSTn=  log((1+SGGP$thetaMAP)/(1-SGGP$thetaMAP))
  thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
  grad0 = SGGP_internal_gneglogpost(thetav,SGGP,y)*(2*(exp(PSTn))/(exp(PSTn)+1)^2)
  for(c in 1:totnumpara){
    rsad = rep(0,totnumpara)
    rsad[c] =10^(-4)
    PSTn=  log((1+SGGP$thetaMAP)/(1-SGGP$thetaMAP)) + rsad
    thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
    H[c,] = (SGGP_internal_gneglogpost(thetav,SGGP,y)*(2*(exp(PSTn))/(exp(PSTn)+1)^2)-grad0 )*10^(4)
  }
  Hmat = H/2+t(H)/2
 # print(Hmat)
 # print(sqrt(diag(solve(Hmat))))
  A = eigen(Hmat)
  cHa = (A$vectors)%*%diag(abs(A$values)^(-1/2))%*%t(A$vectors)
  #print( cHa%*%matrix(rnorm(100*length(SGGP$thetaMAP),0,1),nrow=length(SGGP$thetaMAP)))
  if(laplaceapprox){
    PST= log((1+SGGP$thetaMAP)/(1-SGGP$thetaMAP)) + cHa%*%matrix(rnorm(100*length(SGGP$thetaMAP),0,1),nrow=length(SGGP$thetaMAP))
    SGGP$thetaPostSamples = (exp(PST)-1)/(exp(PST)+1)
  }else{
    U <- function(re){
      PSTn = log((1+SGGP$thetaMAP)/(1-SGGP$thetaMAP))+cHa%*%as.vector(re)
      thetav = (exp(PSTn)-1)/(exp(PSTn)+1)
      return(SGGP_internal_neglogpost(thetav,SGGP,y))
    }
    q = rep(0,totnumpara)
    Uo = U(q)
    scalev = 0.5
    for(i in 1:(100*SGGP$numPostSamples)){
      p = rnorm(length(q),0,1)*scalev
      qp = q + p
      
      Up = U(qp)
      if(runif(1) < exp(Uo-Up)){q=qp;Uo=Up;scalev=exp(log(scalev)+0.9/sqrt(i+4))}else{scalev=exp(log(scalev)-0.1/sqrt(i+4));scalev = max(scalev,1/sqrt(length(q)))}
      
    }
    
    Uo = U(q)
    Bs = matrix(0,ncol=totnumpara,nrow=SGGP$numPostSamples)
    for(i in 1:(100*SGGP$numPostSamples)){
      p = rnorm(length(q),0,1)*scalev
      qp = q + p
      
      Up = U(qp)
      if(runif(1) < exp(Uo-Up)){q=qp;Uo=Up;}
      if((i%%100)==0){Bs[,i/100]=q;}
    }
    PSTn = log((1+SGGP$thetaMAP)/(1-SGGP$thetaMAP))+cHa%*%Bs
    SGGP$thetaPostSamples = (exp(PSTn)-1)/(exp(PSTn)+1)
  }
  
  return(SGGP)
}