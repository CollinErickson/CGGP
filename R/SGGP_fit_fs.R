#' Calculate negative log posterior
#'
#' @param theta Correlation parameters
#' @param SGGP SGGP object
#' @param y Measured values of SGGP$design
#'
#' @return Likelihood
#' @export
#' @useDynLib SGGP
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SGGP_internal_neglogpost(SG$thetaMAP, SG=SG, y=y)
SGGP_internal_neglogpost <- function(theta,SGGP,y) {
  # Return Inf if theta is too large
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
      neglogpost = 1/2*(length(y)*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
    }else{
      neglogpost = 1/2*(dim(y)[1]*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
    }
    
    return(neglogpost)
  }
  
}

#' Gradient of negative log likelihood posterior
#'
#' @param theta Log of correlation parameters
#' @param SGGP SGGP object
#' @param y SGGP$design measured values
#' @param return_lik If yes, it returns a list with lik and glik
#'
#' @return Vector for gradient of likelihood w.r.t. x (theta)
#' @export
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SGGP_internal_gneglogpost(SG$thetaMAP, SG=SG, y=y)
SGGP_internal_gneglogpost <- function(theta, SGGP, y, return_lik=FALSE) {
  
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
    neglogpost = 1/2*(length(y)*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
    gneglogpost = 0.500*(1/(1-theta)-1/(theta+1))+ dlDet+ length(y)*dsigma2_hat / sigma2_hat[1]
    gneglogpost =  gneglogpost/2
  }else{
    neglogpost = 1/2*(dim(y)[1]*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
    gneglogpost = 0.500*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
    for(i in 1:dim(y)[2]){
      gneglogpost = gneglogpost + dim(y)[1]*dsigma2_hat[,i] / sigma2_hat[i]
    }
    gneglogpost =  gneglogpost/2
  }
  
  if(return_lik){
    return(list(neglogpost=neglogpost,gneglogpost=gneglogpost))
  } else {
    return(gneglogpost)
  }
}

#' Calculate theta MLE given data
#'
#' @param SGGP Sparse grid objects
#' @param Y Output values calculated at SGGP$design
#' @param Xs Supplemental X matrix
#' @param Ys Supplemental Y values
#' @param theta0 Initial theta
#' @param laplaceapprox Should Laplace approximation be used?
#' @param lower Lower bound for parameter optimization
#' @param upper Upper bound for parameter optimization
#' @param Ynew Values of `SGGP$design_unevaluated`
# @param method Optimization method, must be "L-BFGS-B" when using lower and upper
# @param tol Relative tolerance for optimization. Can't use absolute tolerance
# since lik can be less than zero.
# @param return_optim If TRUE, return output from optim().
# If FALSE return updated SG.
#' 
#' @importFrom stats optim rnorm runif nlminb
#'
#' @return theta MLE
#' @export
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG=SG, Y=y)
SGGPfit <- function(SGGP, Y, Xs=NULL,Ys=NULL,
                   theta0 = rep(0,SGGP$numpara*SGGP$d),laplaceapprox = TRUE,
                   lower=rep(-1,SGGP$numpara*SGGP$d),upper=rep(1,SGGP$numpara*SGGP$d),
                   Ynew) {
  # If Ynew is given, it is only the points that were added last iteration. Append it to previous Y
  if (!missing(Ynew)) {
    if (!missing(Y)) {stop("Don't give both Y and Ynew, only one")}
    if (is.null(SGGP$Y)) {
      if (is.matrix(Ynew) && nrow(Ynew) != nrow(SGGP$design_unevaluated)) {stop("nrow(Ynew) doesn't match")}
      if (!is.matrix(Ynew) && length(Ynew) != nrow(SGGP$design_unevaluated)) {stop("length(Ynew) doesn't match")}
      Y <- Ynew
    } else if (is.matrix(SGGP$Y)) {
      if (!is.matrix(Ynew)) {stop("Ynew should be a matrix")}
      if (nrow(Ynew) != nrow(SGGP$design_unevaluated)) {stop("Ynew is wrong size")}
      Y <- rbind(SGGP$Y, Ynew)
    } else { # is numeric vector
      if (length(Ynew) != nrow(SGGP$design_unevaluated)) {stop("Ynew is wrong size")}
      Y <- c(SGGP$Y, Ynew)
    }
  }

  if ((is.matrix(Y) && nrow(Y) == nrow(SGGP$design)) || (length(Y) == nrow(SGGP$design))) {
    SGGP$design_unevaluated <- NULL
  } else {
    stop("SGGP$design and Y have different length")
  }
  
  #first do the pre-processing
  #for cleanness: Y is always the user input, y is after transformation
  SGGP$Y = Y
  if(is.null(Xs)){
    SGGP$supplemented = FALSE
    
    if(!is.matrix(Y)){
      SGGP$mu = mean(Y)
      y = Y-SGGP$mu
    }else{
      SGGP$mu = colMeans(Y)
      SGGP$st = (colMeans(Y^2)- colMeans(Y)^2)^(1/6) #somewhat arbitrary power, but seems to work. 1/2 is standard
      Y_centered = (Y - matrix(rep(SGGP$mu,each=dim(Y)[1]), ncol=dim(Y)[2], byrow=FALSE))%*%diag(1/SGGP$st)
      SigV = 1/dim(Y)[1]*t(Y_centered)%*%Y_centered
      Eigen_result =eigen(SigV+10^(-12)*diag(length(SGGP$mu)))
      # Had an error with small negative eigenvalues, so set those to zero force
      nonneg_Eigen_result_values <- pmax(Eigen_result$values, 0)
      percent_explained = cumsum(sqrt(nonneg_Eigen_result_values))/sum(sqrt(nonneg_Eigen_result_values))
      # percent_explained = cumsum(sqrt(Eigen_result$values))/sum(sqrt(Eigen_result$values))
      num_PC = max(min(which(percent_explained>0.99999)),1)
      SGGP$M = t(Eigen_result$vectors[,1:num_PC])%*%diag(SGGP$st)
      y = Y_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
      
      Y_recovered =   matrix(rep(SGGP$mu,each=dim(SGGP$design)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ y%*%(SGGP$M)
      SGGP$leftover_variance = colMeans((Y-Y_recovered)^2)
    }
    SGGP$y = y
    
  } else{
    SGGP$supplemented = TRUE
    SGGP$Xs = Xs
    SGGP$Ys = Ys
    
    if(!is.matrix(Y)){
      SGGP$mu = mean(Ys)
      y = Y-SGGP$mu
      ys = Ys-SGGP$mu
    } else{ 
      if (dim(SGGP$Xs)[1] > 2*dim(Ys)[2]){
        SGGP$mu = colMeans(Ys)
        SGGP$st = (colMeans(Ys^2)- colMeans(Ys)^2)^(1/6) #somewhat arbitrary power, but seems to work. 1/2 is standard
        Ys_centered = (Ys - matrix(rep(SGGP$mu,each=dim(Ys)[1]), ncol=dim(Ys)[2], byrow=FALSE))%*%diag(1/SGGP$st)
        SigV = 1/dim(Y)[1]*t(Ys_centered)%*%Ys_centered
        Eigen_result =eigen(SigV+10^(-12)*diag(length(SGGP$mu)))
        percent_explained = cumsum(sqrt(Eigen_result$values))/sum(sqrt(Eigen_result$values))
        num_PC = max(min(which(percent_explained>0.99999)),1)
        SGGP$M = t(Eigen_result$vectors[,1:num_PC])%*%diag(SGGP$st)
        ys = Ys_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
        
        Ys_recovered =   matrix(rep(SGGP$mu,each=dim(Xs)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ ys%*%(SGGP$M)
        SGGP$leftover_variance = colMeans((Ys-Ys_recovered)^2)
        
        Y_centered = (Y - matrix(rep(SGGP$mu,each=dim(Y)[1]), ncol=dim(Y)[2], byrow=FALSE))%*%diag(1/SGGP$st)
        y = Y_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
      } else {
        SGGP$mu = colMeans(Y)
        SGGP$st = (colMeans(Y^2)- colMeans(Y)^2)^(1/6) #somewhat arbitrary power, but seems to work. 1/2 is standard
        Y_centered = (Y - matrix(rep(SGGP$mu,each=dim(Y)[1]), ncol=dim(Y)[2], byrow=FALSE))%*%diag(1/SGGP$st)
        SigV = 1/dim(Y)[1]*t(Y_centered)%*%Y_centered
        Eigen_result =eigen(SigV+10^(-12)*diag(length(SGGP$mu)))
        percent_explained = cumsum(sqrt(Eigen_result$values))/sum(sqrt(Eigen_result$values))
        num_PC = max(min(which(percent_explained>0.99999)),1)
        SGGP$M = t(Eigen_result$vectors[,1:num_PC])%*%diag(SGGP$st)
        y = Y_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
        
        Y_recovered =   matrix(rep(SGGP$mu,each=dim(SGGP$design)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ y%*%(SGGP$M)
        SGGP$leftover_variance = colMeans((Y-Y_recovered)^2)
        
        Ys_centered = (Ys - matrix(rep(SGGP$mu,each=dim(Ys)[1]), ncol=dim(Ys)[2], byrow=FALSE))%*%diag(1/SGGP$st)
        ys = Ys_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
      }
    }
    SGGP$y = y
    SGGP$ys = ys
  }
  
  opt.out = nlminb(
    theta0,
    objective = SGGP_internal_neglogpost,
    gradient = SGGP_internal_gneglogpost,
    lower = lower, 
    upper = upper,
    y = y,
    SGGP = SGGP,
    #method = method, #"L-BFGS-B", #"BFGS", Only L-BFGS-B can use upper/lower
    #hessian = TRUE,
    control = list(rel.tol = 1e-8,iter.max = 500)#reltol=1e-4)#abstol = tol)
  )
  # for(i in 1:5){
  #   opt.out = nlminb(
  #     runif(SGGP$numpara*SGGP$d, -0.75,0.75),
  #     objective = SGGP_internal_neglogpost,
  #     gradient = SGGP_internal_gneglogpost,
  #     lower = 0.75*lower, 
  #     upper = 0.75*upper,
  #     y = y,
  #     SGGP = SGGP,
  #     #method = method, #"L-BFGS-B", #"BFGS", Only L-BFGS-B can use upper/lower
  #     #hessian = TRUE,
  #     control = list(rel.tol = 1e-8,iter.max = 500))#reltol=1e-4)#abstol = tol)
  #   }
  
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
  
  if(SGGP$supplemented){
    Cs = matrix(1,dim(SGGP$Xs)[1],SGGP$ss)
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      V = SGGP$CorrMat(SGGP$Xs[,dimlcv], SGGP$xb, SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
      Cs = Cs*V[,SGGP$designindex[,dimlcv]]
    }
    
    Sigma_t = matrix(1,dim(SGGP$Xs)[1],dim(SGGP$Xs)[1])
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      V = SGGP$CorrMat(SGGP$Xs[,dimlcv], SGGP$Xs[,dimlcv], SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
      Sigma_t = Sigma_t*V
    }
    
    MSE_s = list(matrix(0,dim(SGGP$Xs)[1],dim(SGGP$Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1)) 
    for (dimlcv in 1:SGGP$d) {
      for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
        MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]] =(-SGGP_internal_postvarmatcalc(SGGP$Xs[,dimlcv],SGGP$Xs[,dimlcv],
                                                                                  SGGP$xb[1:SGGP$sizest[levellcv]],SGGP$thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],CorrMat=SGGP$CorrMat))
      }
    }
    
    for (blocklcv in 1:SGGP$uoCOUNT) {
      ME_s = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
      for (dimlcv in 1:SGGP$d) {
        levelnow = SGGP$uo[blocklcv,dimlcv]
        ME_s = ME_s*MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
      }
      Sigma_t = Sigma_t-SGGP$w[blocklcv]*(ME_s)
    }
    
    yhats = Cs%*%SGGP$pw
    
    Sti_resid = solve(Sigma_t,ys-yhats)
    SGGP$Sti = solve(Sigma_t)
    SGGP$sigma2MAP = (SGGP$sigma2MAP*dim(SGGP$design)[1]+colSums((ys-yhats)*Sti_resid))/(dim(SGGP$design)[1]+dim(Xs)[1])
    
    pw_adj_y = t(Cs)%*%Sti_resid
    pw_adj <- SGGP_internal_calcpw(SGGP=SGGP, y=pw_adj_y, theta=SGGP$thetaMAP)
    
    SGGP$pw_uppadj = SGGP$pw-pw_adj
    SGGP$supppw = Sti_resid
  }
  return(SGGP)
}



#' ????????????
#'
#' @param x1 Points at which to calculate MSE
#' @param x2 Levels along dimension, vector???
#' @param xo No clue what this is
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#' for vector of 1D points.
#'
#' @return MSE predictions
#' @export
#'
#' @examples
#' SGGP_internal_postvarmatcalc(c(.4,.52), c(0,.25,.5,.75,1),
#'              xo=c(.11), theta=c(.1,.2,.3),
#'              CorrMat=SGGP_internal_CorrMatCauchySQ)
SGGP_internal_postvarmatcalc <- function(x1, x2, xo, theta, CorrMat) {
  S = CorrMat(xo, xo, theta)
  n = length(xo)
  cholS = chol(S)
  
  C1o = CorrMat(x1, xo, theta)
  CoinvC1o = backsolve(cholS,backsolve(cholS,t(C1o), transpose = TRUE))
  C2o = CorrMat(x2, xo, theta)
  Sigma_mat = - t(CoinvC1o)%*%t(C2o)  
  return(Sigma_mat)
}
