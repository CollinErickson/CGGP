

#' Update SGGP model given data
#'
#' @param SGGP Sparse grid objects
#' @param Y Output values calculated at SGGP$design
#' @param Xs Supplemental X matrix
#' @param Ys Supplemental Y values
#' @param theta0 Initial theta
#' @param lower Lower bound for parameter optimization
#' @param upper Upper bound for parameter optimization
#' @param Ynew Values of `SGGP$design_unevaluated`
# @param method Optimization method, must be "L-BFGS-B" when using lower and upper
# @param tol Relative tolerance for optimization. Can't use absolute tolerance
# since lik can be less than zero.
# @param return_optim If TRUE, return output from optim().
# If FALSE return updated SG.
#' @param separateoutputparameterdimensions If multiple output dimensions,
#' should separate parameters be fit to each dimension?
#' @param use_progress_bar If using MCMC sampling, should a progress bar be
#' displayed?
#' @param HandlingSuppData How should supplementary data be handled?
#' * Correct: full likelihood with grid and supplemental data
#' * Only: only use supplemental data
#' * Ignore: ignore supplemental data
#' @param corr Will update correlation function, if left missing it will be
#' same as last time.
#' @param ... Forces you to name arguments.
#' 
#' @importFrom stats optim rnorm runif nlminb
#'
#' @return Updated SGGP object fit to data given
#' @export
#' @family SGGP core functions
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG=SG, Y=y)
SGGPfit <- function(SGGP, Y, ..., Xs=NULL,Ys=NULL,
                    theta0 = SGGP$thetaMAP, #rep(0,SGGP$numpara*SGGP$d),
                    HandlingSuppData=SGGP$HandlingSuppData,
                    lower=rep(-1,SGGP$numpara*SGGP$d),upper=rep(1,SGGP$numpara*SGGP$d),
                    separateoutputparameterdimensions=is.matrix(SGGP$thetaMAP),
                    use_progress_bar=TRUE,
                    corr,
                    Ynew) {
  if (length(list(...))>0) {stop("Unnamed arguments given to SGGPcreate")}
  
  # If different correlation function is given, update it
  if (!missing(corr)) {
    message("Changing correlation function")
    SGGP <- SGGP_internal_set_corr(SGGP, corr)
  }
  
  # If Y or Ynew is matrix with 1 column, convert it to vector to avoid issues
  if (!missing(Y) && is.matrix(Y) && ncol(Y)==1) {
    Y <- c(Y)
  }
  if (!missing(Ynew) && is.matrix(Ynew) && ncol(Ynew)==1) {
    Ynew <- c(Ynew)
  }
  
  # If Ynew is given, it is only the points that were added last iteration. Append it to previous Y
  if (!missing(Ynew)) {
    if (!missing(Y)) {stop("Don't give both Y and Ynew, only one")}
    if (is.null(SGGP$Y) || length(SGGP$Y)==0) {
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
  if(is.null(Xs)){ # No supplemental data
    SGGP$supplemented = FALSE
    
    if(!is.matrix(Y)){
      SGGP$mu = mean(Y)
      y = Y-SGGP$mu
    }else{ # Y is matrix, PCA no longer an option
      SGGP$mu = colMeans(Y)
      y <- sweep(Y, 2, SGGP$mu)
      # Need to set SGGP$M somewhere so that it doesn't use transformation
    }
    SGGP$y = y
    
  } else{ # Has supplemental data, used for prediction but not for fitting params
    # stop("Not working for supp")
    SGGP$supplemented = TRUE
    SGGP$Xs = Xs
    SGGP$Ys = Ys
    
    if(!is.matrix(Y)){
      SGGP$mu = mean(Ys)
      y = Y-SGGP$mu
      ys = Ys-SGGP$mu
    } else{ # PCA no longer an option
      SGGP$mu = colMeans(Ys) # Could use Y, or colMeans(rbind(Y, Ys)), or make sure Ys is big enough for this
      y <- sweep(Y, 2, SGGP$mu)
      ys <- sweep(Ys, 2, SGGP$mu)
    }
    SGGP$y = y
    SGGP$ys = ys
  }
  
  # nopd is numberofoutputparameterdimensions
  nopd <- if (separateoutputparameterdimensions && is.matrix(y)) {
    ncol(y)
  } else {
    1
  }
  
  # Fix theta0
  if (nopd > 1) {
    if (is.vector(theta0)) {
      theta0 <- matrix(theta0, nrow=length(theta0), ncol=nopd, byrow=F)
    }
  }
  
  # Can get an error for theta0 if number of PCA dimensions has changed
  if (is.matrix(theta0) && (ncol(theta0) != nopd)) {
    if (ncol(theta0) > nopd) {
      theta0 <- theta0[,1:nopd]
    } else {
      theta0 <- cbind(theta0, matrix(0,nrow(theta0), nopd-ncol(theta0)))
    }
  }
  
  
  for (opdlcv in 1:nopd) { # output parameter dimension
    
    y.thisloop <- if (nopd==1) {y} else {y[,opdlcv]} # All of y or single column
    if (!is.null(Ys)) {ys.thisloop <- if (nopd==1) {ys} else {ys[,opdlcv]}} # All of y or single column
    else {ys.thisloop <- NULL}
    theta0.thisloop <- if (nopd==1) {theta0} else {theta0[,opdlcv]}
    
    if (is.null(SGGP$Xs)){ # No supp data, just optimize
      opt.out = nlminb(
        theta0.thisloop,
        objective = SGGP_internal_neglogpost,
        gradient = SGGP_internal_gneglogpost,
        lower = lower, 
        upper = upper,
        y = y.thisloop,
        SGGP = SGGP,
        ys = ys.thisloop,
        Xs = Xs,
        HandlingSuppData=HandlingSuppData,
        control = list(rel.tol = 1e-4,iter.max = 500)
      )
    } else { # W/ supp data, optimize on grid first, then with both
      # Only grid data b/c it's fast
      opt.out = nlminb(
        theta0.thisloop,
        objective = SGGP_internal_neglogpost,
        gradient = SGGP_internal_gneglogpost,
        lower = lower, 
        upper = upper,
        y = y.thisloop,
        SGGP = SGGP,
        HandlingSuppData="Ignore", # Never supp data here, so set to Ignore
        #  regardless of user setting
        control = list(rel.tol = 1e-2,iter.max = 500)
      )
      
      # Then use best point as initial point with supp data
      opt.out = nlminb(
        opt.out$par,
        objective = SGGP_internal_neglogpost,
        gradient = SGGP_internal_gneglogpost,
        lower = lower,
        upper = upper,
        y = y.thisloop,
        ys = ys.thisloop,
        Xs = Xs,
        SGGP = SGGP,
        HandlingSuppData = HandlingSuppData,
        control = list(rel.tol = 1e-4,iter.max = 500)
      )
      
    }
    
    
    
    # Set new theta
    thetaMAP <- opt.out$par
    sigma2MAP <- SGGP_internal_calcsigma2anddsigma2(SGGP=SGGP, y=y.thisloop, theta=thetaMAP, return_lS=FALSE)$sigma2
    # If one value, it gives it as matrix. Convert it to scalar
    if (length(sigma2MAP) == 1) {sigma2MAP <- sigma2MAP[1,1]}
    
    lik_stuff <- SGGP_internal_faststuff1(SGGP=SGGP, y.thisloop, theta=thetaMAP)
    cholSs = lik_stuff$cholS
    pw <- lik_stuff$pw
    totnumpara = length(thetaMAP)
    
    # H is the Hessian at thetaMAP with reverse transformation
    H = matrix(0,nrow=totnumpara,ncol=totnumpara)
    # Transform so instead of -1 to 1 it is -Inf to Inf. Mostly in -5 to 5 though.
    PSTn=  log((1+thetaMAP)/(1-thetaMAP))
    # Reverse transformation
    thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
    # Grad of reverse transformation function
    grad0 = SGGP_internal_gneglogpost(thetav,SGGP,y.thisloop,
                                      Xs=Xs, ys=ys.thisloop,
                                      HandlingSuppData=HandlingSuppData) *
      (2*(exp(PSTn))/(exp(PSTn)+1)^2)
    for(c in 1:totnumpara){
      rsad = rep(0,totnumpara)
      rsad[c] =10^(-3)
      PSTn=  log((1+thetaMAP)/(1-thetaMAP)) + rsad
      thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
      H[c,] = (SGGP_internal_gneglogpost(thetav,SGGP,y.thisloop,
                                         Xs=Xs, ys=ys.thisloop,
                                         HandlingSuppData=HandlingSuppData) *
                 (2*(exp(PSTn))/(exp(PSTn)+1)^2)-grad0 )*10^(3)
    }
    Hmat = H/2+t(H)/2
    A = eigen(Hmat)
    cHa = (A$vectors)%*%diag(abs(A$values)^(-1/2))%*%t(A$vectors)

    # Get posterior samples using Laplace approximation
    PST= log((1+thetaMAP)/(1-thetaMAP)) +
      cHa%*%matrix(rnorm(SGGP$numPostSamples*length(thetaMAP),0,1),
                   nrow=length(thetaMAP))
    thetaPostSamples = (exp(PST)-1)/(exp(PST)+1)
    
    
    if(SGGP$supplemented){
      Cs = matrix(1,dim(SGGP$Xs)[1],SGGP$ss)
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(SGGP$Xs[,dimlcv], SGGP$xb, thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Cs = Cs*V[,SGGP$designindex[,dimlcv]]
      }
      
      Sigma_t = matrix(1,dim(SGGP$Xs)[1],dim(SGGP$Xs)[1])
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(SGGP$Xs[,dimlcv], SGGP$Xs[,dimlcv], thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Sigma_t = Sigma_t*V
      }
      
      MSE_s = list(matrix(0,dim(SGGP$Xs)[1],dim(SGGP$Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1)) 
      for (dimlcv in 1:SGGP$d) {
        for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
          MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]] =
            (-SGGP_internal_postvarmatcalc(SGGP$Xs[,dimlcv],SGGP$Xs[,dimlcv],
                                           SGGP$xb[1:SGGP$sizest[levellcv]],
                                           thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
                                           CorrMat=SGGP$CorrMat))
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
      yhats = Cs%*%pw
      
      
      Sti_resid = solve(Sigma_t,ys.thisloop-yhats)
      Sti = solve(Sigma_t)
      sigma2MAP = (sigma2MAP*dim(SGGP$design)[1]+colSums((ys.thisloop-yhats)*Sti_resid))/(dim(SGGP$design)[1]+dim(Xs)[1])
      
      pw_adj_y = t(Cs)%*%Sti_resid
      pw_adj <- SGGP_internal_calcpw(SGGP=SGGP, y=pw_adj_y, theta=thetaMAP)
      
      pw_uppadj = pw-pw_adj
      supppw = Sti_resid
    }
    
    
    # Add all new variables to SGGP that are needed
    if (nopd==1) { # Only 1 output parameter dim, so just set them
      SGGP$thetaMAP <- thetaMAP
      SGGP$sigma2MAP <- sigma2MAP
      SGGP$pw <- pw
      SGGP$thetaPostSamples <- thetaPostSamples
      SGGP$cholSs = cholSs
      if (SGGP$supplemented) {
        SGGP$pw_uppadj <- pw_uppadj
        SGGP$supppw <- supppw
        SGGP$Sti = Sti
        SGGP$sigma2MAP <- sigma2MAP
      }
    } else { # More than 1 opd, so need to set as columns of matrix
      if (opdlcv==1) { # First time, initialize matrix/array for all
        
        SGGP$thetaMAP <- matrix(NaN, length(thetaMAP), nopd)
        if (length(sigma2MAP) != 1) {
          stop("ERROR HERE, sigma2map can be matrix??? It is always a 1x1 matrix from what I've seen before.")
        }
        SGGP$sigma2MAP <- numeric(nopd)
        SGGP$pw <- matrix(NaN, length(pw), nopd) 
        # thetaPostSamples is matrix, so this is 3dim array below
        SGGP$thetaPostSamples <- array(data = NaN, dim=c(dim(thetaPostSamples), nopd))
        # SGGP$cholSs <- array(data = NaN, dim=c(dim(cholSs), nopd))
        SGGP$cholSs <- vector("list", nopd) #array(data = NaN, dim=c(dim(cholSs), nopd))
      }
      SGGP$thetaMAP[,opdlcv] <- thetaMAP
      SGGP$sigma2MAP[opdlcv] <- sigma2MAP
      
      # browser()
      SGGP$pw[,opdlcv] <- pw
      SGGP$thetaPostSamples[,,opdlcv] <- thetaPostSamples
      # SGGP$cholSs[,,opdlcv] <- cholSs
      SGGP$cholSs[[opdlcv]] <- cholSs
      if (SGGP$supplemented) {
        if (opdlcv==1) { # First time initialize all
          
          SGGP$pw_uppadj <- matrix(NaN, nrow(pw_uppadj), nopd)
          SGGP$supppw <- matrix(NaN, nrow(supppw), nopd)
          SGGP$Sti = array(NaN, dim=c(dim(Sti), nopd)) # Sti is matrix, so this is 3 dim array
        }
        SGGP$pw_uppadj[,opdlcv] <- pw_uppadj
        SGGP$supppw[,opdlcv] <- supppw
        SGGP$Sti[,,opdlcv] = Sti
      }
    }
  }
  
  return(SGGP)
}



#' Calculate posterior variance
#'
#' @param x1 Points at which to calculate MSE
#' @param x2 Levels along dimension, vector???
#' @param xo No clue what this is
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#' for vector of 1D points.
#' @param ... Placeholder
#' @param returndPVMC Should dPVMC be returned?
#' @param returndiagonly Should only the diagonal be returned?
#'
#' @return Variance posterior
#' @export
#'
#' @examples
#' SGGP_internal_postvarmatcalc(c(.4,.52), c(0,.25,.5,.75,1),
#'              xo=c(.11), theta=c(.1,.2,.3),
#'              CorrMat=SGGP_internal_CorrMatCauchySQT)
SGGP_internal_postvarmatcalc <- function(x1, x2, xo, theta, CorrMat,...,returndPVMC = FALSE,returndiagonly=FALSE) {
  if(!returndiagonly){
    if(!returndPVMC){
      S = CorrMat(xo, xo, theta)
      n = length(xo)
      cholS = chol(S)
      
      C1o = CorrMat(x1, xo, theta)
      CoinvC1o = backsolve(cholS,backsolve(cholS,t(C1o), transpose = TRUE))
      C2o = CorrMat(x2, xo, theta)
      Sigma_mat = - t(CoinvC1o)%*%t(C2o)  
      return(Sigma_mat)
    }else{
      Sall = CorrMat(xo, xo, theta, return_dCdtheta = TRUE)
      S = Sall$C
      dS = Sall$dCdtheta
      n = length(xo)
      cholS = chol(S)
      
      
      Sall = CorrMat(x1, xo, theta, return_dCdtheta = TRUE)
      C1o = Sall$C
      dC1o = Sall$dCdtheta
      
      Sall = CorrMat(x2, xo, theta, return_dCdtheta = TRUE)
      C2o = Sall$C
      dC2o = Sall$dCdtheta
      
      CoinvC1o = backsolve(cholS,backsolve(cholS,t(C1o), transpose = TRUE))
      Sigma_mat = - t(CoinvC1o)%*%t(C2o)
      dSigma_mat = matrix(0,dim(Sigma_mat)[1],dim(Sigma_mat)[2]*length(theta))
      for(k in 1:length(theta)){
        CoinvC1oE = as.matrix(dS[,(n*(k-1)+1):(k*n)])%*%CoinvC1o
        dCoinvC1o = -backsolve(cholS,backsolve(cholS,CoinvC1oE, transpose = TRUE))
        dCoinvC1o = dCoinvC1o + backsolve(cholS,backsolve(cholS,t(dC1o[,(n*(k-1)+1):(k*n)]), transpose = TRUE))
        dSigma_mat[,((k-1)*dim(Sigma_mat)[2]+1):(k*dim(Sigma_mat)[2])] =
          -t(dCoinvC1o)%*%t(C2o)-t(CoinvC1o)%*%t(dC2o[,(n*(k-1)+1):(k*n)])
      }
      return(list("Sigma_mat"= Sigma_mat,"dSigma_mat" = dSigma_mat))
    }
  }else {
    if(!returndPVMC){
      S = CorrMat(xo, xo, theta)
      n = length(xo)
      cholS = chol(S)
      
      C1o = CorrMat(x1, xo, theta)
      CoinvC1o = backsolve(cholS,backsolve(cholS,t(C1o), transpose = TRUE))
      C2o = CorrMat(x2, xo, theta)
      Sigma_mat = -rowSums(t(CoinvC1o)*(C2o))
      return(Sigma_mat)
    }else{
      Sall = CorrMat(xo, xo, theta, return_dCdtheta = TRUE)
      S = Sall$C
      dS = Sall$dCdtheta
      n = length(xo)
      cholS = chol(S)
      
      Sall = CorrMat(x1, xo, theta, return_dCdtheta = TRUE)
      C1o = Sall$C
      dC1o = Sall$dCdtheta
      
      Sall = CorrMat(x2, xo, theta, return_dCdtheta = TRUE)
      C2o = Sall$C
      dC2o = Sall$dCdtheta
      
      CoinvC1o = backsolve(cholS,backsolve(cholS,t(C1o), transpose = TRUE))
      Sigma_mat = -rowSums(t(CoinvC1o)*(C2o))
      
      dSigma_mat = matrix(0,length(Sigma_mat),length(theta))
      for(k in 1:length(theta)){
        CoinvC1oE = as.matrix(dS[,(n*(k-1)+1):(k*n)])%*%CoinvC1o
        dCoinvC1o = -backsolve(cholS,backsolve(cholS,CoinvC1oE, transpose = TRUE))
        dCoinvC1o = dCoinvC1o + backsolve(cholS,backsolve(cholS,t(dC1o[,(n*(k-1)+1):(k*n)]), transpose = TRUE))
        dSigma_mat[,k] =-rowSums(t(dCoinvC1o)*(C2o))-rowSums(t(CoinvC1o)*(dC2o[,(n*(k-1)+1):(k*n)]))
      }
      return(list("Sigma_mat"= Sigma_mat,"dSigma_mat" = dSigma_mat))
    }
  }
  
}

