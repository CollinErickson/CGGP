#' Fit CGGP to only supplemental data
#' 
#' If only supplemental data is given, and there is no grid data,
#' then this function is used.
#' It basically does the standard GP calculations since there is
#' nothing from the grid.
#'
#' @param CGGP CGGP object
#' @param Xs X supp matrix
#' @param Ys Y supp data
#' @param separateoutputparameterdimensions Should output dimensions be fit separately?
#' @param theta0 Initial theta0 for optimization
#' @param numPostSamples How many posterior samples should be calculated?
#'
#' @return CGGP object
## @export
#' @noRd
#' @examples
#' d <- 3
#' n <- 30
#' Xs <- matrix(runif(d*n), n, d)
#' Ys <- apply(Xs, 1, function(x){x[1]/(x[3]+.2)+exp(x[3])*cos(x[2]^2)})
#' cg <- CGGPcreate(d, Xs=Xs, Ys=Ys, batchsize=0)
CGGP_internal_fitwithonlysupp <- function(CGGP, Xs, Ys,
                                          separateoutputparameterdimensions=TRUE,
                                          theta0=rep(0,CGGP$numpara*CGGP$d),
                                          numPostSamples=NULL
) {
  # ===============
  #  Set Values
  # ===============
  CGGP$supplemented = TRUE
  CGGP$Xs = Xs
  CGGP$Ys = Ys
  if (!is.null(numPostSamples)) {
    CGGP$numPostSamples <- numPostSamples
  }
  rm(numPostSamples)
  
  if(!is.matrix(Ys)){ # 1D output
    CGGP$mu = mean(Ys)
    # y = Y-CGGP$mu
    ys = Ys-CGGP$mu
  } else{ # Multiple outputs
    CGGP$mu = colMeans(Ys) # Could use Y, or colMeans(rbind(Y, Ys)), or make sure Ys is big enough for this
    # y <- sweep(Y, 2, CGGP$mu)
    ys <- sweep(Ys, 2, CGGP$mu)
  }
  # CGGP$y = y
  CGGP$ys = ys
  
  
  
  
  
  # nnn is numberofoutputparameterdimensions
  nnn <- if (separateoutputparameterdimensions && is.matrix(ys)) {
    ncol(ys)
  } else {
    1
  }
  

  # Fix theta0
  if (nnn > 1) {
    if (is.vector(theta0)) {
      theta0 <- matrix(theta0, nrow=length(theta0), ncol=nnn, byrow=F)
    }
  }
  
  # Can get an error for theta0 if number of PCA dimensions has changed
  if (is.matrix(theta0) && (ncol(theta0) != nnn)) {
    if (ncol(theta0) > nnn) {
      theta0 <- theta0[,1:nnn]
    } else {
      theta0 <- cbind(theta0, matrix(0,nrow(theta0), nnn-ncol(theta0)))
    }
  }
  
  
  # ===========================================================
  #  Optimize parameters for each output parameter dimension
  # ===========================================================
  
  for (opdlcv in 1:nnn) { # output parameter dimension
    
    ys.thisloop <- if (nnn==1) {ys} else {ys[,opdlcv]}
    theta0.thisloop <- if (nnn==1) {theta0} else {theta0[,opdlcv]}
    
    opt.out = nlminb(
      theta0.thisloop,
      objective = CGGP_internal_neglogpost,
      gradient = CGGP_internal_gneglogpost,
      y = NULL,
      HandlingSuppData="Only",
      Xs=Xs,
      ys=ys.thisloop,
      CGGP = CGGP,
      control = list(rel.tol = 1e-8,iter.max = 500)
    )

    # Set new theta
    thetaMAP <- opt.out$par
    
    # ==========================
    #  Get posterior samples
    # ==========================
        
    totnumpara = length(thetaMAP)
    
    # H is the Hessian at thetaMAP with reverse transformation
    H = matrix(0,nrow=totnumpara,ncol=totnumpara)
    # Transform so instead of -1 to 1 it is -Inf to Inf. Mostly in -5 to 5 though.
    PSTn=  log((1+thetaMAP)/(1-thetaMAP))
    # Reverse transformation
    thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
    # Grad of reverse transformation function
    grad0 = CGGP_internal_gneglogpost(thetav,CGGP,y=NULL, Xs=Xs, ys=ys.thisloop, HandlingSuppData="Only")*(2*(exp(PSTn))/(exp(PSTn)+1)^2)
    for(c in 1:totnumpara){
      rsad = rep(0,totnumpara)
      rsad[c] =10^(-3)
      PSTn=  log((1+thetaMAP)/(1-thetaMAP)) + rsad
      thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
      H[c,] = (CGGP_internal_gneglogpost(thetav,CGGP,y=NULL, Xs=Xs, ys=ys.thisloop, HandlingSuppData="Only")*(2*(exp(PSTn))/(exp(PSTn)+1)^2)-grad0 )*10^(3)
    }
    Hmat = H/2+t(H)/2
    A = eigen(Hmat)
    cHa = (A$vectors)%*%diag(abs(A$values)^(-1/2))%*%t(A$vectors)
    
    # Get posterior samples
    if(laplaceapprox){
      PST= log((1+thetaMAP)/(1-thetaMAP)) + cHa%*%matrix(rnorm(CGGP$numPostSamples*length(thetaMAP),0,1),nrow=length(thetaMAP))
      thetaPostSamples = (exp(PST)-1)/(exp(PST)+1)
    }else{ # MCMC Metropolis-Hastings
      U <- function(re){
        PSTn = log((1+thetaMAP)/(1-thetaMAP))+cHa%*%as.vector(re)
        thetav = (exp(PSTn)-1)/(exp(PSTn)+1)
        return(CGGP_internal_neglogpost(thetav,CGGP,y=NULL, Xs=Xs, ys=ys.thisloop, HandlingSuppData="Only"))
      }
      q = rep(0,totnumpara) # Initial point is 0, this is thetaMAP after transform
      Uo = U(q)
      if (is.infinite(Uo)) {warning("starting Uo is Inf, this is bad")}
      scalev = 0.5
      # This is just burn in period, find good scalev
      # MCMCtracker is for debugging
      if (use_progress_bar) {
        progress_bar <- progress::progress_bar$new(
          total=100*CGGP$numPostSamples*2,
          format = " MCMC [:bar] :percent eta: :eta"
        )
      }
      for(i in 1:(100*CGGP$numPostSamples)){
        p = rnorm(length(q),0,1)*scalev
        qp = q + p
        
        Up = U(qp)
        
        # Sometimes Uo and Up equal Inf, so exp(Uo-Up)=NaN,
        #  and it gives an error: 
        #  Error in if (runif(1) < exp(Uo - Up)) { :
        #   missing value where TRUE/FALSE needed
        # Will just set Uo-Up to 0 when this happens
        exp_Uo_minus_Up <- exp(Uo-Up)
        if (is.nan(exp_Uo_minus_Up)) {
          warning("Uo and Up are both Inf, this shouldn't happen #82039")
          exp_Uo_minus_Up <- exp(0)
        }
        
        # if(runif(1) < exp(Uo-Up)){
        if(runif(1) < exp_Uo_minus_Up){ # accept the sample
          q=qp;Uo=Up;scalev=exp(log(scalev)+0.9/sqrt(i+4))
          if (is.infinite(Up)) {warning("Up is Inf, this shouldn't happen")}
        }else{ # Reject sample, update scalev
          scalev=exp(log(scalev)-0.1/sqrt(i+4))
          scalev = max(scalev,1/sqrt(length(q)))
        }
        if (use_progress_bar) {progress_bar$tick()}
      }
      
      Uo = U(q)
      Bs = matrix(0,nrow=totnumpara,ncol=CGGP$numPostSamples)
      
      for(i in 1:(100*CGGP$numPostSamples)){
        p = rnorm(length(q),0,1)*scalev
        qp = q + p # Next candidate is random normal step from q
        
        Up = U(qp)
        # Again need to avoid NaN problem
        exp_Uo_minus_Up <- exp(Uo-Up)
        if (is.nan(exp_Uo_minus_Up)) {
          warning("Uo and Up are both Inf, this shouldn't happen #42920")
          exp_Uo_minus_Up <- exp(0)
        }
        # if(runif(1) < exp(Uo-Up)){q=qp;Uo=Up;}
        if(runif(1) < exp_Uo_minus_Up){q=qp;Uo=Up;}
        # Only save every 100 samples for good mixing
        if((i%%100)==0){Bs[,i/100]=q;}
        if (use_progress_bar) {progress_bar$tick()}
      }
      if (use_progress_bar) {rm(progress_bar)}
      PSTn = log((1+thetaMAP)/(1-thetaMAP))+cHa%*%Bs
      thetaPostSamples = (exp(PSTn)-1)/(exp(PSTn)+1)
    }
    
    Sigma_t = matrix(1,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1])
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$Xs[,dimlcv], thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
      Sigma_t = Sigma_t*V
    }
    
    Sti_chol <- chol(Sigma_t + diag(1e-10, nrow(Sigma_t), ncol(Sigma_t)))
    Sti <- chol2inv(Sti_chol)
    
    supppw <- Sti %*% ys.thisloop
    if (is.matrix(supppw) && ncol(supppw)==1) {supppw <- as.vector(supppw)}
    
    sigma2MAP <- (t(ys.thisloop) %*% supppw) / nrow(Xs)
    if (is.matrix(sigma2MAP)) {sigma2MAP <- diag(sigma2MAP)}
    
    # Add all new variables to CGGP that are needed
    if (nnn==1) { # Only 1 output parameter dim, so just set them
      CGGP$thetaMAP <- thetaMAP
      CGGP$sigma2MAP <- sigma2MAP
      # CGGP$pw <- pw
      CGGP$thetaPostSamples <- thetaPostSamples
      if (CGGP$supplemented) {
        # CGGP$pw_uppadj <- pw_uppadj
        CGGP$supppw <- supppw
        CGGP$Sti = Sti
        CGGP$sigma2MAP <- sigma2MAP
      }
    } else { # More than 1 opd, so need to set as columns of matrix
      if (opdlcv==1) { # First time, initialize matrix/array for all
        
        CGGP$thetaMAP <- matrix(NaN, length(thetaMAP), nnn)
        if (length(sigma2MAP) != 1) {
          stop("ERROR HERE, sigma2map can be matrix??? It is always a 1x1 matrix from what I've seen before.")
        }
        CGGP$sigma2MAP <- numeric(nnn)
        # CGGP$pw <- matrix(NaN, length(pw), nnn) 
        # thetaPostSamples is matrix, so this is 3dim array below
        CGGP$thetaPostSamples <- array(data = NaN, dim=c(dim(thetaPostSamples), nnn))
      }
      CGGP$thetaMAP[,opdlcv] <- thetaMAP
      CGGP$sigma2MAP[opdlcv] <- sigma2MAP
      # CGGP$pw[,opdlcv] <- pw
      CGGP$thetaPostSamples[,,opdlcv] <- thetaPostSamples
      if (CGGP$supplemented) {
        if (opdlcv==1) { # First time initialize all
          
          # CGGP$pw_uppadj <- matrix(NaN, nrow(pw_uppadj), nnn)
          CGGP$supppw <- matrix(NaN, length(supppw), nnn) # Could be nrow supppw if 1 col matrix
          CGGP$Sti = array(NaN, dim=c(dim(Sti), nnn)) # Sti is matrix, so this is 3 dim array
        }
        # CGGP$pw_uppadj[,opdlcv] <- pw_uppadj
        CGGP$supppw[,opdlcv] <- supppw
        CGGP$Sti[,,opdlcv] = Sti
      }
    }
  }
  
  CGGP
}
