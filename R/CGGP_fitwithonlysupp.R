CGGP_internal_fitwithonlysupp <- function(CGGP, Xs, Ys,
                                          separateoutputparameterdimensions=TRUE,
                                          lower=rep(-1,CGGP$numpara*CGGP$d),upper=rep(1,CGGP$numpara*CGGP$d),
                                          theta0=rep(0,CGGP$numpara*CGGP$d),
                                          laplaceapprox=TRUE,
                                          use_progress_bar=TRUE,
                                          numPostSamples=NULL
) {
  
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
  
  
  
  # This is all copied from fit
  
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
  
  
  
  # browser()
  if (is.matrix(Ys)) {CGGP$M <- diag(ncol(ys))} # Use identity transformation instead of PCA
  for (opdlcv in 1:nnn) { # output parameter dimension
    
    # y.thisloop <- if (nnn==1) {y} else {y[,opdlcv]} # All of y or single column
    # if (!is.null(Ys)) {
    ys.thisloop <- if (nnn==1) {ys} else {ys[,opdlcv]}
    # } # All of y or single column
    theta0.thisloop <- if (nnn==1) {theta0} else {theta0[,opdlcv]}
    
    opt.out = nlminb(
      theta0.thisloop,
      objective = CGGP_internal_neglogpost,
      gradient = CGGP_internal_gneglogpost,
      lower = lower, 
      upper = upper,
      y = NULL,
      HandlingSuppData="Only",
      Xs=Xs,
      ys=ys.thisloop,
      CGGP = CGGP,
      control = list(rel.tol = 1e-8,iter.max = 500)
    )
    # browser()
    # Set new theta
    thetaMAP <- opt.out$par
    
    # if (F) {
    #   # sigma2MAP <- CGGP_internal_calcsigma2anddsigma2(CGGP=CGGP, y=y.thisloop, theta=thetaMAP, return_lS=FALSE)$sigma2
    #   
    #   
    #   Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1])
    #   for (dimlcv in 1:CGGP$d) { # Loop over dimensions
    #     V = CGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
    #     Sigma_t = Sigma_t*V
    #   }
    #   Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
    #   
    #   Sigma_chol = chol(Sigma_t)
    #   
    #   Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys.thisloop,transpose = TRUE))
    #   sigma2_hat_supp = colSums(as.matrix(ys.thisloop*Sti_resid))/dim(Xs)[1]
    #   sigma2MAP <- sigma2_hat_supp
    #   
    #   supppw = Sti_resid
    # }    
    
    # If one value, it gives it as matrix. Convert it to scalar
    # if (length(sigma2MAP) == 1) {sigma2MAP <- sigma2MAP[1,1]}
    # pw <- CGGP_internal_calcpw(CGGP=CGGP, y.thisloop, theta=thetaMAP)
    
    
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
    # print(Hmat)
    # print(sqrt(diag(solve(Hmat))))
    A = eigen(Hmat)
    cHa = (A$vectors)%*%diag(abs(A$values)^(-1/2))%*%t(A$vectors)
    #print( cHa%*%matrix(rnorm(100*length(CGGP$thetaMAP),0,1),nrow=length(CGGP$thetaMAP)))
    
    # Get posterior samples
    if(laplaceapprox){
      PST= log((1+thetaMAP)/(1-thetaMAP)) + cHa%*%matrix(rnorm(CGGP$numPostSamples*length(thetaMAP),0,1),nrow=length(thetaMAP))
      thetaPostSamples = (exp(PST)-1)/(exp(PST)+1)
    }else{ # MCMC Metropolis-Hastings
      # browser()
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
        # print(Up)
        # MCMCtracker[i,13] <<- Up
        # MCMCtracker[i,1:12] <<- qp
        # MCMCtracker[i,14] <<- Uo
        
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
          # MCMCtracker[i,15] <<- 1
          q=qp;Uo=Up;scalev=exp(log(scalev)+0.9/sqrt(i+4))
          if (is.infinite(Up)) {warning("Up is Inf, this shouldn't happen")}
        }else{ # Reject sample, update scalev
          # MCMCtracker[i,15] <<- 0
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
    # browser("nothing done below here")
    
    # It is always supplemented
    # if(CGGP$supplemented){
    # Cs = matrix(1,dim(CGGP$Xs)[1],CGGP$ss)
    # for (dimlcv in 1:CGGP$d) { # Loop over dimensions
    #   V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$xb, thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
    #   Cs = Cs*V[,CGGP$designindex[,dimlcv]]
    # }
    
    Sigma_t = matrix(1,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1])
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$Xs[,dimlcv], thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
      Sigma_t = Sigma_t*V
    }
    
    # MSE_s = list(matrix(0,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1]),(CGGP$d+1)*(CGGP$maxlevel+1)) 
    # for (dimlcv in 1:CGGP$d) {
    #   for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
    #     MSE_s[[(dimlcv)*CGGP$maxlevel+levellcv]] =(-CGGP_internal_postvarmatcalc(CGGP$Xs[,dimlcv],CGGP$Xs[,dimlcv],
    #                                                                              CGGP$xb[1:CGGP$sizest[levellcv]],
    #                                                                              thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
    #                                                                              CorrMat=CGGP$CorrMat))
    #   }
    # }
    
    # for (blocklcv in 1:CGGP$uoCOUNT) {
    #   ME_s = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
    #   for (dimlcv in 1:CGGP$d) {
    #     levelnow = CGGP$uo[blocklcv,dimlcv]
    #     ME_s = ME_s*MSE_s[[(dimlcv)*CGGP$maxlevel+levelnow]]
    #   }
    #   Sigma_t = Sigma_t-CGGP$w[blocklcv]*(ME_s)
    # }
    # yhats = Cs%*%pw
    
    
    # Sti_resid = solve(Sigma_t,ys.thisloop-yhats)
    # Sti = solve(Sigma_t)
    Sti_chol <- chol(Sigma_t + diag(1e-10, nrow(Sigma_t), ncol(Sigma_t)))
    Sti <- chol2inv(Sti_chol)
    # sigma2MAP = (sigma2MAP*dim(CGGP$design)[1]+colSums((ys.thisloop-yhats)*Sti_resid))/(dim(CGGP$design)[1]+dim(Xs)[1])
    
    # pw_adj_y = t(Cs)%*%Sti_resid
    # pw_adj <- CGGP_internal_calcpw(CGGP=CGGP, y=pw_adj_y, theta=thetaMAP)
    
    # pw_uppadj = pw-pw_adj
    # supppw = solve(Sigma_t, ys)
    supppw <- Sti %*% ys.thisloop
    if (is.matrix(supppw) && ncol(supppw)==1) {supppw <- as.vector(supppw)}
    # print(expect_equal(supppw, solve(Sigma_t, ys)))
    sigma2MAP <- (t(ys.thisloop) %*% supppw) / nrow(Xs)
    if (is.matrix(sigma2MAP)) {sigma2MAP <- diag(sigma2MAP)}
    # }
    # browser()
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