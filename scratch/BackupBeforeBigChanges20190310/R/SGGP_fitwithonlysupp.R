SGGP_internal_fitwithonlysupp <- function(SGGP, Xs, Ys, use_PCA=TRUE,
                                          separateoutputparameterdimensions=TRUE,
                                          lower=rep(-1,SGGP$numpara*SGGP$d),upper=rep(1,SGGP$numpara*SGGP$d),
                                          theta0=rep(0,SGGP$numpara*SGGP$d),
                                          laplaceapprox=TRUE,
                                          use_progress_bar=TRUE,
                                          numPostSamples=NULL
                                          ) {
  
  SGGP$supplemented = TRUE
  SGGP$Xs = Xs
  SGGP$Ys = Ys
  SGGP$use_PCA <- use_PCA; rm(use_PCA)
  if (!is.null(numPostSamples)) {
    SGGP$numPostSamples <- numPostSamples
  }
  rm(numPostSamples)
  
  if(!is.matrix(Ys)){ # 1D output
    SGGP$mu = mean(Ys)
    # y = Y-SGGP$mu
    ys = Ys-SGGP$mu
  } else{ # Multiple outputs
    if (SGGP$use_PCA) {#browser()
      # Have to use Ys for data
      # if (dim(SGGP$Xs)[1] > 2*dim(Ys)[2]){
      SGGP$mu = colMeans(Ys)
      SGGP$st = (colMeans(Ys^2)- colMeans(Ys)^2)^(1/6) #somewhat arbitrary power, but seems to work. 1/2 is standard
      Ys_centered = (Ys - matrix(rep(SGGP$mu,each=dim(Ys)[1]), ncol=dim(Ys)[2], byrow=FALSE))%*%diag(1/SGGP$st)
      # SigV = 1/dim(Y)[1]*t(Ys_centered)%*%Ys_centered
      SigV = 1/dim(Ys)[1]*t(Ys_centered)%*%Ys_centered
      Eigen_result =eigen(SigV+10^(-12)*diag(length(SGGP$mu)))
      percent_explained = cumsum(sqrt(Eigen_result$values))/sum(sqrt(Eigen_result$values))
      # browser()
      num_PC = max(min(which(percent_explained>0.99999)),1)
      SGGP$M = t(Eigen_result$vectors[,1:num_PC])%*%diag(SGGP$st)
      ys = Ys_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
      
      Ys_recovered =   matrix(rep(SGGP$mu,each=dim(Xs)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ ys%*%(SGGP$M)
      SGGP$leftover_variance = colMeans((Ys-Ys_recovered)^2)
      
      # Y_centered = (Y - matrix(rep(SGGP$mu,each=dim(Y)[1]), ncol=dim(Y)[2], byrow=FALSE))%*%diag(1/SGGP$st)
      # y = Y_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
      # } else {
      #   SGGP$mu = colMeans(Y)
      #   SGGP$st = (colMeans(Y^2)- colMeans(Y)^2)^(1/6) #somewhat arbitrary power, but seems to work. 1/2 is standard
      #   Y_centered = (Y - matrix(rep(SGGP$mu,each=dim(Y)[1]), ncol=dim(Y)[2], byrow=FALSE))%*%diag(1/SGGP$st)
      #   SigV = 1/dim(Y)[1]*t(Y_centered)%*%Y_centered
      #   Eigen_result =eigen(SigV+10^(-12)*diag(length(SGGP$mu)))
      #   percent_explained = cumsum(sqrt(Eigen_result$values))/sum(sqrt(Eigen_result$values))
      #   num_PC = max(min(which(percent_explained>0.99999)),1)
      #   SGGP$M = t(Eigen_result$vectors[,1:num_PC])%*%diag(SGGP$st)
      #   y = Y_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
      #   
      #   Y_recovered =   matrix(rep(SGGP$mu,each=dim(SGGP$design)[1]), ncol=dim(SGGP$M)[2], byrow=FALSE)+ y%*%(SGGP$M)
      #   SGGP$leftover_variance = colMeans((Y-Y_recovered)^2)
      #   
      #   Ys_centered = (Ys - matrix(rep(SGGP$mu,each=dim(Ys)[1]), ncol=dim(Ys)[2], byrow=FALSE))%*%diag(1/SGGP$st)
      #   ys = Ys_centered%*%diag(1/SGGP$st)%*%t(SGGP$M)
      # }
    } else { # no PCA
      SGGP$mu = colMeans(Ys) # Could use Y, or colMeans(rbind(Y, Ys)), or make sure Ys is big enough for this
      # y <- sweep(Y, 2, SGGP$mu)
      ys <- sweep(Ys, 2, SGGP$mu)
    }
  }
  # SGGP$y = y
  SGGP$ys = ys
  
  
  
  
  
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
  if (is.matrix(Ys) && !SGGP$use_PCA) {SGGP$M <- diag(ncol(ys))} # Use identity transformation instead of PCA
  for (opdlcv in 1:nnn) { # output parameter dimension
    
    # y.thisloop <- if (nnn==1) {y} else {y[,opdlcv]} # All of y or single column
    # if (!is.null(Ys)) {
    ys.thisloop <- if (nnn==1) {ys} else {ys[,opdlcv]}
    # } # All of y or single column
    theta0.thisloop <- if (nnn==1) {theta0} else {theta0[,opdlcv]}
    
    opt.out = nlminb(
      theta0.thisloop,
      objective = SGGP_internal_neglogpost,
      gradient = SGGP_internal_gneglogpost,
      lower = lower, 
      upper = upper,
      y = NULL,
      HandlingSuppData="Only",
      Xs=Xs,
      ys=ys.thisloop,
      SGGP = SGGP,
      control = list(rel.tol = 1e-8,iter.max = 500)
    )
    # browser()
    # Set new theta
    thetaMAP <- opt.out$par
    
    # if (F) {
    #   # sigma2MAP <- SGGP_internal_calcsigma2anddsigma2(SGGP=SGGP, y=y.thisloop, theta=thetaMAP, return_lS=FALSE)$sigma2
    #   
    #   
    #   Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1])
    #   for (dimlcv in 1:SGGP$d) { # Loop over dimensions
    #     V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
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
    # pw <- SGGP_internal_calcpw(SGGP=SGGP, y.thisloop, theta=thetaMAP)
    
    
    totnumpara = length(thetaMAP)
    
    # H is the Hessian at thetaMAP with reverse transformation
    H = matrix(0,nrow=totnumpara,ncol=totnumpara)
    # Transform so instead of -1 to 1 it is -Inf to Inf. Mostly in -5 to 5 though.
    PSTn=  log((1+thetaMAP)/(1-thetaMAP))
    # Reverse transformation
    thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
    # Grad of reverse transformation function
    grad0 = SGGP_internal_gneglogpost(thetav,SGGP,y=NULL, Xs=Xs, ys=ys.thisloop, HandlingSuppData="Only")*(2*(exp(PSTn))/(exp(PSTn)+1)^2)
    for(c in 1:totnumpara){
      rsad = rep(0,totnumpara)
      rsad[c] =10^(-3)
      PSTn=  log((1+thetaMAP)/(1-thetaMAP)) + rsad
      thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
      H[c,] = (SGGP_internal_gneglogpost(thetav,SGGP,y=NULL, Xs=Xs, ys=ys.thisloop, HandlingSuppData="Only")*(2*(exp(PSTn))/(exp(PSTn)+1)^2)-grad0 )*10^(3)
    }
    Hmat = H/2+t(H)/2
    # print(Hmat)
    # print(sqrt(diag(solve(Hmat))))
    A = eigen(Hmat)
    cHa = (A$vectors)%*%diag(abs(A$values)^(-1/2))%*%t(A$vectors)
    #print( cHa%*%matrix(rnorm(100*length(SGGP$thetaMAP),0,1),nrow=length(SGGP$thetaMAP)))
    
    # Get posterior samples
    if(laplaceapprox){
      PST= log((1+thetaMAP)/(1-thetaMAP)) + cHa%*%matrix(rnorm(SGGP$numPostSamples*length(thetaMAP),0,1),nrow=length(thetaMAP))
      thetaPostSamples = (exp(PST)-1)/(exp(PST)+1)
    }else{ # MCMC Metropolis-Hastings
      # browser()
      U <- function(re){
        PSTn = log((1+thetaMAP)/(1-thetaMAP))+cHa%*%as.vector(re)
        thetav = (exp(PSTn)-1)/(exp(PSTn)+1)
        return(SGGP_internal_neglogpost(thetav,SGGP,y=NULL, Xs=Xs, ys=ys.thisloop, HandlingSuppData="Only"))
      }
      q = rep(0,totnumpara) # Initial point is 0, this is thetaMAP after transform
      Uo = U(q)
      if (is.infinite(Uo)) {warning("starting Uo is Inf, this is bad")}
      scalev = 0.5
      # This is just burn in period, find good scalev
      # MCMCtracker is for debugging
      if (use_progress_bar) {
        progress_bar <- progress::progress_bar$new(
          total=100*SGGP$numPostSamples*2,
          format = " MCMC [:bar] :percent eta: :eta"
        )
      }
      for(i in 1:(100*SGGP$numPostSamples)){
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
      Bs = matrix(0,nrow=totnumpara,ncol=SGGP$numPostSamples)
      
      for(i in 1:(100*SGGP$numPostSamples)){
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
    # if(SGGP$supplemented){
      # Cs = matrix(1,dim(SGGP$Xs)[1],SGGP$ss)
      # for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      #   V = SGGP$CorrMat(SGGP$Xs[,dimlcv], SGGP$xb, thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
      #   Cs = Cs*V[,SGGP$designindex[,dimlcv]]
      # }
      
      Sigma_t = matrix(1,dim(SGGP$Xs)[1],dim(SGGP$Xs)[1])
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(SGGP$Xs[,dimlcv], SGGP$Xs[,dimlcv], thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Sigma_t = Sigma_t*V
      }
      
      # MSE_s = list(matrix(0,dim(SGGP$Xs)[1],dim(SGGP$Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1)) 
      # for (dimlcv in 1:SGGP$d) {
      #   for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      #     MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]] =(-SGGP_internal_postvarmatcalc(SGGP$Xs[,dimlcv],SGGP$Xs[,dimlcv],
      #                                                                              SGGP$xb[1:SGGP$sizest[levellcv]],
      #                                                                              thetaMAP[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
      #                                                                              CorrMat=SGGP$CorrMat))
      #   }
      # }
      
      # for (blocklcv in 1:SGGP$uoCOUNT) {
      #   ME_s = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
      #   for (dimlcv in 1:SGGP$d) {
      #     levelnow = SGGP$uo[blocklcv,dimlcv]
      #     ME_s = ME_s*MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
      #   }
      #   Sigma_t = Sigma_t-SGGP$w[blocklcv]*(ME_s)
      # }
      # yhats = Cs%*%pw
      
      
      # Sti_resid = solve(Sigma_t,ys.thisloop-yhats)
      # Sti = solve(Sigma_t)
      Sti_chol <- chol(Sigma_t + diag(1e-10, nrow(Sigma_t), ncol(Sigma_t)))
      Sti <- chol2inv(Sti_chol)
      # sigma2MAP = (sigma2MAP*dim(SGGP$design)[1]+colSums((ys.thisloop-yhats)*Sti_resid))/(dim(SGGP$design)[1]+dim(Xs)[1])
      
      # pw_adj_y = t(Cs)%*%Sti_resid
      # pw_adj <- SGGP_internal_calcpw(SGGP=SGGP, y=pw_adj_y, theta=thetaMAP)
      
      # pw_uppadj = pw-pw_adj
      # supppw = solve(Sigma_t, ys)
      supppw <- Sti %*% ys.thisloop
      if (is.matrix(supppw) && ncol(supppw)==1) {supppw <- as.vector(supppw)}
      # print(expect_equal(supppw, solve(Sigma_t, ys)))
      sigma2MAP <- (t(ys.thisloop) %*% supppw) / nrow(Xs)
      if (is.matrix(sigma2MAP)) {sigma2MAP <- diag(sigma2MAP)}
    # }
    # browser()
    # Add all new variables to SGGP that are needed
    if (nnn==1) { # Only 1 output parameter dim, so just set them
      SGGP$thetaMAP <- thetaMAP
      SGGP$sigma2MAP <- sigma2MAP
      # SGGP$pw <- pw
      SGGP$thetaPostSamples <- thetaPostSamples
      if (SGGP$supplemented) {
        # SGGP$pw_uppadj <- pw_uppadj
        SGGP$supppw <- supppw
        SGGP$Sti = Sti
        SGGP$sigma2MAP <- sigma2MAP
      }
    } else { # More than 1 opd, so need to set as columns of matrix
      if (opdlcv==1) { # First time, initialize matrix/array for all
        
        SGGP$thetaMAP <- matrix(NaN, length(thetaMAP), nnn)
        if (length(sigma2MAP) != 1) {
          stop("ERROR HERE, sigma2map can be matrix??? It is always a 1x1 matrix from what I've seen before.")
        }
        SGGP$sigma2MAP <- numeric(nnn)
        # SGGP$pw <- matrix(NaN, length(pw), nnn) 
        # thetaPostSamples is matrix, so this is 3dim array below
        SGGP$thetaPostSamples <- array(data = NaN, dim=c(dim(thetaPostSamples), nnn))
      }
      SGGP$thetaMAP[,opdlcv] <- thetaMAP
      SGGP$sigma2MAP[opdlcv] <- sigma2MAP
      # SGGP$pw[,opdlcv] <- pw
      SGGP$thetaPostSamples[,,opdlcv] <- thetaPostSamples
      if (SGGP$supplemented) {
        if (opdlcv==1) { # First time initialize all
          
          # SGGP$pw_uppadj <- matrix(NaN, nrow(pw_uppadj), nnn)
          SGGP$supppw <- matrix(NaN, length(supppw), nnn) # Could be nrow supppw if 1 col matrix
          SGGP$Sti = array(NaN, dim=c(dim(Sti), nnn)) # Sti is matrix, so this is 3 dim array
        }
        # SGGP$pw_uppadj[,opdlcv] <- pw_uppadj
        SGGP$supppw[,opdlcv] <- supppw
        SGGP$Sti[,,opdlcv] = Sti
      }
    }
  }
  
  SGGP
}