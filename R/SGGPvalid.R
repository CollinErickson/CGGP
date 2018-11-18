#' Estimate correlation parameters using validation
#'
#' @inheritParams lik
#' @param xval Validation X matrix
#' @param yval Validation y vector
#'
#' @return logtheta estimates
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2+rnorm(1,0,.01)}
#' y <- apply(SG$design, 1, f1)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' validation(c(.1,.1,.1), SG=SG, y=y, xval=Xval, yval=Yval)
validation <- function(logtheta, SG, y,xval,yval) {

  # Return Inf if theta is too large. Why????
  if (max(logtheta) >= (4 - 10 ^ (-6))) {
    return(Inf)
  } else{
    
    yval = yval-mean(y)
    y = y - mean(y)
    
    # Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max value of all blocks
    # # Now going to store choleskys instead of inverses for stability
    # #CiS = list(matrix(1,1,1),Q*SG$d) # A list of matrices, Q for each dimension
    # CCS = list(matrix(1,1,1),Q*SG$d)
    # lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Save log determinant of matrices
    # # Loop over each dimension
    # for (lcv2 in 1:SG$d) {
    #   # Loop over each possible needed correlation matrix
    #   for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
    #     Xbrn = SG$xb[1:SG$sizest[lcv1]] # xb are the possible points
    #     Xbrn = Xbrn[order(Xbrn)] # Sort them low to high, is this necessary? Probably just needs to be consistent.
    #     S = SG$CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
    #     diag(S) = diag(S) + SG$nugget
    #     # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
    #     solvetry <- try({
    #       #CiS[[(lcv2-1)*Q+lcv1]] = solve(S)
    #       CCS[[(lcv2-1)*Q+lcv1]] = chol(S)
    #     })
    #     if (inherits(solvetry, "try-error")) {return(Inf)}
    #     #lS[lcv1, lcv2] = sum(log(eigen(S)$values))
    #     lS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
    #   }
    # }
    # 
    # # pw is the predictive weights, i.e., Sigma^{-1} * y
    # pw = rep(0, length(y)) # For each point
    # # Loop over blocks selected
    # for (lcv1 in 1:SG$uoCOUNT) {
    #   B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    #   for (e in SG$d:1) {
    #     if(SG$gridsizest[lcv1,e] > 1.5){
    #       B <- matrix(as.vector(B),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
    #       B <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B, transpose = TRUE))
    #       #B <-  solve(CS[[((e-1)*Q+SG$uo[lcv1,e])]],B)
    #       B <- t(B)
    #     }
    #     else{
    #       B = as.vector(B)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
    #     }
    #   }
    #   
    #   pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
    #     SG$w[lcv1] * B
    # }
    
    pw <- calculate_pw(SG=SG, y=y, logtheta=logtheta)
    
    Cp = matrix(1,dim(xval)[1],SG$ss)
    for (e in 1:SG$d) { # Loop over dimensions
      V = SG$CorrMat(xval[,e], SG$xb, logtheta=logtheta[e])
      Cp = Cp*V[,SG$designindex[,e]]
    }
    
    Ciy = Cp%*%pw
    yhatp = Ciy
    
    MSE_v = array(0, c(SG$d, 9,dim(xval)[1]))
    for (lcv1 in 1:SG$d) {
      MSE_v[lcv1, 1,] = SG$diag_corrMat(xval[,lcv1], logtheta=logtheta[lcv1], nugget=SG$nugget)
    }
    for (lcv1 in 1:SG$d) {
      for (lcv2 in 1:8) {
        MSE_v[lcv1, lcv2+1,] = abs(MSEpred_calc(xval[,lcv1],SG$xb[1:SG$sizest[lcv2]],
                                                logtheta=logtheta[lcv1], nugget=SG$nugget,
                                                CorrMat=SG$CorrMat,
                                                diag_corrMat=SG$diag_corrMat))
      }
    }
    
    ME_t = prod(MSE_v[,1,],1)
    for (lcv1 in 1:SG$uoCOUNT) {
      ME_v = rep(1,dim(xval)[1])
      for (e in 1:SG$d) {
        levelnow = SG$uo[lcv1,e]
        ME_v = ME_v*(MSE_v[e,levelnow,]-MSE_v[e,levelnow+1,]+10^(-12))
      }
      ME_t = ME_t-ME_v
    }
    
    #print(ME_t)
    sigma_hat = mean((yhatp-yval)^2/ME_t)
    
    # Where does sum(theta^2) come from? Looks like regularization? Or from coordinate transformation
    # This next line is really wrong? The paranthese closes off the return before including the lDet.
    pred_score = log(sigma_hat)+mean(log(ME_t))
  }
  
  pred_score
}


#' Calculate gradient of validation with respect to logtheta
#'
#' @inheritParams validation
#'
#' @return Vector, gradient of `validation`
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2+rnorm(1,0,.01)}
#' y <- apply(SG$design, 1, f1)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' gvalidation(c(.1,.1,.1), SG=SG, y=y, xval=Xval, yval=Yval)
gvalidation <- function(logtheta, SG, y,xval,yval) {

  yval = yval-mean(y)
  y = y - mean(y)
  
  # Calculate with function now instead.
  # Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max level of all blocks
  # # Now storing choleskys instead of inverses
  # CCS = list(matrix(1,1,1),Q*SG$d) # To store choleskys
  # dCS = list(matrix(1,1,1),Q*SG$d)
  # #CiS = list(matrix(1,1,1),Q*SG$d) # To store correlation matrices
  # #dCiS = list(matrix(1,1,1),Q*SG$d) # To store derivatives of corr mats
  # lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store log det S
  # dlS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store deriv of log det S
  # dlS2 = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # ???? added for chols
  # 
  # # Loop over each dimension
  # for (lcv2 in 1:SG$d) {
  #   # Loop over depth of each dim
  #   for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
  #     Xbrn = SG$xb[1:SG$sizest[lcv1]]
  #     Xbrn = Xbrn[order(Xbrn)]
  #     S = SG$CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
  #     diag(S) = diag(S) + SG$nugget
  #     dS = SG$dCorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
  #     
  #     CCS[[(lcv2-1)*Q+lcv1]] = chol(S)#-CiS[[(lcv2-1)*Q+lcv1]]  %*% 
  #     dCS[[(lcv2-1)*Q+lcv1]] = dS
  #     lS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
  #     V = solve(S,dS)
  #     dlS[lcv1, lcv2] = sum(diag(V))
  #   }
  # }
  # 
  # pw = rep(0, length(y)) # predictive weights
  # 
  # dpw = matrix(0, nrow = length(y), ncol = SG$d) # derivative of predictive weights
  # 
  # 
  # for (lcv1 in 1:SG$uoCOUNT) {
  #   
  #   B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
  #   for (e in SG$d:1) {
  #     if(SG$gridsizest[lcv1,e] > 1.5){
  #       B <- matrix(as.vector(B),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
  #       B <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B, transpose = TRUE))
  #       B <- t(B)
  #     }
  #     else{
  #       B = as.vector(B)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
  #     }
  #   }
  #   
  #   pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
  #     SG$w[lcv1] * B
  #   
  #   
  #   B3 = B
  #   for (e in  SG$d:1) {
  #     if(SG$gridsizest[lcv1,e] > 1.5){
  #       B3 <- matrix(as.vector(B3),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
  #     }
  #     
  #     B2 = B3
  #     
  #     if(SG$gridsizest[lcv1,e] > 1.5){
  #       B3 <- t(B3)
  #     } else{
  #       B3 = as.vector(B3)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
  #     }
  #     
  #     
  #     if(SG$gridsizest[lcv1,e] > 1.5){
  #       B2 <- -dCS[[((e-1)*Q+SG$uo[lcv1,e])]]%*%B2
  #       B2 <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B2, transpose = TRUE))
  #       B2 = t(B2)
  #     }else{
  #       B2 =-as.vector(dCS[[((e-1)*Q+SG$uo[lcv1,e])]])*as.vector(B2)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
  #     }
  #     
  #     if(e>1.5){
  #       for (e2 in  (e-1):1) {
  #         if(SG$gridsizest[lcv1,e2] > 1.5){
  #           B2 <- matrix(as.vector(B2),SG$gridsizest[lcv1,e2],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e2])
  #           B2 <- t(B2)
  #         }
  #         else{
  #           B2= as.vector(B2)
  #         }
  #       }
  #     }
  #     
  #     dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],e] = dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],e] +
  #       SG$w[lcv1] * B2
  #   }
  #   
  # }
  pw_dpw <- calculate_pw_and_dpw(SG=SG, y=y, logtheta=logtheta)  
  pw <- pw_dpw$pw
  dpw <- pw_dpw$dpw
  
  Cp = matrix(1,dim(xval)[1],SG$ss)
  for (e in 1:SG$d) { # Loop over dimensions
    #Cp = Cp*CorrMat(xval[,e], SG$design[,e], logtheta=logtheta[e]) # Multiply correlation from each dimension
    V = SG$CorrMat(xval[,e], SG$xb, logtheta=logtheta[e])
    Cp = Cp*V[,SG$designindex[,e]]
  }
  
  dCp = array(1,c(dim(xval)[1],SG$ss,SG$d))
  for (e in 1:SG$d) { # Loop over dimensions
    #Cp = Cp*CorrMat(xval[,e], SG$design[,e], logtheta=logtheta[e]) # Multiply correlation from each dimension
    V = SG$CorrMat(xval[,e], SG$xb, logtheta=logtheta[e])
    dV = SG$dCorrMat(xval[,e], SG$xb,logtheta=logtheta[e])
    dCp[,,e] = Cp/V[,SG$designindex[,e]]*dV[,SG$designindex[,e]]
  }
  
  Ciy = Cp%*%pw
  
  dCiy = array(1,c(dim(Ciy)[1],SG$d))
  for (e in 1:SG$d) { # Loop over dimensions
    dCiy[,e] = as.matrix(dCp[,,e])%*%pw+Cp%*%dpw[,e]
  }
  
  yhatp = Ciy
  dyhatp = dCiy
  
  MSE_v = array(0, c(SG$d, 9,dim(xval)[1]))
  dMSE_v = array(0, c(SG$d, 9,dim(xval)[1]))
  for (lcv1 in 1:SG$d) {
    MSE_v[lcv1, 1,] = SG$diag_corrMat(xval[,lcv1], logtheta=logtheta[lcv1], nugget=SG$nugget)
    dMSE_v[lcv1, 1,] = SG$ddiag_corrMat(xval[,lcv1], logtheta=logtheta[lcv1], nugget=SG$nugget)
  }
  for (lcv1 in 1:SG$d) {
    for (lcv2 in 1:8) {
      MSE_v[lcv1, lcv2+1,] = abs(MSEpred_calc(xval[,lcv1],SG$xb[1:SG$sizest[lcv2]],
                                              logtheta=logtheta[lcv1], nugget=SG$nugget,
                                              CorrMat=SG$CorrMat,
                                              diag_corrMat=SG$diag_corrMat))
      dMSE_v[lcv1, lcv2+1,] = abs(dMSEpred_calc(xval[,lcv1],SG$xb[1:SG$sizest[lcv2]],
                                                logtheta=logtheta[lcv1], nugget=SG$nugget,
                                                CorrMat=SG$CorrMat,
                                                dCorrMat=SG$dCorrMat,
                                                diag_corrMat=SG$diag_corrMat,
                                                ddiag_corrMat=SG$ddiag_corrMat))
    }
  }
  
  
  #print(MSE_v)
  # print(dMSE_v)
  
  ME_t = apply(t(MSE_v[,1,]),1,prod)
  dME_t = matrix(0,dim(dMSE_v)[3],SG$d)
  for (e in 1:SG$d) {
    dME_t[,e] = ME_t/MSE_v[e,1,]*dMSE_v[e,1,]
  }
  
  for (lcv1 in 1:SG$uoCOUNT) {
    ME_v = rep(1,dim(xval)[1])
    dME_v = matrix(1,dim(xval)[1],SG$d)
    for (e in 1:SG$d) {
      levelnow = SG$uo[lcv1,e]
      ME_v = ME_v*(MSE_v[e,levelnow,]-MSE_v[e,levelnow+1,]+10^(-12))
      #   for (e2 in 1:SG$d) {
      #     levelnow = SG$uo[lcv1,e]
      #     if(e!=e2){
      #       dME_v[,e2] = dME_v[,e2]*(MSE_v[e,levelnow,]-MSE_v[e,levelnow+1,])
      #     } else{
      #       dME_v[,e2] = dME_v[,e2]*(dMSE_v[e,levelnow,]-dMSE_v[e,levelnow+1,])
      #     }
      #}
    }
    for (e in 1:SG$d) {
      levelnow = SG$uo[lcv1,e]
      dME_v[,e] = ME_v/(MSE_v[e,levelnow,]-MSE_v[e,levelnow+1,]+10^(-12))*(dMSE_v[e,levelnow,]-dMSE_v[e,levelnow+1,])
    }
    
    
    ME_t = ME_t-ME_v
    dME_t = dME_t+dME_v
  }
  
  
  yhatmat = matrix(yhatp,nrow=length(yhatp),ncol=SG$d)
  ymat = matrix(as.matrix(yval),nrow=length(yval),ncol=SG$d)
  
  sigma_hat = mean((yhatp-yval)^2/ME_t)
  ME_tmat = matrix(ME_t,nrow=length(yval),ncol=SG$d)
  dsigma_hat =apply(-((yhatmat-ymat)^2/ME_tmat^2)*dME_t+2*(yhatmat-ymat)/(ME_tmat)*dyhatp,2,mean)
  
  # Where does sum(theta^2) come from? Looks like regularization? Or from coordinate transformation
  # This next line is really wrong? The paranthese closes off the return before including the lDet.
  # dpred_score = dsigma_hat/sigma_hat+
  dpred_score = dsigma_hat/sigma_hat+apply(dME_t/ME_tmat,2,mean)
  
  
  return(dpred_score)
  
}


#' Calculate correlation parameters using validation data.
#' 
#' Picks the logthetas that will minimize the score when predicting
#' on validation data.
#' The validation data should be representative of the entire region.
#'
#' @param SG SGGP object
#' @param y Output at SG$design
#' @param xval X matrix of validation points
#' @param yval Output of validation points
#' @param logtheta0 Initial point for optimization
#' @param tol Relative tolerance, not used
#'
#' @return Vector, optimal logtheta parameters
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2+rnorm(1,0,.01)}
#' y <- apply(SG$design, 1, f1)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' logthetaVALID(SG=SG, y=y, xval=Xval, yval=Yval)
logthetaVALID <- function(SG, y,xval,yval, logtheta0 = rep(0,SG$d),tol=1e-4, return_optim=FALSE) {
  opt.out = optim(
    logtheta0,
    fn = validation,
    gr = gvalidation,
    lower = rep(-2, SG$d),
    upper = rep(2.9, SG$d),
    y = y,
    yval = yval,
    xval = xval,
    SG = SG,
    method = "L-BFGS-B", #"BFGS",
    hessian = FALSE,
    control = list()#reltol=1e-4)#abstol = tol)
    # Is minimizing, default option of optim.
  )
  
  if (return_optim) {
    return(opt.out)
  }
  
  SG$logtheta <- opt.out$par
  SG
}
