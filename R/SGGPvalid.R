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
#' @param return_optim If TRUE, return output from optim().
#' If FALSE return updated SG.
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
  
  # Set new logtheta
  SG$logtheta <- opt.out$par
  
  # Save pw with SG
  SG$pw <- calculate_pw(SG=SG, y=y-mean(y), logtheta=SG$logtheta)
  
  SG
}
