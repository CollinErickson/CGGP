#' Calculate negative log posterior
#'
#' @param theta Correlation parameters
#' @param SGGP SGGP object
#' @param y Measured values of SGGP$design
#' @param ... Just a placeholder
#' @param Xs Supplementary input data
#' @param ys Supplementary output data
#' @param HandlingSuppData How should supplementary data be handled?
#' * Correct: full likelihood with grid and supplemental data
#' * Only: only use supplemental data
#' * Ignore: ignore supplemental data
#' @md
#'
#' @return Likelihood
#' @export
#' @useDynLib SGGP
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=20)
#' Y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG, Y)
#' SGGP_internal_neglogpost(SG$thetaMAP, SG=SG, y=SG$y)
SGGP_internal_neglogpost <- function(theta,SGGP,y,...,ys=NULL,Xs=NULL,HandlingSuppData = "Correct") {
  # Return Inf if theta is too large
  epsssV = 0
  
  if (max(theta) >= 1 || min(theta) <= -1) {
    return(Inf)
  } else{
    if (runif(1)<.001) message(paste("HandlingSuppData is:", HandlingSuppData))
    if (!(HandlingSuppData %in% c("Correct", "Only", "Ignore"))) {
      stop(paste("HandlingSuppData in SGGP_internal_neglogpost must be one of",
                 "Correct, Only, Ignore"))
    }
    
    if(!(is.null(ys) || length(ys)==0) && (is.null(y) || length(y)==0)){
      # Message user if actually changing it
      if (HandlingSuppData != "Only") {
        if (runif(1)<.01) message(paste("Changing HandlingSuppData to Only from", HandlingSuppData))
      }
      HandlingSuppData = "Only"
    }else if((is.null(ys) || length(ys)==0) && !(is.null(y) || length(y)==0)){
      # If making change, message user
      if (HandlingSuppData != "Ignore") {
        if (runif(1)<.01) message(paste("Changing HandlingSuppData to Ignore from", HandlingSuppData))
      }
      HandlingSuppData = "Ignore"
    }else if((is.null(ys) || length(ys)==0) && (is.null(y) || length(y)==0)){
      stop(paste("You have given no y or ys to SGGP_internal_neglogpost"))
    }
    
    
    
    if(HandlingSuppData == "Only" || 
       HandlingSuppData == "Correct"){
      Sigma_t = matrix(0,dim(Xs)[1],dim(Xs)[1])
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],returnlogs=TRUE)
        Sigma_t = Sigma_t+V
      }
      Sigma_t = exp(Sigma_t)
      
      Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
    }
    
    if(HandlingSuppData == "Only"){
      Sigma_chol = chol(Sigma_t)
      
      Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose = TRUE))
      sigma2_hat_supp = colSums(as.matrix(ys*Sti_resid))/dim(Xs)[1]
      lDet_supp = 2*sum(log(diag(Sigma_chol)))
    }
    
    if(HandlingSuppData == "Ignore"){
      sigma2anddsigma2 <- SGGP_internal_calcsigma2(SGGP=SGGP, y=y, theta=theta, return_lS=TRUE)
      lS <- sigma2anddsigma2$lS
      sigma2_hat_grid = sigma2anddsigma2$sigma2
      
      lDet_grid = 0
      for (blocklcv in 1:SGGP$uoCOUNT) {
        nv = SGGP$gridsize[blocklcv]/SGGP$gridsizes[blocklcv,]
        uonow = SGGP$uo[blocklcv,]
        for (dimlcv in which(uonow>1.5)) {
          lDet_grid = lDet_grid + (lS[uonow[dimlcv], dimlcv] - lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
        }
      }
    }
    
    if(HandlingSuppData == "Correct"){
      lik_stuff <- SGGP_internal_faststuff1(SGGP=SGGP, y=y, theta=theta)
      cholS = lik_stuff$cholS
      lS <- lik_stuff$lS
      sigma2_hat_grid = lik_stuff$sigma2
      pw = lik_stuff$pw
      
      lDet_grid = 0 # Not needed for glik, only for lik
      
      for (blocklcv in 1:SGGP$uoCOUNT) {
        nv = SGGP$gridsize[blocklcv]/SGGP$gridsizes[blocklcv,]
        uonow = SGGP$uo[blocklcv,]
        for (dimlcv in which(uonow>1.5)) {
          lDet_grid = lDet_grid + (lS[uonow[dimlcv], dimlcv] - lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
        }
      }
      
      #For these three, I need Cs,pw,allChols
      Cs = (matrix(0,dim(Xs)[1],SGGP$ss))
      GGGG = list(matrix(1,dim(Xs)[1],length(SGGP$xb)),SGGP$d)
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb[1:SGGP$sizest[max(SGGP$uo[,dimlcv])]],theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],returnlogs=TRUE)
        GGGG[[dimlcv]] = exp(V)
        Cs = Cs+V[,SGGP$designindex[,dimlcv]]
      }
      Cs = exp(Cs)
      yhats = Cs%*%pw
    }
    
    if(HandlingSuppData == "Correct" ){
      
      Sigma_t = matrix(0,dim(Xs)[1],dim(Xs)[1])
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],returnlogs=TRUE)
        Sigma_t = Sigma_t+V
      }
      Sigma_t = exp(Sigma_t)
      
      MSE_s = matrix(NaN,nrow=dim(Xs)[1]*dim(Xs)[1],ncol=(SGGP$d)*(SGGP$maxlevel))
      Q  = max(SGGP$uo[1:SGGP$uoCOUNT,])
      for (dimlcv in 1:SGGP$d) {
        gg = (dimlcv-1)*Q
        TT1 = GGGG[[dimlcv]]
        for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
          INDSN = 1:SGGP$sizest[levellcv]
          INDSN = INDSN[sort(SGGP$xb[1:SGGP$sizest[levellcv]],index.return = TRUE)$ix]
          REEALL =(SGGP_internal_postvarmatcalcfaster(TT1,
                                                                                        c(),
                                                                                        as.matrix(cholS[[gg+levellcv]]),
                                                                                        c(),
                                                                                        INDSN,
                                                                                        SGGP$numpara))
          MSE_s[,(dimlcv-1)*SGGP$maxlevel+levellcv]  =  as.vector(REEALL)
        }
      }
      
      Sigma_t2 = as.vector(Sigma_t)
      rcpp_fastmatclcr(SGGP$uo[1:SGGP$uoCOUNT,], SGGP$w[1:SGGP$uoCOUNT], MSE_s,Sigma_t2,SGGP$maxlevel)
      Sigma_t = matrix(Sigma_t2,nrow=dim(Xs)[1] , byrow = FALSE)
      
      Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
      
      Sigma_chol = chol(Sigma_t)
      
      Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys-yhats,transpose = TRUE))
      sigma2_hat_supp = colSums((ys-yhats)*Sti_resid)/dim(Xs)[1]
      lDet_supp = 2*sum(log(diag(Sigma_chol)))
    }
    
    neglogpost = 0
    
    if(HandlingSuppData == "Ignore"){
      if(!is.matrix(y)){
        neglogpost = neglogpost+1/2*((length(y))*log(sigma2_hat_grid[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet_grid)
      }else{
        neglogpost = neglogpost+1/2*((dim(y)[1])*sum(log(c(sigma2_hat_grid)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet_grid)
      }
    }
    
    if(HandlingSuppData == "Only"){
      if(!is.matrix(y)){
        neglogpost = neglogpost+1/2*((length(ys))*log(sigma2_hat_supp[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet_supp)
      }else{
        neglogpost = neglogpost+1/2*((dim(ys)[1])*sum(log(c(sigma2_hat_supp)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet_supp)
      }
    }
    
    
    if(HandlingSuppData =="Correct"){
      sigma2_hat = sigma2_hat_grid*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
      
      lDet = lDet_grid+lDet_supp
      if(!is.matrix(y)){
        neglogpost = 1/2*((length(y)+length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
      }else{
        neglogpost = 1/2*((dim(y)[1]+dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
      }
    }
    
    return(neglogpost)
  }
}
