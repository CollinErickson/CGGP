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
#' * Mixture: sum of grid LLH and supplemental LLH, not statistically valid
#' * MarginalValidation: a validation shortcut
#' * FullValidation: a validation shortcut
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
    
    if (!(HandlingSuppData %in% c("Correct", "Only", "Ignore", "Mixture", "MarginalValidation", "FullValidation"))) {
      stop(paste("HandlingSuppData in SGGP_internal_neglogpost must be one of",
                 "Correct, Only, Ignore, Mixture, MarginalValidation, FullValidation"))
    }
    
    if(!is.null(ys) && is.null(y)){
      HandlingSuppData = "Only"
    }else if(is.null(ys) && !is.null(y)){
      HandlingSuppData = "Ignore"
    }else if(is.null(ys) && is.null(y)){
      stop(paste("You have given no y or ys to SGGP_internal_neglogpost"))
    }
    
    
    
    if(HandlingSuppData == "Only" || 
       HandlingSuppData == "Mixture" || 
       HandlingSuppData == "FullValidation" || 
       HandlingSuppData == "Correct"){
      Sigma_t = matrix(0,dim(Xs)[1],dim(Xs)[1])
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],returnlogs=TRUE)
        Sigma_t = Sigma_t+V
      }
      Sigma_t = exp(Sigma_t)
      
      Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
    }
    
    if(HandlingSuppData == "Only" || HandlingSuppData == "Mixture"){
      Sigma_chol = chol(Sigma_t)
      
      Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose = TRUE))
      sigma2_hat_supp = colSums(as.matrix(ys*Sti_resid))/dim(Xs)[1]
      lDet_supp = 2*sum(log(diag(Sigma_chol)))
    }
    
    if(HandlingSuppData == "Ignore" || HandlingSuppData == "Mixture"){
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
    
    if(HandlingSuppData == "FullValidation" ||
       HandlingSuppData == "Correct" ||
       HandlingSuppData == "MarginalValidation"){
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
    
    if(HandlingSuppData == "FullValidation" ||
       HandlingSuppData == "Correct" ){
      
      Sigma_t = matrix(0,dim(Xs)[1],dim(Xs)[1])
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],returnlogs=TRUE)
        Sigma_t = Sigma_t+V
      }
      Sigma_t = exp(Sigma_t)
      
      MSE_s = list(matrix(0,dim(Xs)[1],dim(Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1)) 
      Q  = max(SGGP$uo[1:SGGP$uoCOUNT,])
      for (dimlcv in 1:SGGP$d) {
        gg = (dimlcv-1)*Q
        for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
          INDSN = 1:SGGP$sizest[levellcv]
          INDSN = INDSN[sort(SGGP$xb[1:SGGP$sizest[levellcv]],index.return = TRUE)$ix]
          MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]] =(SGGP_internal_postvarmatcalcfaster(GGGG[[dimlcv]],
                                                                                        c(),
                                                                                        as.matrix(cholS[[gg+levellcv]]),
                                                                                        c(),
                                                                                        INDSN,
                                                                                        SGGP$numpara))
        }
      }
      
      # ME_s = matrix(0,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
      for (blocklcv in 1:SGGP$uoCOUNT) {
        if(abs(SGGP$w[blocklcv]) > 0.5){
        ME_s=matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
        for (dimlcv in 1:SGGP$d) {
          levelnow = SGGP$uo[blocklcv,dimlcv]
          ME_s = ME_s*MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
        }
        Sigma_t = Sigma_t-SGGP$w[blocklcv]*(ME_s)
        }
      }
      Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
      
      Sigma_chol = chol(Sigma_t)
      
      Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys-yhats,transpose = TRUE))
      sigma2_hat_supp = colSums((ys-yhats)*Sti_resid)/dim(Xs)[1]
      lDet_supp = 2*sum(log(diag(Sigma_chol)))
    }
    
    
    if(HandlingSuppData == "MarginalValidation"){
      Sigma_t = matrix(1,dim(Xs)[1],1)
      
      MSE_s = list(matrix(0,dim(Xs)[1],dim(Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1)) 
      for (dimlcv in 1:SGGP$d) {
        for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
          Q  = max(SGGP$uo[1:SGGP$uoCOUNT,])
          gg = (dimlcv-1)*Q
          INDSN = 1:SGGP$sizest[levellcv]
          INDSN = INDSN[sort(SGGP$xb[1:SGGP$sizest[levellcv]],index.return = TRUE)$ix]
          MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]] = (SGGP_internal_postvarmatcalcfaster(GGGG[[dimlcv]],
                                                                                         c(),
                                                                                         as.matrix(cholS[[gg+levellcv]]),
                                                                                         c(),
                                                                                         INDSN,
                                                                                         SGGP$numpara,
                                                                                         returndiag=TRUE))
        }
      }
      
      for (blocklcv in 1:SGGP$uoCOUNT) {
        if(abs(SGGP$w[blocklcv]) > 0.5){
        ME_s = matrix(1,nrow=dim(Xs)[1],1)
        for (dimlcv in 1:SGGP$d) {
          levelnow = SGGP$uo[blocklcv,dimlcv]
          ME_s = ME_s*MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
        }
        Sigma_t = Sigma_t-SGGP$w[blocklcv]*(ME_s)
        }
      }
      
      Sigma_t = (1-epsssV)*Sigma_t+epsssV
      lDet_supp = sum(log(Sigma_t))
      if(is.matrix(ys)){
        Sigma_t = t(matrix( rep( Sigma_t , dim(ys)[2] ) , ncol = ncol(t(Sigma_t)) , byrow = TRUE ))
      }
      sigma2_hat_supp = colSums((ys-yhats)^2/Sigma_t)/dim(Xs)[1]
    }
    
    neglogpost = 0
    
    if(HandlingSuppData == "Ignore" || HandlingSuppData == "Mixture"){
      if(!is.matrix(y)){
        neglogpost = neglogpost+1/2*((length(y))*log(sigma2_hat_grid[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet_grid)
      }else{
        neglogpost = neglogpost+1/2*((dim(y)[1])*sum(log(c(sigma2_hat_grid)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet_grid)
      }
    }
    
    if(HandlingSuppData == "Only" || HandlingSuppData == "Mixture"){
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
    
    if(HandlingSuppData =="FullValidation" || HandlingSuppData =="MarginalValidation"){  
      sigma2_hat = sigma2_hat_grid*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
      
      lDet = lDet_supp
      if(!is.matrix(y)){
        neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)+1/2*length(ys)*sigma2_hat_supp[1]/sigma2_hat[1]
      }else{
        neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet)+1/2*dim(ys)[1]*sum(sigma2_hat_supp/sigma2_hat)
      }
    }
    
    return(neglogpost)
  }
}
