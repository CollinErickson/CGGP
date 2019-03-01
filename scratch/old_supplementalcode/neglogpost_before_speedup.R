# On 2/20 Matt redid SGGP_internal_neglogpost, SGGP_internal_gneglogpost, SGGP_internal_postvarmatcalc


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
#' SG <- SGGPcreate(d=3, batchsize=2000)
#' Y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG, Y)
#' SGGP_internal_neglogpost(SG$thetaMAP, SG=SG, y=SG$y)
SGGP_internal_neglogpostOLD <- function(theta,SGGP,y,...,ys=NULL,Xs=NULL,HandlingSuppData = "Correct") {
  # Return Inf if theta is too large
  epsssV = 0
  
  if (max(theta) >= 1 || min(theta) <= -1) {
    return(Inf)
  } else{
    
    if (!(HandlingSuppData %in% c("Correct", "Only", "Ignore", "Mixture", "MarginalValidation", "FullValidation"))) {
      stop(paste("HandlingSuppData in SGGP_internal_neglogpost must be one of",
                 "Correct, Only, Ignore, Mixture, MarginalValidation, FullValidation"))
    }
    
    if(!is.null(ys) && ( HandlingSuppData == "Only" || HandlingSuppData == "Mixture")){
      Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1])
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Sigma_t = Sigma_t*V
      }
      Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
      
      Sigma_chol = chol(Sigma_t)
      
      Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose = TRUE))
      sigma2_hat_supp = colSums(as.matrix(ys*Sti_resid))/dim(Xs)[1]
      lDet_supp = 2*sum(log(diag(Sigma_chol)))
      
      if( HandlingSuppData == "Only"){
        if(!is.matrix(y)){
          neglogpost = 1/2*((length(ys))*log(sigma2_hat_supp[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet_supp)
        }else{
          neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat_supp)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet_supp)
        }
      }
    }
    
    if(!is.null(y) && HandlingSuppData != "Only"){
      sigma2anddsigma2 <- SGGP_internal_calcsigma2anddsigma2(SGGP=SGGP, y=y, theta=theta, return_lS=TRUE)
      
      lS <- sigma2anddsigma2$lS
      dlS <-sigma2anddsigma2$dlS
      
      sigma2_hat_grid = sigma2anddsigma2$sigma2
      dsigma2_hat_grid = sigma2anddsigma2$dsigma2
      
      lDet_grid = 0 # Not needed for glik, only for lik
      
      for (blocklcv in 1:SGGP$uoCOUNT) {
        nv = SGGP$gridsize[blocklcv]/SGGP$gridsizes[blocklcv,]
        uonow = SGGP$uo[blocklcv,]
        for (dimlcv in which(uonow>1.5)) {
          lDet_grid = lDet_grid + (lS[uonow[dimlcv], dimlcv] - lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
        }
      }
    }
    
    
    if(HandlingSuppData == "Ignore" || HandlingSuppData == "Mixture" || is.null(ys)){
      if(!is.null(y)){
        if(!is.matrix(y)){
          neglogpost = 1/2*((length(y))*log(sigma2_hat_grid[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet_grid)
        }else{
          neglogpost = 1/2*((dim(y)[1])*sum(log(c(sigma2_hat_grid)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet_grid)
        }
      }else{
        neglogpost = 0
      }
      
      if ( HandlingSuppData == "Mixture" && !is.null(ys)){
        if(!is.matrix(y)){
          neglogpost = neglogpost+1/2*((length(ys))*log(sigma2_hat_supp[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet_supp)
        }else{
          neglogpost = neglogpost+1/2*((dim(ys)[1])*sum(log(c(sigma2_hat_supp)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet_supp)
        }
      }
    }
    
    
    if(!is.null(Xs) && (HandlingSuppData == "FullValidation" || HandlingSuppData == "Correct" || HandlingSuppData == "MarginalValidation" )){
      Cs = matrix(1,dim(Xs)[1],SGGP$ss)
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb,theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Cs = Cs*V[,SGGP$designindex[,dimlcv]]
      }
      
      pw <- SGGP_internal_calcpw(SGGP, y, theta)
      yhats = Cs%*%pw
      
      if(HandlingSuppData =="Correct" || HandlingSuppData =="FullValidation"){
        Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1])
        for (dimlcv in 1:SGGP$d) { # Loop over dimensions
          V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
          Sigma_t = Sigma_t*V
        }
        
        MSE_s = list(matrix(0,dim(Xs)[1],dim(Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1)) 
        for (dimlcv in 1:SGGP$d) {
          for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
            MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]] =(-SGGP_internal_postvarmatcalc(Xs[,dimlcv],SGGP$Xs[,dimlcv],
                                                                                     SGGP$xb[1:SGGP$sizest[levellcv]],
                                                                                     theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
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
        Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
        
        Sigma_chol = chol(Sigma_t)
        
        Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys-yhats,transpose = TRUE))
        sigma2_hat_supp = colSums((ys-yhats)*Sti_resid)/dim(Xs)[1]
        lDet_supp = 2*sum(log(diag(Sigma_chol)))
      }else{
        
        Sigma_t = matrix(1,dim(Xs)[1],1)
        
        MSE_s = list(matrix(0,dim(Xs)[1],1),(SGGP$d+1)*(SGGP$maxlevel+1)) 
        for (dimlcv in 1:SGGP$d) {
          for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
            MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]] =(-SGGP_internal_postvarmatcalc(Xs[,dimlcv],SGGP$Xs[,dimlcv],
                                                                                     SGGP$xb[1:SGGP$sizest[levellcv]],
                                                                                     theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
                                                                                     CorrMat=SGGP$CorrMat,returndiagonly=TRUE))
          }
        }
        
        for (blocklcv in 1:SGGP$uoCOUNT) {
          ME_s = matrix(1,nrow=dim(Xs)[1],1)
          for (dimlcv in 1:SGGP$d) {
            levelnow = SGGP$uo[blocklcv,dimlcv]
            ME_s = ME_s*MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
          }
          Sigma_t = Sigma_t-SGGP$w[blocklcv]*(ME_s)
        }
        Sigma_t = (1-epsssV)*Sigma_t+epsssV
        lDet_supp = -sum(log(Sigma_t))
        if(is.matrix(ys)){
          Sigma_t = t(matrix( rep( Sigma_t , dim(ys)[2] ) , ncol = ncol(t(Sigma_t)) , byrow = TRUE ))
        }
        sigma2_hat_supp = colSums((ys-yhats)^2/Sigma_t)/dim(Xs)[1]
      }
      
      
      
      if(HandlingSuppData =="Correct"){
        sigma2_hat = sigma2_hat_grid*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
        
        lDet = lDet_grid+lDet_supp
        if(!is.matrix(y)){
          neglogpost = 1/2*((length(y)+length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
        }else{
          neglogpost = 1/2*((dim(y)[1]+dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
        }
      }else if (HandlingSuppData =="FullValidation" || HandlingSuppData =="MarginalValidation"){  
        sigma2_hat = sigma2_hat_grid*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
        
        lDet = lDet_supp
        if(!is.matrix(y)){
          neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)+1/2*length(ys)*sigma2_hat_supp[1]/sigma2_hat[1]
        }else{
          neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet)+1/2*dim(ys)[1]*sum(sigma2_hat_supp/sigma2_hat)
        }
      }else{
        stop("ERROR")
      }
    }
  }
  
  return(neglogpost)
}


#' Gradient of negative log likelihood posterior
#'
#' @param theta Log of correlation parameters
#' @param SGGP SGGP object
#' @param y SGGP$design measured values
#' @param return_lik If yes, it returns a list with lik and glik
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
#'
#' @return Vector for gradient of likelihood w.r.t. x (theta)
#' @export
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=20)
#' Y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG, Y)
#' SGGP_internal_gneglogpost(SG$thetaMAP, SG=SG, y=SG$y)
SGGP_internal_gneglogpostOLD <- function(theta, SGGP, y,..., return_lik=FALSE,ys=NULL,Xs=NULL,HandlingSuppData = "Correct") {
  
  if (!(HandlingSuppData %in% c("Correct", "Only", "Ignore", "Mixture", "MarginalValidation", "FullValidation"))) {
    stop(paste("HandlingSuppData in SGGP_internal_neglogpost must be one of",
               "Correct, Only, Ignore, Mixture, MarginalValidation, FullValidation"))
  }
  
  epsssV = 0 # 10^(-14)
  
  if(is.null(y) || HandlingSuppData == "Only" || HandlingSuppData == "Mixture"){
    Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1])
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
      Sigma_t = Sigma_t*V
    }
    
    
    Cmat1 = t(matrix( rep( Sigma_t , SGGP$numpara ) , ncol = ncol(t(Sigma_t)) , byrow = TRUE ))
    dSigma_to = list(Sigma_t,SGGP$d) 
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      Vall = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta=TRUE)
      G = Vall$C
      Cmat2 =  t(matrix( rep( G , SGGP$numpara ) , ncol = ncol(t(G)) , byrow = TRUE ))
      dSigma_to[[dimlcv]] = Vall$dC*(Cmat1/Cmat2)
    }
    
    Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
    for (dimlcv in 1:SGGP$d) {
      dSigma_to[[dimlcv]] = (1-epsssV)*dSigma_to[[dimlcv]]
    }
    
    Sigma_chol = chol(Sigma_t)
    
    tempvec1 = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose = TRUE))
    sigma2_hat_new = colSums(as.matrix(ys*tempvec1))/dim(Xs)[1]
    if(is.matrix(ys)){
      dsigma2_hat_new = matrix(0,SGGP$d*SGGP$numpara,ncol=dim(ys)[2])
    }else{
      dsigma2_hat_new = rep(0,SGGP$d*SGGP$numpara)
    }
    dlDet_new=rep(0,SGGP$d*SGGP$numpara)
    for (dimlcv in 1:SGGP$d) {
      for(paralcv in 1:SGGP$numpara){
        dSigma_now = as.matrix((dSigma_to[[dimlcv]])[,((paralcv-1)*dim(Sigma_chol)[2]+1):(paralcv*dim(Sigma_chol)[2])])
        tempvec2= dSigma_now%*%tempvec1
        if(is.matrix(dsigma2_hat_new )){
          if(dim(dsigma2_hat_new)[1]>1.5){
            dsigma2_hat_new[(dimlcv-1)*SGGP$numpara+paralcv,] = -colSums(as.matrix(tempvec1*tempvec2))
          }else{
            dsigma2_hat_new[,(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(as.matrix(tempvec1*tempvec2))
          }
        }else{
          dsigma2_hat_new[(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(as.matrix(tempvec1*tempvec2))
        }
        dlDet_new[(dimlcv-1)*SGGP$numpara+paralcv] = sum(diag(backsolve(Sigma_chol,backsolve(Sigma_chol,dSigma_now,transpose = TRUE))))
      }
    }
    
    lDet_new = 2*sum(log(diag(Sigma_chol)))
    if( is.null(y) ||  HandlingSuppData == "Only" ){
      sigma2_hat = sigma2_hat_new
      dsigma2_hat = dsigma2_hat_new/dim(Xs)[1]
      dlDet = dlDet_new
      lDet = lDet_new
      if(!is.matrix(ys)){
        neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
        gneglogpost =0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*length(ys)*dsigma2_hat / sigma2_hat[1]
      }else{
        neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet)
        gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(ys)[2]*dlDet
        for(i in 1:dim(ys)[2]){
          gneglogpost = gneglogpost + (dim(ys)[1])*dsigma2_hat[,i] / sigma2_hat[i]
        }
        gneglogpost = gneglogpost/2
      }
    }
    
  }
  
  if(!is.null(y) && HandlingSuppData != "Only"){
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
    
    if(!is.null(Xs) && (HandlingSuppData == "FullValidation" || HandlingSuppData == "Correct" || HandlingSuppData == "MarginalValidation" )){
      Cs = matrix(1,dim(Xs)[1],SGGP$ss)
      dCs = array(1,dim=c(dim(Xs)[1],SGGP$numpara*SGGP$ss,dim(Xs)[2]))
      dsigma2_hat_old = dsigma2_hat
      sigma2_hat_old = sigma2_hat
      dlDet_old = dlDet
      lDet_old = lDet
      
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        Vall = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb,theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Cs = Cs*Vall[,SGGP$designindex[,dimlcv]]
      }
      
      pw <- SGGP_internal_calcpw(SGGP, y, theta)
      yhats = Cs%*%pw
      
      Cmat1 = matrix( rep(Cs , SGGP$numpara ) , nrow = nrow(Cs) , byrow = FALSE )
      Cmat2 = Cmat1
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        Vall = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb,theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta=TRUE)
        A = round(as.vector(outer(SGGP$designindex[,dimlcv],length(SGGP$xb)*(0:(SGGP$numpara-1)),FUN=function(x,y) x+y)))
        G = Vall$C[,SGGP$designindex[,dimlcv]]
        Cmat2 =  matrix( rep(G , SGGP$numpara ) , nrow = nrow(G) , byrow = FALSE )
        # print(dim(Cmat1))
        # print(dim(Cmat2))
        # print(dim(Vall$dC[,A]))
        dCs[,,dimlcv] = (Cmat1/Cmat2)*(Vall$dC[,A])
      }
      #print(Cs[1:10,1:10])
      
      if(HandlingSuppData =="Correct" || HandlingSuppData =="FullValidation"){
        
        Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1]) 
        dSigma_to = list(Sigma_t,SGGP$d) 
        for (dimlcv in 1:SGGP$d) { # Loop over dimensions
          Vall = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
          Sigma_t =Sigma_t*Vall
        }
        
        Cmat1 = matrix( rep(Sigma_t , SGGP$numpara ) , nrow = nrow(Sigma_t) , byrow = FALSE )
        dSigma_to = list(Sigma_t,SGGP$d) 
        for (dimlcv in 1:SGGP$d) { # Loop over dimensions
          Vall = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta=TRUE)
          G = Vall$C
          Cmat2 =  matrix( rep(G , SGGP$numpara ) , nrow = nrow(G) , byrow = FALSE )
          dSigma_to[[dimlcv]] = Vall$dC*(Cmat1/Cmat2)
        }
        
        MSE_s = list(matrix(0,dim(Xs)[1],dim(Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1))
        dMSE_s = list(matrix(0,dim(Xs)[1],dim(Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1))
        for (dimlcv in 1:SGGP$d) {
          for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
            REEALL =SGGP_internal_postvarmatcalc(SGGP$Xs[,dimlcv],SGGP$Xs[,dimlcv],SGGP$xb[1:SGGP$sizest[levellcv]],theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],CorrMat=SGGP$CorrMat,returndPVMC = TRUE)
            MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]]  = -REEALL$Sigma_mat
            dMSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]]   = -REEALL$dSigma_mat
          }
        }
        
        
        dsigma2_hat_part1 = 0*dsigma2_hat
        dsigma2_hat_part2 = 0*dsigma2_hat
        dsigma2_hat_part3 = 0*dsigma2_hat
        dlDet_new = 0*dlDet
        pw <- SGGP_internal_calcpw(SGGP, y, theta)
        for (blocklcv in 1:SGGP$uoCOUNT) {
          ME_s = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
          for (dimlcv in 1:SGGP$d) {
            levelnow = SGGP$uo[blocklcv,dimlcv]
            ME_s = ME_s*MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
          }
          Sigma_t = Sigma_t-SGGP$w[blocklcv]*(ME_s)
          
          
          for (dimlcv in 1:SGGP$d) {
            levelnow = SGGP$uo[blocklcv,dimlcv]
            ME_n = MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
            MEE_n = matrix( rep(ME_s/ME_n , SGGP$numpara ) , nrow = nrow(ME_s) , byrow = FALSE )
            dME_n = dMSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
            dSigma_to[[dimlcv]] = dSigma_to[[dimlcv]]-SGGP$w[blocklcv]*(MEE_n*dME_n)
          }
        }
        Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
        for (dimlcv in 1:SGGP$d) {
          dSigma_to[[dimlcv]] = (1-epsssV)*dSigma_to[[dimlcv]]
        }
        
        Sigma_chol = chol(Sigma_t)
        tempvec1= backsolve(Sigma_chol,backsolve(Sigma_chol,ys-yhats,transpose = TRUE))
        
        
        for (dimlcv in 1:SGGP$d) {
          for(paralcv in 1:SGGP$numpara){
            dCpn = as.matrix(dCs[,((paralcv-1)*dim(Cs)[2]+1):(paralcv*dim(Cs)[2]),dimlcv])
            if(is.matrix(dsigma2_hat)){
              if(dim(dsigma2_hat)[1]>1.5){
                dsigma2_hat_part2[(dimlcv-1)*SGGP$numpara+paralcv,] = -2*colSums((tempvec1)*(dCpn%*%pw))
              }else{
                dsigma2_hat_part2[,(dimlcv-1)*SGGP$numpara+paralcv] =  -2*colSums((tempvec1)*(dCpn%*%pw))
              }
            }else{
              dsigma2_hat_part2[(dimlcv-1)*SGGP$numpara+paralcv] = -2*colSums((tempvec1)*(dCpn%*%pw))
            }
            
            
            dSigma_now = as.matrix((dSigma_to[[dimlcv]])[,((paralcv-1)*dim(Sigma_chol)[2]+1):(paralcv*dim(Sigma_chol)[2])])
            tempvec2= dSigma_now%*%tempvec1
            if(is.matrix(dsigma2_hat)){
              if(dim(dsigma2_hat)[1]>1.5){
                dsigma2_hat_part1[(dimlcv-1)*SGGP$numpara+paralcv,] = -colSums(tempvec1*tempvec2)
              }else{
                dsigma2_hat_part1[,(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(tempvec1*tempvec2)
              }
            }else{
              dsigma2_hat_part1[(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(tempvec1*tempvec2)
            }
            
            
            dlDet_new[(dimlcv-1)*SGGP$numpara+paralcv] = sum(diag(backsolve(Sigma_chol,backsolve(Sigma_chol,dSigma_now,transpose = TRUE))))
          }
        }
        
        if(is.vector(y)){
          temp4 = as.vector(t(Cs)%*%tempvec1)
        }else{
          temp4 = t(Cs)%*%tempvec1
        }
        ytemp = y
        dsigma2_hat_part3 =  -2*(SGGP_internal_calcusedforsupp(SGGP,ytemp,temp4, theta)$dsigma2)
        
        sigma2_hat_new = colSums((ys-yhats)*tempvec1)/dim(Xs)[1]
        lDet_new = 2*sum(log(diag(Sigma_chol)))
        
        dsigma2_hat_new = (dsigma2_hat_part1+dsigma2_hat_part2+dsigma2_hat_part3)/dim(Xs)[1]
      }else{
        Sigma_t = matrix(1,dim(Xs)[1],1)
        
        Cmat1 = 0*matrix( rep( (Sigma_t),SGGP$numpara) , nrow = nrow(Cs) , byrow = FALSE )
        dSigma_to = list(Cmat1,SGGP$d) 
        for (dimlcv in 1:SGGP$d) {
          dSigma_to[[dimlcv]] = 0*Cmat1
        }
        
        MSE_s = list(matrix(0,dim(Xs)[1],1),(SGGP$d+1)*(SGGP$maxlevel+1))
        dMSE_s = list(matrix(0,dim(Xs)[1],1),(SGGP$d+1)*(SGGP$maxlevel+1))
        for (dimlcv in 1:SGGP$d) {
          for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
            REEALL =SGGP_internal_postvarmatcalc(SGGP$Xs[,dimlcv],SGGP$Xs[,dimlcv],SGGP$xb[1:SGGP$sizest[levellcv]],theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],CorrMat=SGGP$CorrMat,returndPVMC = TRUE,returndiagonly=TRUE)
            MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]]  = -(REEALL$Sigma_mat)
            dMSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]]   = -(REEALL$dSigma_mat)
          }
        }
        
        
        dsigma2_hat_part1 = 0*dsigma2_hat
        dsigma2_hat_part2 = 0*dsigma2_hat
        dsigma2_hat_part3 = 0*dsigma2_hat
        dlDet_new = 0*dlDet
        for (blocklcv in 1:SGGP$uoCOUNT) {
          ME_s = matrix(1,nrow=dim(Xs)[1],ncol=1)
          for (dimlcv in 1:SGGP$d) {
            levelnow = SGGP$uo[blocklcv,dimlcv]
            ME_s = ME_s*MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
          }
          Sigma_t = Sigma_t-SGGP$w[blocklcv]*(ME_s)
          
          ME_so1 = t(matrix( rep( (ME_s) , SGGP$numpara ) , ncol = ncol(t(ME_s)) , byrow = TRUE ))
          ME_so2 = ME_so1
          for (dimlcv in 1:SGGP$d) {
            levelnow = SGGP$uo[blocklcv,dimlcv]
            ME_n = (MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]])
            ME_so2 =  t(matrix( rep( (ME_n) , SGGP$numpara ) , ncol = ncol(t(ME_n)) , byrow = TRUE ))
            dME_n = dMSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
            dSigma_to[[dimlcv]] = dSigma_to[[dimlcv]]-SGGP$w[blocklcv]*(ME_so1*dME_n/ME_so2)
          }
        }
        
        Sigma_t = (1-epsssV)*Sigma_t+epsssV
        for (dimlcv in 1:SGGP$d) {
          dSigma_to[[dimlcv]] = (1-epsssV)*dSigma_to[[dimlcv]]
        }
        
        for (dimlcv in 1:SGGP$d) {
          for(paralcv in 1:SGGP$numpara){
            dSigma_now = as.matrix((dSigma_to[[dimlcv]])[,paralcv])
            dlDet_new[(dimlcv-1)*SGGP$numpara+paralcv] =-colSums(dSigma_now/Sigma_t)
          }
        }
        if(is.matrix(ys)){
          Sigma_t = t(matrix( rep( Sigma_t , dim(ys)[2] ) , ncol = ncol(t(Sigma_t)) , byrow = TRUE ))
        }
        
        tempvec1= (ys-yhats)/Sigma_t
        
        for (dimlcv in 1:SGGP$d) {
          for(paralcv in 1:SGGP$numpara){
            dCpn = as.matrix(dCs[,((paralcv-1)*dim(Cs)[2]+1):(paralcv*dim(Cs)[2]),dimlcv])
            if(is.matrix(dsigma2_hat)){
              if(dim(dsigma2_hat)[1]>1.5){
                dsigma2_hat_part2[(dimlcv-1)*SGGP$numpara+paralcv,] = -2*colSums((tempvec1)*(dCpn%*%pw))
              }else{
                dsigma2_hat_part2[,(dimlcv-1)*SGGP$numpara+paralcv] = -2*colSums((tempvec1)*(dCpn%*%pw))
              }
            }else{
              dsigma2_hat_part2[(dimlcv-1)*SGGP$numpara+paralcv] = -2*colSums((tempvec1)*(dCpn%*%pw))
            }
            
            
            dSigma_now = as.matrix((dSigma_to[[dimlcv]])[,paralcv])
            if(is.matrix(ys)){
              dSigma_now  = t(matrix( rep( dSigma_now  , dim(ys)[2] ) , ncol = ncol(t(dSigma_now)) , byrow = TRUE ))
            }
            
            tempvec2= dSigma_now*tempvec1
            if(is.matrix(dsigma2_hat)){
              if(dim(dsigma2_hat)[1]>1.5){
                dsigma2_hat_part1[(dimlcv-1)*SGGP$numpara+paralcv,] = -colSums(tempvec1*tempvec2)
              }else{
                dsigma2_hat_part1[,(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(tempvec1*tempvec2)
              }
            }else{
              dsigma2_hat_part1[(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(tempvec1*tempvec2)
            }
            
            
          }
        }
        if(is.vector(y)){
          temp4 = as.vector(t(Cs)%*%tempvec1)
        }else{
          temp4 = t(Cs)%*%tempvec1
        }
        ytemp = y
        dsigma2_hat_part3 =  -2*(SGGP_internal_calcusedforsupp(SGGP,ytemp,temp4, theta)$dsigma2)
        
        sigma2_hat_new = colSums((ys-yhats)*tempvec1)/dim(Xs)[1]
        lDet_new = -sum(log(Sigma_t))
        
        dsigma2_hat_new = (dsigma2_hat_part1+dsigma2_hat_part2+dsigma2_hat_part3)/dim(Xs)[1]
      }
      
      
      if(HandlingSuppData =="Correct"){
        sigma2_hat = sigma2_hat_old*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_new*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
        dsigma2_hat = dsigma2_hat_old*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+dsigma2_hat_new*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
        dlDet = dlDet_old+dlDet_new
        lDet = lDet_old+lDet_new
        
        if(!is.matrix(y)){
          neglogpost = 1/2*((dim(SGGP$design)[1]+dim(Xs)[1])*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
          gneglogpost = 0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*(length(y)+length(ys))*dsigma2_hat / sigma2_hat[1]
        }else{
          neglogpost = 1/2*((dim(SGGP$design)[1]+dim(Xs)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
          gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
          for(i in 1:dim(y)[2]){
            gneglogpost = gneglogpost + (dim(SGGP$design)[1]+dim(Xs)[1])*dsigma2_hat[,i] / sigma2_hat[i]
          }
          gneglogpost =  gneglogpost/2
        }
      }else if (HandlingSuppData =="FullValidation" || HandlingSuppData =="MarginalValidation"){
        sigma2_hat = sigma2_hat_old*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_new*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
        dsigma2_hat = dsigma2_hat_old*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+dsigma2_hat_new*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
        dlDet = dlDet_new
        lDet = lDet_new
        
        if(!is.matrix(y)){
          neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)+1/2*length(ys)*(sigma2_hat_new/sigma2_hat)
          gneglogpost = 0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*(length(ys))*dsigma2_hat / sigma2_hat[1]+ 1/2*(length(ys))*dsigma2_hat_new / sigma2_hat[1]- 1/2*(length(ys))*sigma2_hat_new*dsigma2_hat/sigma2_hat[1]^2
        }else{
          neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet)+1/2*dim(ys)[1]*sum(sigma2_hat_new/sigma2_hat)
          gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
          for(i in 1:dim(y)[2]){
            gneglogpost = gneglogpost +  dim(ys)[1]*dsigma2_hat[,i] / sigma2_hat[i]+ dim(ys)[1]*dsigma2_hat_new[,i]/sigma2_hat[i]- dim(ys)[1]*sigma2_hat_new[i]*dsigma2_hat[,i]/sigma2_hat[i]^2
          }
          gneglogpost =  gneglogpost/2
        }
        
      }}
    else{
      if(!is.matrix(y)){
        neglogpost = 1/2*(length(y)*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
        gneglogpost = 0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*(length(y))*dsigma2_hat / sigma2_hat[1]
      }else{
        neglogpost = 1/2*(dim(y)[1]*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
        gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
        for(i in 1:dim(y)[2]){
          gneglogpost = gneglogpost + dim(y)[1]*dsigma2_hat[,i] / sigma2_hat[i]
        }
        gneglogpost =  gneglogpost/2
      }
      
      if(HandlingSuppData == "Mixture" && !is.null(ys)){   
        dsigma2_hat_new = dsigma2_hat_new/dim(Xs)[1]
        
        neglogpost_old = neglogpost
        gneglogpost_old = gneglogpost
        
        if(!is.matrix(y)){
          neglogpost_new = 1/2*((length(ys))*log(sigma2_hat_new[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet_new)
          gneglogpost_new = 0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet_new+ 1/2*(length(ys))*dsigma2_hat_new / sigma2_hat_new[1]
        }else{
          neglogpost_new = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat_new)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet_new)
          gneglogpost_new = 0.5*(1/(1-theta)-1/(theta+1))+dim(ys)[2]*dlDet_new
          for(i in 1:dim(y)[2]){
            gneglogpost_new = gneglogpost_new + (dim(ys)[1])*dsigma2_hat_new[,i] / sigma2_hat_new[i]
          }
          gneglogpost_new =  gneglogpost_new/2
        }
        
        gneglogpost = gneglogpost_old + gneglogpost_new
        neglogpost = neglogpost_old + neglogpost_new
      }
    }
  }
  if(return_lik){
    return(list(neglogpost=neglogpost,gneglogpost=gneglogpost))
  } else {
    return(gneglogpost)
  }
}

