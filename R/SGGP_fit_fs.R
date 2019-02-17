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
#' SG <- SGGPcreate(d=3, batchsize=20)
#' Y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG, Y)
#' SGGP_internal_neglogpost(SG$thetaMAP, SG=SG, y=SG$y)
SGGP_internal_neglogpost <- function(theta,SGGP,y,...,ys=NULL,Xs=NULL,HandlingSuppData = "Correct") {
  # Return Inf if theta is too large
  if (max(theta) >= 1 || min(theta) <= -1) {
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
    
    if(!is.null(Xs) || HandlingSuppData == "Ignore"){
      Cs = matrix(1,dim(Xs)[1],SGGP$ss)
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        V = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb,theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Cs = Cs*V[,SGGP$designindex[,dimlcv]]
      }
      # print(Cs[1:10,1:10])
      
      pw <- SGGP_internal_calcpw(SGGP, y, theta)
      yhats = Cs%*%pw
      sigma2_hat_old = sigma2_hat
      
      if(HandlingSuppData =="Correct" || HandlingSuppData =="ValidationFull"){
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
        Sigma_t = (1-10^(-12))*Sigma_t+diag(dim(Sigma_t)[1])*10^(-12)
        
        Sigma_chol = chol(Sigma_t)
        
        Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys-yhats,transpose = TRUE))
        sigma2_hat_new = colSums((ys-yhats)*Sti_resid)/dim(Xs)[1]
        lDet_new = 2*sum(log(diag(Sigma_chol)))
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
        Sigma_t = (1-10^(-12))*Sigma_t+10^(-12)
        lDet_new = -sum(log(Sigma_t))
        if(is.matrix(ys)){
        Sigma_t = t(matrix( rep( Sigma_t , dim(ys)[2] ) , ncol = ncol(t(Sigma_t)) , byrow = TRUE ))
        }
        sigma2_hat_new = colSums((ys-yhats)^2/Sigma_t)/dim(Xs)[1]
      }
      
      
      
      if(HandlingSuppData =="Correct"){
        sigma2_hat = sigma2_hat_old*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_new*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
        
        lDet = lDet+lDet_new
        if(!is.matrix(y)){
          neglogpost = 1/2*((length(y)+length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
        }else{
          neglogpost = 1/2*((dim(y)[1]+dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
        }
        
      }else if (HandlingSuppData =="ValidationFull"){
        sigma2_hat = sigma2_hat_new
        lDet = lDet_new
        if(!is.matrix(y)){
          neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
        }else{
          neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
        }
      }else if (HandlingSuppData =="ValidationMarginal"){
        sigma2_hat = sigma2_hat_new
        lDet = lDet_new
        if(!is.matrix(y)){
          neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
        }else{
          neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
        }
      }
    }else{
      if(!is.matrix(y)){
        neglogpost = 1/2*(length(y)*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
      }else{
        neglogpost = 1/2*(dim(y)[1]*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
      }
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
#' SG <- SGGPcreate(d=3, batchsize=20)
#' Y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG, Y)
#' SGGP_internal_gneglogpost(SG$thetaMAP, SG=SG, y=SG$y)
SGGP_internal_gneglogpost <- function(theta, SGGP, y,..., return_lik=FALSE,ys=NULL,Xs=NULL,HandlingSuppData = "Correct") {
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
  
  if(!is.null(Xs) || HandlingSuppData == "Ignore"){
    Cs = matrix(1,dim(Xs)[1],SGGP$ss)
    dCs = array(1,dim=c(dim(Xs)[1],SGGP$numpara*SGGP$ss,dim(Xs)[2]))
    dsigma2_hat_old = dsigma2_hat
    sigma2_hat_old = sigma2_hat
    dlDet_old = dlDet
    
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      Vall = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb,theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
      Cs = Cs*Vall[,SGGP$designindex[,dimlcv]]
    }
    
    pw <- SGGP_internal_calcpw(SGGP, y, theta)
    yhats = Cs%*%pw
    
    Cmat1 = t(matrix( rep(Cs , SGGP$numpara ) , ncol = ncol(t(Cs)) , byrow = TRUE ))
    Cmat2 = Cmat1
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      Vall = SGGP$CorrMat(Xs[,dimlcv], SGGP$xb,theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta=TRUE)
      A = round(as.vector(outer(SGGP$designindex[,dimlcv],length(SGGP$xb)*(0:(SGGP$numpara-1)),FUN=function(x,y) x+y)))
      G = Vall$C[,SGGP$designindex[,dimlcv]]
      Cmat2 =  t(matrix( rep( G , SGGP$numpara ) , ncol = ncol(t(G)) , byrow = TRUE ))
      # print(dim(Cmat1))
      # print(dim(Cmat2))
      # print(dim(Vall$dC[,A]))
      dCs[,,dimlcv] = (Cmat1/Cmat2)*(Vall$dC[,A])
    }
    #print(Cs[1:10,1:10])
    
    if(HandlingSuppData =="Correct" || HandlingSuppData =="ValidationFull"){
      
      Sigma_t = matrix(1,dim(Xs)[1],dim(Xs)[1]) 
      dSigma_to = list(Sigma_t,SGGP$d) 
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        Vall = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara])
        Sigma_t =Sigma_t*Vall
      }
      
      Cmat1 = t(matrix( rep( Sigma_t , SGGP$numpara ) , ncol = ncol(t(Sigma_t)) , byrow = TRUE ))
      dSigma_to = list(Sigma_t,SGGP$d) 
      for (dimlcv in 1:SGGP$d) { # Loop over dimensions
        Vall = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta=TRUE)
        G = Vall$C
        Cmat2 =  t(matrix( rep( G , SGGP$numpara ) , ncol = ncol(t(G)) , byrow = TRUE ))
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
        
        ME_so1 = t(matrix( rep( (ME_s) , SGGP$numpara ) , ncol = ncol(t(ME_s)) , byrow = TRUE ))
        ME_so2 = ME_so1
        for (dimlcv in 1:SGGP$d) {
          levelnow = SGGP$uo[blocklcv,dimlcv]
          ME_n = MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
          ME_so2 =  t(matrix( rep( (ME_n) , SGGP$numpara ) , ncol = ncol(t(ME_n)) , byrow = TRUE ))
          dME_n = dMSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
          dSigma_to[[dimlcv]] = dSigma_to[[dimlcv]]-SGGP$w[blocklcv]*(ME_so1*dME_n/ME_so2)
        }
      }
      Sigma_t = (1-10^(-12))*Sigma_t+diag(dim(Sigma_t)[1])*10^(-12)
      for (dimlcv in 1:SGGP$d) {
        dSigma_to[[dimlcv]] = (1-10^(-12))*dSigma_to[[dimlcv]]
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
      
      dsigma2_hat_part3 =  -2*(SGGP_internal_calcusedforsupp(SGGP,t(Cs)%*%tempvec1,y, theta)$dsigma2)
      
      sigma2_hat_new = colSums((ys-yhats)*tempvec1)/dim(Xs)[1]
      
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
      
      Sigma_t = (1-10^(-12))*Sigma_t+10^(-12)
      for (dimlcv in 1:SGGP$d) {
        dSigma_to[[dimlcv]] = (1-10^(-12))*dSigma_to[[dimlcv]]
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
      dsigma2_hat_part3 =  -2*(SGGP_internal_calcusedforsupp(SGGP,t(Cs)%*%tempvec1,y, theta)$dsigma2)
      
      sigma2_hat_new = colSums((ys-yhats)*tempvec1)/dim(Xs)[1]
      
      dsigma2_hat_new = (dsigma2_hat_part1+dsigma2_hat_part2+dsigma2_hat_part3)/dim(Xs)[1]
    }
    
    
    if(HandlingSuppData =="Correct"){
      sigma2_hat = sigma2_hat_old*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_new*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
      dsigma2_hat = dsigma2_hat_old*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+dsigma2_hat_new*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
      dlDet = dlDet_old+dlDet_new
      
      if(!is.matrix(y)){
        neglogpost = 1/2*((length(y)+length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
        gneglogpost = 0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*(length(y)+length(ys))*dsigma2_hat / sigma2_hat[1]
      }else{
        neglogpost = 1/2*((dim(y)[1]+dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
        gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
        for(i in 1:dim(y)[2]){
          gneglogpost = gneglogpost + (dim(y)[1]+dim(ys)[1])*dsigma2_hat[,i] / sigma2_hat[i]
        }
        gneglogpost =  gneglogpost/2
      }
    }else if (HandlingSuppData =="ValidationFull"){
      sigma2_hat = sigma2_hat_new
      dsigma2_hat = dsigma2_hat_new
      dlDet = dlDet_new
      
      if(!is.matrix(y)){
        neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
      }else{
        neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet)
        gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(ys)[2]*dlDet
        for(i in 1:dim(y)[2]){
          gneglogpost = gneglogpost + (dim(ys)[1])*dsigma2_hat[,i] / sigma2_hat[i]
        }
        gneglogpost =  gneglogpost/2
      }
    }else if (HandlingSuppData =="ValidationMarginal"){
      sigma2_hat = sigma2_hat_new
      dsigma2_hat = dsigma2_hat_new
      dlDet = dlDet_new
      
      if(!is.matrix(y)){
        neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
        gneglogpost =0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*length(ys)*dsigma2_hat / sigma2_hat[1]
      }else{
        neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet)
        gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(ys)[2]*dlDet
        for(i in 1:dim(y)[2]){
          gneglogpost = gneglogpost + (dim(ys)[1])*dsigma2_hat[,i] / sigma2_hat[i]
        }
        gneglogpost = gneglogpost/2
      }
      
    }}else{  if(!is.matrix(y)){
      neglogpost = 1/2*(length(y)*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
      gneglogpost = 0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*(length(y)+length(ys))*dsigma2_hat / sigma2_hat[1]
    }else{
      neglogpost = 1/2*(dim(y)[1]*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
      gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
      for(i in 1:dim(y)[2]){
        gneglogpost = gneglogpost + dim(y)[1]*dsigma2_hat[,i] / sigma2_hat[i]
      }
      gneglogpost =  gneglogpost/2
    }
    }
  
  if(return_lik){
    return(list(neglogpost=neglogpost,gneglogpost=gneglogpost))
  } else {
    return(gneglogpost)
  }
}

#' Update SGGP model given data
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
#' @param use_PCA Should PCA be used if output is multivariate?
#' @param separateoutputparameterdimensions If multiple output dimensions,
#' should separate parameters be fit to each dimension?
#' @param use_progress_bar If using MCMC sampling, should a progress bar be
#' displayed?
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
SGGPfit <- function(SGGP, Y, Xs=NULL,Ys=NULL,
                    theta0 = SGGP$thetaMAP, #rep(0,SGGP$numpara*SGGP$d),
                    laplaceapprox = TRUE,
                    lower=rep(-1,SGGP$numpara*SGGP$d),upper=rep(1,SGGP$numpara*SGGP$d),
                    use_PCA=SGGP$use_PCA,
                    separateoutputparameterdimensions=is.matrix(SGGP$thetaMAP),
                    use_progress_bar=TRUE,
                    Ynew) {
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
  
  if (is.matrix(Y)) {
    SGGP$use_PCA <- use_PCA
    rm(use_PCA) # Just to make sure it uses fixed version
    if (is.null(SGGP$use_PCA)) {
      SGGP$use_PCA <- TRUE
    }
  }
  
  #first do the pre-processing
  #for cleanness: Y is always the user input, y is after transformation
  SGGP$Y = Y
  if(is.null(Xs)){ # No supplemental data
    SGGP$supplemented = FALSE
    
    if(!is.matrix(Y)){
      SGGP$mu = mean(Y)
      y = Y-SGGP$mu
    }else{ # Y is matrix
      SGGP$mu = colMeans(Y)
      if (SGGP$use_PCA) { # Use PCA
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
      } else { # No PCA
        y <- sweep(Y, 2, SGGP$mu)
        # Need to set SGGP$M somewhere so that it doesn't use transformation
      }
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
    } else{
      if (SGGP$use_PCA) {
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
      } else { # no PCA
        SGGP$mu = colMeans(Ys) # Could use Y, or colMeans(rbind(Y, Ys)), or make sure Ys is big enough for this
        y <- sweep(Y, 2, SGGP$mu)
        ys <- sweep(Ys, 2, SGGP$mu)
      }
    }
    SGGP$y = y
    SGGP$ys = ys
  }
  
  # nnn is numberofoutputparameterdimensions
  nnn <- if (separateoutputparameterdimensions && is.matrix(y)) {
    ncol(y)
  } else {
    1
  }
  
  # Fix theta0
  if (nnn > 1) {
    if (is.vector(theta0)) {
      theta0 <- matrix(theta0, nrow=length(theta0), ncol=nnn, byrow=F)
    }
  }
  
  
  if (is.matrix(Y) && !SGGP$use_PCA) {SGGP$M <- diag(ncol(y))} # Use identity transformation instead of PCA
  for (opdlcv in 1:nnn) { # output parameter dimension
    
    y.thisloop <- if (nnn==1) {y} else {y[,opdlcv]} # All of y or single column
    if (!is.null(Ys)) {ys.thisloop <- if (nnn==1) {ys} else {ys[,opdlcv]}} # All of y or single column
    theta0.thisloop <- if (nnn==1) {theta0} else {theta0[,opdlcv]}
    
    opt.out = nlminb(
      theta0.thisloop,
      objective = SGGP_internal_neglogpost,
      gradient = SGGP_internal_gneglogpost,
      lower = lower, 
      upper = upper,
      y = y.thisloop,
      SGGP = SGGP,
      control = list(rel.tol = 1e-8,iter.max = 500)
    )
    
    
    # opt.out = nlminb(
    #   opt.out$par,
    #   objective = SGGP_internal_neglogpost,
    #   gradient = SGGP_internal_gneglogpost,
    #   lower = lower, 
    #   upper = upper,
    #   y = y.thisloop,
    #   ys = ys.thisloop,
    #   Xs = Xs,
    #   SGGP = SGGP,
    #   control = list(rel.tol = 1e-8,iter.max = 500)
    # )
    
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
    thetaMAP <- opt.out$par
    sigma2MAP <- SGGP_internal_calcsigma2anddsigma2(SGGP=SGGP, y=y.thisloop, theta=thetaMAP, return_lS=FALSE)$sigma2
    # If one value, it gives it as matrix. Convert it to scalar
    if (length(sigma2MAP) == 1) {sigma2MAP <- sigma2MAP[1,1]}
    pw <- SGGP_internal_calcpw(SGGP=SGGP, y.thisloop, theta=thetaMAP)
    totnumpara = length(thetaMAP)
    
    # H is the Hessian at thetaMAP with reverse transformation
    H = matrix(0,nrow=totnumpara,ncol=totnumpara)
    # Transform so instead of -1 to 1 it is -Inf to Inf. Mostly in -5 to 5 though.
    PSTn=  log((1+thetaMAP)/(1-thetaMAP))
    # Reverse transformation
    thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
    # Grad of reverse transformation function
    grad0 = SGGP_internal_gneglogpost(thetav,SGGP,y.thisloop)*(2*(exp(PSTn))/(exp(PSTn)+1)^2)
    for(c in 1:totnumpara){
      rsad = rep(0,totnumpara)
      rsad[c] =10^(-12)
      PSTn=  log((1+thetaMAP)/(1-thetaMAP)) + rsad
      thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
      H[c,] = (SGGP_internal_gneglogpost(thetav,SGGP,y.thisloop)*(2*(exp(PSTn))/(exp(PSTn)+1)^2)-grad0 )*10^(4)
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
      U <- function(re){
        PSTn = log((1+thetaMAP)/(1-thetaMAP))+cHa%*%as.vector(re)
        thetav = (exp(PSTn)-1)/(exp(PSTn)+1)
        return(SGGP_internal_neglogpost(thetav,SGGP,y.thisloop))
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
    
    if(SGGP$supplemented){
      thetaMAP = thetaMAP*0.2
      for (k in 1:(SGGP$d * SGGP$numpara)){
        thetaADJ = 0*thetaMAP
        thetaADJ[k] = 10^(-3)

        print(k)
        print(10^3*(SGGP_internal_neglogpost(thetaMAP+thetaADJ,SGGP,y.thisloop)-SGGP_internal_neglogpost(thetaMAP,SGGP,y.thisloop)))
        print(SGGP_internal_gneglogpost(thetaMAP,SGGP,y.thisloop)[k])
        print(10^3*(SGGP_internal_neglogpost(thetaMAP+thetaADJ,SGGP,y.thisloop,Xs=Xs,ys=ys.thisloop,HandlingSuppData="Correct")
                    -SGGP_internal_neglogpost(thetaMAP,SGGP,y.thisloop,Xs=Xs,ys=ys.thisloop,HandlingSuppData="Correct")))
        print(SGGP_internal_gneglogpost(thetaMAP,SGGP,y.thisloop,Xs=Xs,ys=ys.thisloop,HandlingSuppData="Correct")[k])
      }
      
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
          MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]] =(-SGGP_internal_postvarmatcalc(SGGP$Xs[,dimlcv],SGGP$Xs[,dimlcv],
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
    if (nnn==1) { # Only 1 output parameter dim, so just set them
      SGGP$thetaMAP <- thetaMAP
      SGGP$sigma2MAP <- sigma2MAP
      SGGP$pw <- pw
      SGGP$thetaPostSamples <- thetaPostSamples
      if (SGGP$supplemented) {
        SGGP$pw_uppadj <- pw_uppadj
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
        SGGP$pw <- matrix(NaN, length(pw), nnn) 
        # thetaPostSamples is matrix, so this is 3dim array below
        SGGP$thetaPostSamples <- array(data = NaN, dim=c(dim(thetaPostSamples), nnn))
      }
      SGGP$thetaMAP[,opdlcv] <- thetaMAP
      SGGP$sigma2MAP[opdlcv] <- sigma2MAP
      SGGP$pw[,opdlcv] <- pw
      SGGP$thetaPostSamples[,,opdlcv] <- thetaPostSamples
      if (SGGP$supplemented) {
        if (opdlcv==1) { # First time initialize all
          
          SGGP$pw_uppadj <- matrix(NaN, nrow(pw_uppadj), nnn)
          SGGP$supppw <- matrix(NaN, nrow(supppw), nnn)
          SGGP$Sti = array(NaN, dim=c(dim(Sti), nnn)) # Sti is matrix, so this is 3 dim array
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
      # browser()
      # print(dim(dS))
      dSigma_mat = matrix(0,dim(Sigma_mat)[1],dim(Sigma_mat)[2]*length(theta))
      for(k in 1:length(theta)){
        CoinvC1oE = as.matrix(dS[,(n*(k-1)+1):(k*n)])%*%CoinvC1o
        dCoinvC1o = -backsolve(cholS,backsolve(cholS,CoinvC1oE, transpose = TRUE))
        dCoinvC1o = dCoinvC1o + backsolve(cholS,backsolve(cholS,t(dC1o[,(n*(k-1)+1):(k*n)]), transpose = TRUE))
        dSigma_mat[,((k-1)*dim(Sigma_mat)[2]+1):(k*dim(Sigma_mat)[2])] =-t(dCoinvC1o)%*%t(C2o)-t(CoinvC1o)%*%t(dC2o[,(n*(k-1)+1):(k*n)])
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

