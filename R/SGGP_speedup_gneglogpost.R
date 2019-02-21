

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
SGGP_internal_gneglogpost <- function(theta, SGGP, y,..., return_lik=FALSE,ys=NULL,Xs=NULL,HandlingSuppData = "Correct") {
  
  if (!(HandlingSuppData %in% c("Correct", "Only", "Ignore", "Mixture", "MarginalValidation", "FullValidation"))) {
    stop(paste("HandlingSuppData in SGGP_internal_neglogpost must be one of",
               "Correct, Only, Ignore, Mixture, MarginalValidation, FullValidation"))
  }
  
  epsssV = 0 # 10^(-14)
  
  
  if(!is.null(Xs) && is.null(y)){
    HandlingSuppData = "Only"
  }else if(is.null(Xs) && !is.null(y)){
    HandlingSuppData = "Ignore"
  }else if(is.null(Xs) && is.null(y)){
    stop(paste("You have given no y or ys to SGGP_internal_neglogpost"))
  }
  
  
  if(HandlingSuppData == "Only" || 
     HandlingSuppData == "Mixture" || 
     HandlingSuppData == "FullValidation" || 
     HandlingSuppData == "Correct"){
    
    Sigma_t = matrix(0,dim(Xs)[1],dim(Xs)[1])
    dSigma_to = list(Sigma_t,SGGP$d) 
    RTR = list(Sigma_t,SGGP$d) 
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      V = SGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv], theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta=TRUE,returnlogs=TRUE)
      dSigma_to[[dimlcv]] = V$dCdtheta
      Sigma_t = Sigma_t+V$C
    }
    Sigma_t = exp(Sigma_t)
    
    Cmat1 = matrix( rep(Sigma_t , SGGP$numpara ) , nrow = nrow(Sigma_t) , byrow = FALSE )
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      dSigma_to[[dimlcv]] =Cmat1*dSigma_to[[dimlcv]]
    }
    
    Sigma_t = (1-epsssV)*Sigma_t+diag(dim(Sigma_t)[1])*epsssV
    for (dimlcv in 1:SGGP$d) {
      dSigma_to[[dimlcv]] = (1-epsssV)*dSigma_to[[dimlcv]]
    }
  }
  
  
  if(HandlingSuppData == "Only" || 
     HandlingSuppData == "Mixture"){
    Sigma_chol = chol(Sigma_t)
    tempvec1 = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose = TRUE))
    sigma2_hat_supp = colSums(as.matrix(ys*tempvec1))/dim(Xs)[1]
    if(is.matrix(ys)){
      dsigma2_hat_supp = matrix(0,SGGP$d*SGGP$numpara,ncol=dim(ys)[2])
    }else{
      dsigma2_hat_supp = rep(0,SGGP$d*SGGP$numpara)
    }
    
    dlDet_supp=rep(0,SGGP$d*SGGP$numpara)
    for (dimlcv in 1:SGGP$d) {
      for(paralcv in 1:SGGP$numpara){
        dSigma_supp = as.matrix((dSigma_to[[dimlcv]])[,((paralcv-1)*dim(Sigma_chol)[2]+1):(paralcv*dim(Sigma_chol)[2])])
        tempvec2= dSigma_supp%*%tempvec1
        if(is.matrix(dsigma2_hat_supp )){
          if(dim(dsigma2_hat_supp)[1]>1.5){
            dsigma2_hat_supp[(dimlcv-1)*SGGP$numpara+paralcv,] = -colSums(as.matrix(tempvec1*tempvec2))/dim(Xs)[1]
          }else{
            dsigma2_hat_supp[,(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(as.matrix(tempvec1*tempvec2))/dim(Xs)[1]
          }
        }else{
          dsigma2_hat_supp[(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(as.matrix(tempvec1*tempvec2))/dim(Xs)[1]
        }
        dlDet_supp[(dimlcv-1)*SGGP$numpara+paralcv] = sum(diag(backsolve(Sigma_chol,backsolve(Sigma_chol,dSigma_supp,transpose = TRUE))))
      }
    }
    lDet_supp = 2*sum(log(diag(Sigma_chol)))
    
  }
  
  if(HandlingSuppData == "Ignore" ||
     HandlingSuppData == "Mixture" ){
    
    sigma2anddsigma2 <- SGGP_internal_calcsigma2anddsigma2(SGGP=SGGP, y=y, theta=theta, return_lS=TRUE)
    lS <- sigma2anddsigma2$lS
    dlS <-sigma2anddsigma2$dlS
    
    sigma2_hat_grid = sigma2anddsigma2$sigma2
    dsigma2_hat_grid = sigma2anddsigma2$dsigma2
    
    lDet_grid = 0
    dlDet_grid = rep(0, SGGP$numpara*SGGP$d) 
    
    for (blocklcv in 1:SGGP$uoCOUNT) {
      nv = SGGP$gridsize[blocklcv]/SGGP$gridsizes[blocklcv,]
      uonow = SGGP$uo[blocklcv,]
      for (dimlcv in which(uonow>1.5)) {
        if (return_lik) {
          lDet_grid = lDet_grid + (lS[uonow[dimlcv], dimlcv] - lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
        }
        IS = (dimlcv-1)*SGGP$numpara+1:SGGP$numpara
        dlDet_grid[IS] = dlDet_grid[IS] + (dlS[uonow[dimlcv], IS] - dlS[uonow[dimlcv]-1, IS])*nv[dimlcv]
      }
    }
  }
  if(HandlingSuppData == "FullValidation" || 
     HandlingSuppData == "Correct" || 
     HandlingSuppData == "MarginalValidation" ){
    Cs = matrix(0,dim(Xs)[1],SGGP$ss)
    dCs = array(0,dim=c(dim(Xs)[1],SGGP$numpara*SGGP$ss,dim(Xs)[2]))
    
    GGGG = list(matrix(1,dim(Xs)[1],length(SGGP$xb)),SGGP$d)
    dGGGG1 = list(matrix(1,dim(Xs)[1],length(SGGP$xb)),SGGP$d)
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      Xbn = SGGP$xb[1:SGGP$sizest[max(SGGP$uo[,dimlcv])]]
      V = SGGP$CorrMat(Xs[,dimlcv],  Xbn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta=TRUE,returnlogs=TRUE)
      A = round(as.vector(outer(SGGP$designindex[,dimlcv],length(Xbn)*(0:(SGGP$numpara-1)),FUN=function(x,y) x+y)))
      GGGG[[dimlcv]] = exp(V$C)
      dGGGG1[[dimlcv]] = V$dCdtheta
      dCs[,,dimlcv] = V$dCdtheta[,A]
      Cs= Cs+V$C[,SGGP$designindex[,dimlcv]]
    }
    Cs = exp(Cs)
    
    Cmat1 = matrix( rep(Cs, SGGP$numpara ) , nrow = nrow(Cs) , byrow = FALSE )
    for (dimlcv in 1:SGGP$d) { # Loop over dimensions
      dCs[,,dimlcv] =Cmat1*(dCs[,,dimlcv])
      Cmat2 = matrix(rep(GGGG[[dimlcv]], SGGP$numpara ) ,nrow = nrow(GGGG[[dimlcv]]) , byrow = FALSE )
      dGGGG1[[dimlcv]] =Cmat2*(dGGGG1[[dimlcv]])
    }
    
    lik_stuff <- SGGP_internal_faststuff2(SGGP=SGGP, y=y, theta=theta)
    cholS = lik_stuff$cholS
    dSV = lik_stuff$dMatdtheta
    lS <- lik_stuff$lS
    dlS <- lik_stuff$dlS
    sigma2_hat_grid = lik_stuff$sigma2
    dsigma2_hat_grid = lik_stuff$dsigma2
    pw = lik_stuff$pw
    
    yhats = Cs%*%pw
    
    lDet_grid = 0
    dlDet_grid = rep(0, SGGP$numpara*SGGP$d) 
    
    for (blocklcv in 1:SGGP$uoCOUNT) {
      nv = SGGP$gridsize[blocklcv]/SGGP$gridsizes[blocklcv,]
      uonow = SGGP$uo[blocklcv,]
      for (dimlcv in which(uonow>1.5)) {
        if (return_lik) {
          lDet_grid = lDet_grid + (lS[uonow[dimlcv], dimlcv] - lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
        }
        IS = (dimlcv-1)*SGGP$numpara+1:SGGP$numpara
        dlDet_grid[IS] = dlDet_grid[IS] + (dlS[uonow[dimlcv], IS] - dlS[uonow[dimlcv]-1, IS])*nv[dimlcv]
      }
    }
  }
  
  
  if(HandlingSuppData == "FullValidation" || 
     HandlingSuppData == "Correct"){
    MSE_s = list(matrix(0,dim(Xs)[1],dim(Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1))
    dMSE_s = list(matrix(0,dim(Xs)[1],dim(Xs)[1]),(SGGP$d+1)*(SGGP$maxlevel+1))
    Q  = max(SGGP$uo[1:SGGP$uoCOUNT,])
    for (dimlcv in 1:SGGP$d) {
      gg = (dimlcv-1)*Q
      for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
        INDSN = 1:SGGP$sizest[levellcv]
        INDSN = INDSN[sort(SGGP$xb[1:SGGP$sizest[levellcv]],index.return = TRUE)$ix]
        REEALL = SGGP_internal_postvarmatcalcfaster(GGGG[[dimlcv]],
                                                    dGGGG1[[dimlcv]],
                                                    as.matrix(cholS[[gg+levellcv]]),
                                                    as.matrix(dSV[[gg+levellcv]]),
                                                    INDSN,
                                                    SGGP$numpara,
                                                    returndG = TRUE,
                                                    returnderiratio =TRUE)
        MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]]  =  REEALL$Sigma_mat
        dMSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]]   = REEALL$dSigma_mat
      }
    }
    
    dsigma2_hat_part1 = 0*dsigma2_hat_grid
    dsigma2_hat_part2 = 0*dsigma2_hat_grid
    dsigma2_hat_part3 = 0*dsigma2_hat_grid
    dlDet_supp = 0*dlDet_grid
    
    for (blocklcv in 1:SGGP$uoCOUNT) {
      ME_s = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
      for (dimlcv in 1:SGGP$d) {
        levelnow = SGGP$uo[blocklcv,dimlcv]
        ME_s = ME_s*MSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
      }
      Sigma_t = Sigma_t-SGGP$w[blocklcv]*ME_s
      
      MEE_n = matrix( rep(SGGP$w[blocklcv]*ME_s, SGGP$numpara ) , nrow = nrow(ME_s) , byrow = FALSE )
      
      for (dimlcv in 1:SGGP$d) {
        levelnow = SGGP$uo[blocklcv,dimlcv]
        dME_n = MEE_n*dMSE_s[[(dimlcv)*SGGP$maxlevel+levelnow]]
        dSigma_to[[dimlcv]] = dSigma_to[[dimlcv]]-dME_n
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
        if(is.matrix(dsigma2_hat_part2)){
          if(dim(dsigma2_hat_part2)[1]>1.5){
            dsigma2_hat_part2[(dimlcv-1)*SGGP$numpara+paralcv,] = -2*colSums((tempvec1)*(dCpn%*%pw))
          }else{
            dsigma2_hat_part2[,(dimlcv-1)*SGGP$numpara+paralcv] =  -2*colSums((tempvec1)*(dCpn%*%pw))
          }
        }else{
          dsigma2_hat_part2[(dimlcv-1)*SGGP$numpara+paralcv] = -2*colSums((tempvec1)*(dCpn%*%pw))
        }
        
        
        dSigma_now = as.matrix((dSigma_to[[dimlcv]])[,((paralcv-1)*dim(Sigma_chol)[2]+1):(paralcv*dim(Sigma_chol)[2])])
        tempvec2= dSigma_now%*%tempvec1
        if(is.matrix(dsigma2_hat_part2)){
          if(dim(dsigma2_hat_part2)[1]>1.5){
            dsigma2_hat_part1[(dimlcv-1)*SGGP$numpara+paralcv,] = -colSums(tempvec1*tempvec2)
          }else{
            dsigma2_hat_part1[,(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(tempvec1*tempvec2)
          }
        }else{
          dsigma2_hat_part1[(dimlcv-1)*SGGP$numpara+paralcv] = -colSums(tempvec1*tempvec2)
        }
        
        
        dlDet_supp[(dimlcv-1)*SGGP$numpara+paralcv] = sum(diag(backsolve(Sigma_chol,backsolve(Sigma_chol,dSigma_now,transpose = TRUE))))
      }
    }
    
    if(is.vector(y)){
      temp4 = as.vector(t(Cs)%*%tempvec1)
    }else{
      temp4 = t(Cs)%*%tempvec1
    }
    dsigma2_hat_part3 =  -2*(SGGP_internal_faststuff3(SGGP,y,temp4,cholS,dSV)$dvalo)
    lDet_supp = 2*sum(log(diag(Sigma_chol)))
    sigma2_hat_supp = colSums((ys-yhats)*tempvec1)/dim(Xs)[1]
    dsigma2_hat_supp = (dsigma2_hat_part1+dsigma2_hat_part2+dsigma2_hat_part3)/dim(Xs)[1]
  }
  if(HandlingSuppData == "MarginalValidation"){
    Sigma_t = matrix(1,dim(Xs)[1],1)
    # 
    Cmat1 = 0*matrix( rep( (Sigma_t),SGGP$numpara) , nrow = nrow(Cs) , byrow = FALSE )
    dSigma_to = list(Cmat1,SGGP$d)
    for (dimlcv in 1:SGGP$d) {
      dSigma_to[[dimlcv]] = 0*Cmat1
    }
    # 
    MSE_s = list(matrix(0,dim(Xs)[1],1),(SGGP$d+1)*(SGGP$maxlevel+1))
    dMSE_s = list(matrix(0,dim(Xs)[1],SGGP$numpara),(SGGP$d+1)*(SGGP$maxlevel+1))
    QRE  = max(SGGP$uo[1:SGGP$uoCOUNT,])
    for (dimlcv in 1:SGGP$d) {
      gg = (dimlcv-1)*QRE
      for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
        INDSN = 1:SGGP$sizest[levellcv]
        INDSN = INDSN[sort(SGGP$xb[1:SGGP$sizest[levellcv]],index.return = TRUE)$ix]
        REEALL = SGGP_internal_postvarmatcalcfaster(GGGG[[dimlcv]],
                                                    dGGGG1[[dimlcv]],
                                                    as.matrix(cholS[[gg+levellcv]]),
                                                    as.matrix(dSV[[gg+levellcv]]),
                                                    INDSN,
                                                    SGGP$numpara,
                                                    returndG = TRUE,
                                                    returndiag=TRUE)
        MSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]]  = REEALL$Sigma_mat
        dMSE_s[[(dimlcv)*SGGP$maxlevel+levellcv]]   = REEALL$dSigma_mat
      }
    }
    
    
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
    
    
    dsigma2_hat_part1 = 0*dsigma2_hat_grid
    dsigma2_hat_part2 = 0*dsigma2_hat_grid
    dsigma2_hat_part3 = 0*dsigma2_hat_grid
    dlDet_supp = 0*dlDet_grid
    
    lDet_supp = -sum(log(Sigma_t))
    
    for (dimlcv in 1:SGGP$d) {
      for(paralcv in 1:SGGP$numpara){
        dSigma_now = as.matrix((dSigma_to[[dimlcv]])[,paralcv])
        dlDet_supp[(dimlcv-1)*SGGP$numpara+paralcv] =-colSums(dSigma_now/Sigma_t)
      }
    }
    if(is.matrix(ys)){
      Sigma_t = t(matrix( rep( Sigma_t , dim(ys)[2] ) , ncol = ncol(t(Sigma_t)) , byrow = TRUE ))
    }
    
    tempvec1= (ys-yhats)/Sigma_t
    for (dimlcv in 1:SGGP$d) {
      for(paralcv in 1:SGGP$numpara){
        dCpn = as.matrix(dCs[,((paralcv-1)*dim(Cs)[2]+1):(paralcv*dim(Cs)[2]),dimlcv])
        if(is.matrix(dsigma2_hat_part2)){
          if(dim(dsigma2_hat_part2)[1]>1.5){
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
        if(is.matrix(dsigma2_hat_part2)){
          if(dim(dsigma2_hat_part2)[1]>1.5){
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
    dsigma2_hat_part3 =  -2*(SGGP_internal_faststuff3(SGGP,y,temp4,cholS,dSV)$dvalo)
    
    
    sigma2_hat_supp = colSums((ys-yhats)*tempvec1)/dim(Xs)[1]
    dsigma2_hat_supp = (dsigma2_hat_part1+dsigma2_hat_part2+dsigma2_hat_part3)/dim(Xs)[1]
  }
  
  neglogpost = 0
  gneglogpost = rep(0,length(theta))
  
  if(HandlingSuppData =="Only" || HandlingSuppData == "Mixture"){
    sigma2_hat = sigma2_hat_supp
    dsigma2_hat = dsigma2_hat_supp
    dlDet = dlDet_supp
    lDet = lDet_supp
    
    if(!is.matrix(ys)){
      neglogpost = neglogpost+1/2*(dim(Xs)[1]*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
      gneglogpost = gneglogpost+0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*length(ys)*dsigma2_hat / sigma2_hat[1]
    }else{
      neglogpost = neglogpost+1/2*(dim(Xs)[1]*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet)
      gneglogpostn = 0.5*(1/(1-theta)-1/(theta+1))+dim(ys)[2]*dlDet
      for(i in 1:dim(ys)[2]){
        gneglogpostn = gneglogpostn + dim(Xs)[1]*dsigma2_hat[,i] / sigma2_hat[i]
      }
      gneglogpost = gneglogpost+  gneglogpostn/2
      
    }
  }
  if(HandlingSuppData =="Ignore" || HandlingSuppData == "Mixture"){
    sigma2_hat = sigma2_hat_grid
    dsigma2_hat = dsigma2_hat_grid
    dlDet = dlDet_grid
    lDet = lDet_grid
    
    if(!is.matrix(y)){
      neglogpost = neglogpost+1/2*(dim(SGGP$design)[1]*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)
      gneglogpost = gneglogpost+0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*length(y)*dsigma2_hat / sigma2_hat[1]
    }else{
      neglogpost = neglogpost+1/2*(dim(SGGP$design)[1]*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(y)[2]*lDet)
      gneglogpostn = 0.5*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
      for(i in 1:dim(y)[2]){
        gneglogpostn = gneglogpostn + dim(SGGP$design)[1]*dsigma2_hat[,i] / sigma2_hat[i]
      }
      gneglogpost = gneglogpost+  gneglogpostn/2
    }
  }
  
  if(HandlingSuppData =="Correct"){
    sigma2_hat = sigma2_hat_grid*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
    dsigma2_hat = dsigma2_hat_grid*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+dsigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
    dlDet = dlDet_grid+dlDet_supp
    lDet = lDet_grid+lDet_supp
    
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
  }
  
  if (HandlingSuppData =="FullValidation" || HandlingSuppData =="MarginalValidation"){
    sigma2_hat = sigma2_hat_grid*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+sigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
    dsigma2_hat = dsigma2_hat_grid*dim(SGGP$design)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])+dsigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(SGGP$design)[1])
    dlDet = dlDet_supp
    lDet = lDet_supp
    
    if(!is.matrix(y)){
      neglogpost = 1/2*((length(ys))*log(sigma2_hat[1])-0.500*sum(log(1-theta)+log(theta+1))+lDet)+1/2*length(ys)*(sigma2_hat_supp/sigma2_hat)
      gneglogpost = 0.25*(1/(1-theta)-1/(theta+1))+ 1/2*dlDet+ 1/2*(length(ys))*dsigma2_hat / sigma2_hat[1]+ 1/2*(length(ys))*dsigma2_hat_supp / sigma2_hat[1]- 1/2*(length(ys))*sigma2_hat_supp*dsigma2_hat/sigma2_hat[1]^2
    }else{
      neglogpost = 1/2*((dim(ys)[1])*sum(log(c(sigma2_hat)))-0.500*sum(log(1-theta)+log(theta+1))+dim(ys)[2]*lDet)+1/2*dim(ys)[1]*sum(sigma2_hat_supp/sigma2_hat)
      gneglogpost = 0.5*(1/(1-theta)-1/(theta+1))+dim(y)[2]*dlDet
      for(i in 1:dim(y)[2]){
        gneglogpost = gneglogpost +  dim(ys)[1]*dsigma2_hat[,i] / sigma2_hat[i]+ dim(ys)[1]*dsigma2_hat_supp[,i]/sigma2_hat[i]- dim(ys)[1]*sigma2_hat_supp[i]*dsigma2_hat[,i]/sigma2_hat[i]^2
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