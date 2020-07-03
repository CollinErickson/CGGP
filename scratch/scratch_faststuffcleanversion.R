# Trying to make single function that will only calculate needed things
# Gets super messy, going to ditch it
# It's easier to just have 3 separate functions

CGGP_internal_faststuff123 <- function(CGGP,y,theta, 
                                       revc, cholS, dMatdtheta,
                                       r_sigma2, r_dsigma2, r_pw,
                                       r_cholS, r_lS, r_dlS, r_dMatdtheta,
                                       r_valo) {
  #We need to return pw, sigma2, dsigma2, cholS, dMatdtheta and lS
  
  Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max level of all blocks
  if (r_cholS) {
    cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
  }
  if (r_lS) {
    lS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$d) # Save log determinant of matrices
  }
  if (r_dMatdtheta || r_dlS) {
    dMatdtheta = list(matrix(1,1,1),Q*CGGP$d)
    dlS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$numpara*CGGP$d) 
  }
  
  if (any(r_sigma2, r_dsigma2, r_pw, r_lS, r_dlS, r_dMatdtheta)) { # None of this for 3
    for (dimlcv in 1:CGGP$d) {
      for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
        Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
        Xbrn = Xbrn[order(Xbrn)]
        
        Sstuff = CGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = r_dMatdtheta || r_dlS)
        S <- if (r_dCdt) Sstuff else Sstuff$C
        if (r_cholS) {
          solvetry <- try({cS <- chol(S)})
          if (inherits(solvetry, "try-error")) {return(Inf)}
          cholS[[(dimlcv-1)*Q+levellcv]] = as.matrix(cS+t(cS)-diag(diag(cS)))
        }
        
        if (r_lS) {lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))}
        if (r_dMatdtheta || r_dlS) {
          dMatdtheta[[(dimlcv-1)*Q+levellcv]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
          nv = length(Xbrn)
          for(paralcv in 1:CGGP$numpara){
            dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv] = t(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv])
          }
          for(paralcv in 1:CGGP$numpara){
          if(nv > 1.5){
            dlS[levellcv, CGGP$numpara*(dimlcv-1)+paralcv] = -sum(diag(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv]))
          } else {
            dlS[levellcv, CGGP$numpara*(dimlcv-1)+paralcv] = -dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv]
          }
          }
        }
      }
    }
  }
  
  if(is.matrix(y)){
    numout = dim(y)[2]
    if (r_sigma2) sigma2 = rep(0,numout) # Predictive weight for each measured point
    if (r_dsigma2) dsigma2 = matrix(0,nrow=CGGP$numpara*CGGP$d,ncol=numout) # Predictive weight for each measured point
    if (r_pw) pw = matrix(0,nrow=dim(y)[1],ncol=numout) # Predictive weight for each measured point
    if (r_valo) {
      valo = rep(0,numout) # Predictive weight for each measured point
      dvalo = matrix(0,nrow=CGGP$numpara*CGGP$d,ncol=numout) # Predictive weight for each measured point
    }
    
    # Loop over blocks selected
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        VVV1=unlist(cholS[gg+CGGP$uo[blocklcv,]])
        VVV2=unlist(dMatdtheta[gg+CGGP$uo[blocklcv,]])
        if (num!=1) VVV3=CGGP$gridsizest[blocklcv,]
        for(outdimlcv in 1:numout){
          if (num!=3) B0_12 = y[IS,outdimlcv]
          if (num==3) B0_3 = revc[IS,outdimlcv]
          B_12 = CGGP$w[blocklcv]*B0_12
          if (num!=1) dB = rcpp_gkronDBS(VVV1,VVV2,B,VVV3)
          if (num!=3) pw[IS,outdimlcv] = pw[IS,outdimlcv]+B
          if (num==2) dsigma2[,outdimlcv] = dsigma2[,outdimlcv] + as.vector(dB%*%B0)
          if (num!=3) sigma2[outdimlcv] = sigma2[outdimlcv]  + sum(B0_12*B_12)
          if (num==3) dvalo[,outdimlcv] = dvalo[,outdimlcv] +t(B2)%*%t(dB)
          if (num==3) valo[outdimlcv] = valo[outdimlcv]  + sum(B2*B)+ t(B2)%*%B
        }
      }
    }
    if (num==1) out <- list(sigma2=sigma2/dim(y)[1],pw=pw,cholS=cholS,lS=lS)
    if (num==2) out <- list(sigma2=sigma2/dim(y)[1],dsigma2=dsigma2/dim(y)[1],lS=lS,dlS=dlS,pw=pw,cholS=cholS,dMatdtheta=dMatdtheta)
    if (num==3) out <- list(valo=valo,dvalo=dvalo)
  }else{ # !is.matrix(y)
    if (num!=3) sigma2 = 0 # Predictive weight for each measured point
    if (num==2) dsigma2 = rep(0,nrow=CGGP$d) # Predictive weight for each measured point
    # Loop over blocks selected
    if (num==3) valo= 0 # Predictive weight for each measured point
    if (num==3) dvalo = rep(0,nrow=CGGP$d) # Predictive weight for each measured point
    gg = (1:CGGP$d-1)*Q
    if (num!=3) pw = rep(0, length(y)) # Predictive weight for each measured point
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        if (num!=3) B0 = y[IS]
        if (num==3) B0 = revc[IS]
        if (num==3) B2 = y[IS]
        B = CGGP$w[blocklcv]*B0
        if (num!=1) dB = rcpp_gkronDBS(unlist(cholS[gg+CGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
        
        if (num!=3) pw[IS] = pw[IS]+B
        if (num==2) dsigma2 = dsigma2 +t(B0)%*%t(dB)
        if (num==2) sigma2 = sigma2 + t(B0)%*%B
        if (num==1) sigma2 = sigma2 + sum(B0*B)
        if (num==3) dvalo = dvalo + as.vector(dB%*%B2)
        if (num==3) valo = valo + sum(B2*B)
      }
    }
    if (num==1) out <- list(sigma2=sigma2,pw=pw,cholS=cholS, lS=lS)
    if (num==2) out <- list(sigma2=sigma2/length(y),dsigma2=dsigma2/length(y),lS=lS,dlS=dlS,pw=pw,cholS=cholS,dMatdtheta=dMatdtheta)
    if (num==3) out <- list(valo=valo,dvalo=dvalo)
  }
  out
}


