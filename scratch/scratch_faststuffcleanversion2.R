

SGGP_internal_faststuff123 <- function(SGGP,y,theta, num) {
  #We need to return pw, sigma2, dsigma2, cholS, dMatdtheta and lS
  
  Q  = max(SGGP$uo[1:SGGP$uoCOUNT,]) # Max level of all blocks
  if (num==1 || num==2) {
    cholS = list(matrix(1,1,1),Q*SGGP$d) # To store choleskys
    lS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$d) # Save log determinant of matrices
  }
  if (num==2) {
    dMatdtheta = list(matrix(1,1,1),Q*SGGP$d)
    dlS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$numpara*SGGP$d) 
  }
  
  if (num==1 || num==2) { # None of this for 3
    for (dimlcv in 1:SGGP$d) {
      for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
        Xbrn = SGGP$xb[1:SGGP$sizest[levellcv]]
        Xbrn = Xbrn[order(Xbrn)]
        if (num==1) {
          Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = FALSE)
          S = Sstuff
          # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
          solvetry <- try({
            cS = chol(S)
            cholS[[(dimlcv-1)*Q+levellcv]]= as.matrix(cS+t(cS)-diag(diag(cS))) #store the symmetric version for C code
          })
          if (inherits(solvetry, "try-error")) {return(Inf)}
        }
        if (num==2) {
          nv = length(Xbrn)
          Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = TRUE)
          S = Sstuff$C
          cS = chol(S)
          cholS[[(dimlcv-1)*Q+levellcv]] = as.matrix(cS+t(cS)-diag(diag(cS)))#store the symmetric version for C code
        }
        lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
        if (num==2) {
          dMatdtheta[[(dimlcv-1)*Q+levellcv]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
          for(paralcv in 1:SGGP$numpara){
            dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv] = t(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv])
          }
          for(paralcv in 1:SGGP$numpara){
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
    if (num==1 || num==2) sigma2 = rep(0,numout) # Predictive weight for each measured point
    if (num==2) dsigma2 = matrix(0,nrow=SGGP$numpara*SGGP$d,ncol=numout) # Predictive weight for each measured point
    if (num==1 || num==2) pw = matrix(0,nrow=dim(y)[1],ncol=numout) # Predictive weight for each measured point
    if (num==3) {
      valo = rep(0,numout) # Predictive weight for each measured point
      dvalo = matrix(0,nrow=SGGP$numpara*SGGP$d,ncol=numout) # Predictive weight for each measured point
    }
    
    # Loop over blocks selected
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv])>0.5){
        IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
        VVV1=unlist(cholS[gg+SGGP$uo[blocklcv,]])
        VVV2=unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]])
        if (num!=1) VVV3=SGGP$gridsizest[blocklcv,]
        for(outdimlcv in 1:numout){
          if (num!=3) B0 = y[IS,outdimlcv]
          if (num==3) B0 = revc[IS,outdimlcv]
          B = SGGP$w[blocklcv]*B0
          if (num!=1) dB = rcpp_gkronDBS(VVV1,VVV2,B,VVV3)
          if (num!=3) pw[IS,outdimlcv] = pw[IS,outdimlcv]+B
          if (num==2) dsigma2[,outdimlcv] = dsigma2[,outdimlcv] + as.vector(dB%*%B0)
          if (num!=3) sigma2[outdimlcv] = sigma2[outdimlcv]  + sum(B0*B)
          if (num==3) dvalo[,outdimlcv] = dvalo[,outdimlcv] +t(B2)%*%t(dB)
          if (num==3) valo[outdimlcv] = valo[outdimlcv]  + sum(B2*B)+ t(B2)%*%B
        }
      }
    }
    if (num==1) out <- list(sigma2=sigma2/dim(y)[1],pw=pw,cholS=cholS,lS=lS)
    if (num==2) out <- list(sigma2=sigma2/dim(y)[1],dsigma2=dsigma2/dim(y)[1],lS=lS,dlS=dlS,pw=pw,cholS=cholS,dMatdtheta=dMatdtheta)
    if (num==3) out <- list(valo=valo,dvalo=dvalo)
  }else{
    if (num!=3) sigma2 = 0 # Predictive weight for each measured point
    if (num==2) dsigma2 = rep(0,nrow=SGGP$d) # Predictive weight for each measured point
    # Loop over blocks selected
    if (num==3) valo= 0 # Predictive weight for each measured point
    if (num==3) dvalo = rep(0,nrow=SGGP$d) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    if (num!=3) pw = rep(0, length(y)) # Predictive weight for each measured point
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv])>0.5){
        IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
        if (num!=3) B0 = y[IS]
        if (num==3) B0 = revc[IS]
        if (num==3) B2 = y[IS]
        B = SGGP$w[blocklcv]*B0
        if (num!=1) dB = rcpp_gkronDBS(unlist(cholS[gg+SGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]]), B, SGGP$gridsizest[blocklcv,])
        
        if (num!=3) pw[IS] = pw[IS]+B
        if (num==2) dsigma2 = dsigma2 +t(B0)%*%t(dB)
        if (num==2) sigma2 = sigma2 + t(B0)%*%B
        if (num==1) sigma2 = sigma2 + sum(B0*B)
      }
    }
    if (num==1) out <- list(sigma2=sigma2/length(y),pw=pw,cholS=cholS, lS=lS)
    if (num==2) out <- list(sigma2=sigma2/length(y),dsigma2=dsigma2/length(y),lS=lS,dlS=dlS,pw=pw,cholS=cholS,dMatdtheta=dMatdtheta)
    if (num==3) out <- list(valo=valo,dvalo=dvalo)
  }
  out
}


