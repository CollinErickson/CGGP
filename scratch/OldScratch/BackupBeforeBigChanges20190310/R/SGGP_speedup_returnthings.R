SGGP_internal_faststuff1 <- function(SGGP,y,theta) {
  #We need to return pw, sigma2 and cholS and lS
  
  Q  = max(SGGP$uo[1:SGGP$uoCOUNT,]) # Max value of all blocks
  
  cholS = list(matrix(1,1,1),Q*SGGP$d) # To store choleskys
  lS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$d) # Save log determinant of matrices
  
  # Loop over each dimension
  for (dimlcv in 1:SGGP$d) {
    # Loop over each possible needed correlation matrix
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      Xbrn = SGGP$xb[1:SGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = FALSE)
      S = Sstuff
      # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
      solvetry <- try({
        cS = chol(S)
        cholS[[(dimlcv-1)*Q+levellcv]]= as.matrix(cS+t(cS)-diag(diag(cS))) #store the symmetric version for C code
      })
      if (inherits(solvetry, "try-error")) {return(Inf)}
      lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
    }
  }
  
  if(!is.matrix(y)){
    sigma2 = 0 # Predictive weight for each measured point
    pw = rep(0, length(y)) # Predictive weight for each measured point
    # Loop over blocks selected
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv])>0.5){
        IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
        B0 = y[IS]
        B = SGGP$w[blocklcv]*B0
        VVV1 = unlist(cholS[gg+SGGP$uo[blocklcv,]])
        VVV2 = SGGP$gridsizest[blocklcv,]
        rcpp_kronDBS(VVV1, B, VVV2)
        pw[IS] = pw[IS]+B
        sigma2 = sigma2 + sum(B0*B)
      }
    }
    sigma2=sigma2/length(y)
    
    return(list(sigma2=sigma2,pw=pw,cholS=cholS, lS=lS))
  }else{
    numout = dim(y)[2]
    sigma2 = rep(0,numout) # Predictive weight for each measured point
    numout = dim(y)[2]
    pw = matrix(0,nrow=dim(y)[1],ncol=numout) # Predictive weight for each measured point
    # Loop over blocks selected
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv])>0.5){
        IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
        VVV1 = unlist(cholS[gg+SGGP$uo[blocklcv,]]);
        VVV2 = SGGP$gridsizest[blocklcv,];
        for(outdimlcv in 1:numout){
          B0 = y[IS,outdimlcv]
          B = SGGP$w[blocklcv]*B0
          rcpp_kronDBS(VVV1, B, VVV2)
          pw[IS,outdimlcv] = pw[IS,outdimlcv]+B
          sigma2[outdimlcv] = sigma2[outdimlcv]  + (t(B0)%*%B)
        }
      }
    }
    sigma2=sigma2/dim(y)[1]
    
    return(list(sigma2=sigma2,pw=pw,cholS=cholS,lS=lS))
  }
}


SGGP_internal_faststuff2 <- function(SGGP,y,theta) {
  #We need to return pw, sigma2, dsigma2, cholS, dMatdtheta and lS
  
  Q  = max(SGGP$uo[1:SGGP$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*SGGP$d) # To store choleskys
  dMatdtheta = list(matrix(1,1,1),Q*SGGP$d)
  lS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$d) # Save log determinant of matrices
  dlS = matrix(0, nrow = max(SGGP$uo[1:SGGP$uoCOUNT,]), ncol = SGGP$numpara*SGGP$d) 
  
  for (dimlcv in 1:SGGP$d) {
    for (levellcv in 1:max(SGGP$uo[1:SGGP$uoCOUNT,dimlcv])) {
      Xbrn = SGGP$xb[1:SGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = SGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      cS = chol(S)
      
      cholS[[(dimlcv-1)*Q+levellcv]] = as.matrix(cS+t(cS)-diag(diag(cS)))#store the symmetric version for C code
      dMatdtheta[[(dimlcv-1)*Q+levellcv]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
      for(paralcv in 1:SGGP$numpara){
        dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv] = t(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv])
      }
      lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
      for(paralcv in 1:SGGP$numpara){
        dlS[levellcv, SGGP$numpara*(dimlcv-1)+paralcv] = -sum(diag(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv]))
      }
    }
  }
  
  if(is.matrix(y)){
    numout = dim(y)[2]
    sigma2 = rep(0,numout) # Predictive weight for each measured point
    dsigma2 = matrix(0,nrow=SGGP$numpara*SGGP$d,ncol=numout) # Predictive weight for each measured point
    pw = matrix(0,nrow=dim(y)[1],ncol=numout) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv])>0.5){
        IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
        VVV1=unlist(cholS[gg+SGGP$uo[blocklcv,]])
        VVV2=unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]])
        VVV3=SGGP$gridsizest[blocklcv,]
        for(outdimlcv in 1:numout){
          B0 = y[IS,outdimlcv]
          B = SGGP$w[blocklcv]*B0
          dB = rcpp_gkronDBS(VVV1,VVV2,B,VVV3)
          pw[IS,outdimlcv] = pw[IS,outdimlcv]+B
          dsigma2[,outdimlcv] = dsigma2[,outdimlcv] + as.vector(dB%*%B0)
          sigma2[outdimlcv] = sigma2[outdimlcv]  + sum(B0*B)
        }
      }
    }
    out <- list(sigma2=sigma2/dim(y)[1],dsigma2=dsigma2/dim(y)[1],lS=lS,dlS=dlS,pw=pw,cholS=cholS,dMatdtheta=dMatdtheta)
  }else{
    sigma2 = 0 # Predictive weight for each measured point
    dsigma2 = rep(0,nrow=SGGP$d) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    pw = rep(0, length(y)) # Predictive weight for each measured point
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv])>0.5){
        IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
        B0 = y[IS]
        B = SGGP$w[blocklcv]*B0
        dB = rcpp_gkronDBS(unlist(cholS[gg+SGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]]), B, SGGP$gridsizest[blocklcv,])
        
        pw[IS] = pw[IS]+B
        dsigma2 = dsigma2 +t(B0)%*%t(dB)
        sigma2 = sigma2 + t(B0)%*%B
      }
    }
    out <- list(sigma2=sigma2/length(y),dsigma2=dsigma2/length(y),lS=lS,dlS=dlS,pw=pw,cholS=cholS,dMatdtheta=dMatdtheta)
  }
  out
}



SGGP_internal_faststuff3 <- function(SGGP,revc,y,cholS,dMatdtheta) {
  Q  = max(SGGP$uo[1:SGGP$uoCOUNT,]) # Max level of all blocks
  
  if(is.matrix(y)){
    numout = dim(y)[2]
    valo = rep(0,numout) # Predictive weight for each measured point
    dvalo = matrix(0,nrow=SGGP$numpara*SGGP$d,ncol=numout) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv])>0.5){
        IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
        VVV1=unlist(cholS[gg+SGGP$uo[blocklcv,]])
        VVV2=unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]])
        VVV3=SGGP$gridsizest[blocklcv,]
        for(outdimlcv in 1:numout){
          B2 = y[IS,outdimlcv]
          B0 = revc[IS,outdimlcv]
          B = (SGGP$w[blocklcv])*B0#/dim(y)[1]
          dB = rcpp_gkronDBS(VVV1,VVV2,B,VVV3)
          dvalo[,outdimlcv] = dvalo[,outdimlcv] +t(B2)%*%t(dB)
          valo[outdimlcv] = valo[outdimlcv]  + sum(B2*B)+ t(B2)%*%B
        }
      }
    }
    out <- list(valo=valo,dvalo=dvalo)
  }else{
    valo= 0 # Predictive weight for each measured point
    dvalo = rep(0,nrow=SGGP$d) # Predictive weight for each measured point
    gg = (1:SGGP$d-1)*Q
    for (blocklcv in 1:SGGP$uoCOUNT) {
      if(abs(SGGP$w[blocklcv])>0.5){
        IS = SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]];
        B0 = revc[IS]
        B2 = y[IS]
        B = (SGGP$w[blocklcv])*B0#/length(y)
        dB = rcpp_gkronDBS(unlist(cholS[gg+SGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+SGGP$uo[blocklcv,]]), B, SGGP$gridsizest[blocklcv,])
        
        dvalo = dvalo + as.vector(dB%*%B2)
        valo = valo + sum(B2*B)
      }
    }
    out <- list(valo=valo,dvalo=dvalo)
  }
  out
}