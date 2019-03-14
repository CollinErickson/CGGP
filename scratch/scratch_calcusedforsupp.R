# This was in fastcalcassist, but was never used.
# Moved to scratch.

CGGP_internal_calcusedforsupp <- function(CGGP, revc, y, theta, return_lS=FALSE) {
  Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
  dMatdtheta = list(matrix(1,1,1),Q*CGGP$d)
  
  if(return_lS){
    lS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$d) # Save log determinant of matrices
    dlS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$numpara*CGGP$d) 
  }
  
  # Loop over each dimension
  for (dimlcv in 1:CGGP$d) {
    # Loop over depth of each dim
    for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
      Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = CGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      
      cS = chol(S)
      
      cholS[[(dimlcv-1)*Q+levellcv]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      dMatdtheta[[(dimlcv-1)*Q+levellcv]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
      for(paralcv in 1:CGGP$numpara){
        dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv] = t(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv])
      }
      if(return_lS){
        lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
        for(paralcv in 1:CGGP$numpara){
          dlS[levellcv, CGGP$numpara*(dimlcv-1)+paralcv] = -sum(diag(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv]))
        }
      }
    }
  }
  
  if(is.matrix(y)){
    numout = dim(y)[2]
    sigma2 = rep(0,numout) # Predictive weight for each measured point
    dsigma2 = matrix(0,nrow=CGGP$numpara*CGGP$d,ncol=numout) # Predictive weight for each measured point
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        VVV1=unlist(cholS[gg+CGGP$uo[blocklcv,]])
        VVV2=unlist(dMatdtheta[gg+CGGP$uo[blocklcv,]])
        VVV3=CGGP$gridsizest[blocklcv,]
        for(outdimlcv in 1:numout){
          B2 = y[IS,outdimlcv]
          B0 = revc[IS,outdimlcv]
          B = (CGGP$w[blocklcv])*B0#/dim(y)[1]
          dB = rcpp_gkronDBS(VVV1,VVV2,B,VVV3)
          dsigma2[,outdimlcv] = dsigma2[,outdimlcv] + as.vector(dB%*%B2)
          sigma2[outdimlcv] = sigma2[outdimlcv]  + sum(B2*B)
        }
      }
    }
    out <- list(sigma2=sigma2,
                dsigma2=dsigma2)
    if (return_lS) {
      out$lS <- lS
      out$dlS <- dlS
    }
  }else{
    sigma2 = 0 # Predictive weight for each measured point
    dsigma2 = rep(0,nrow=CGGP$d) # Predictive weight for each measured point
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        B0 = revc[IS]
        B2 = y[IS]
        B = (CGGP$w[blocklcv])*B0#/length(y)
        dB = rcpp_gkronDBS(unlist(cholS[gg+CGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
        
        dsigma2 = dsigma2 +t(B2)%*%t(dB)
        sigma2 = sigma2 + t(B2)%*%B
      }
    }
    out <- list(sigma2=sigma2,
                dsigma2=dsigma2)
    if (return_lS) {
      out$lS <- lS
      out$dlS <- dlS
    }
  }
  out
}
