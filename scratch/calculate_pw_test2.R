cppFunction("
            NumericVector pw_c(NumericMatrix CholC, NumericVector UO, NumericMatrix DIT, NumericVector SizeGrid, NumericVector SizeGrids, NumericVector Ydata,  NumericVector w, NumericVector MatrixRefNumbers, int n, int nblock, int d, int maxLevel, int maxGridSizes){
            
            NumericVector currLevels(d);
            NumericVector x(maxGridSizes);
            NumericVector y(maxGridSizes);
            NumericVector b(maxGridSizes);

            NumericVector pw(n);
            int N = 0;

          
            
            for(int iblock = 0; iblock<nblock;iblock++){
            for(int dim = d-1; dim>=0;dim--) currLevels(dim) = UO(iblock,dim);
            N = SizeGrid(iblock);
            for(int nstart = 0; nstart<N;nstart++) b(nstart) = Ydata(DIT(iblock,nstart)-1);

            for(int dim = d-1; dim>=0;dim--){
            int p0 =  SizeGrids(iblock,dim);
            int n0 = N/p0;
            int sv = MatrixRefNumbers(currLevels(dim)-1);
            
            for(int h = 0; h < n0; h++)
            {
            for ( int i = 0; i < p0; i++)
            {
            x(h*p0+i) = b(h*p0+i);
            for ( int k = i-1; k >=0; k-- ) x(h*p0+i) -= CholC(sv+i*p0+k,dim) * x(h*p0+k);
            x(h*p0+i) /= CholC(sv+i*p0+i,dim);
            }
            
            for ( int i = p0-1; i >=0; i--)
            {
            y(h*p0+i) = x(h*p0+i);
            for ( int k=i+1; k < p0; k++) y(h*p0+i) -= CholC(sv+k*p0+i,dim)  *y(h*p0+k);
            y(h*p0+i) /= CholC(sv+i*p0+i,dim);
            }
            }
            
            int c = 0;
            for ( int i = 0; i < p0; i++)
            {
            for(int h = 0; h < n0; h++)
            {
            b(c) = y(h*p0+i);
            c++;
            }
            }
            }

            for(int nstart = 0; nstart<N;nstart++) {
            pw(DIT(iblock,nstart)-1) += w(iblock)*b(nstart);
            }
            }
            
            return pw;
            } ")



#' Calculate predictive weights for SGGP
#' 
#' Predictive weights are Sigma^{-1}*y in standard GP.
#' This calculation is much faster since we don't need to
#' solve the full system of equations.
#'
#' @param SG SGGP object
#' @param y Measured values for SG$design
#' @param logtheta Log of correlation parameters
#' @param return_lS Should lS be returned?
#'
#' @return Vector with predictive weights
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' calculate_pw(SG=SG, y=y, logtheta=c(-.1,.1,.3))
calculate_pw2 <- function(SG, y, logtheta, return_lS=FALSE) {
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max value of all blocks
  # Now going to store choleskys instead of inverses for stability
  #CiS = list(matrix(1,1,1),Q*SG$d) # A list of matrices, Q for each dimension
  CCS = list(matrix(1,1,1),Q*SG$d)
  lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Save log determinant of matrices
  # Loop over each dimension
  
  maxSGlevel = max(SG$uo[1:SG$uoCOUNT,])
  
  SGlevelMatSizes = (SG$sizest[1:max(SG$uo[1:SG$uoCOUNT,])])^2
  sumCholSave = sum(SGlevelMatSizes)
  cumsumCholSave = c(0,cumsum(SGlevelMatSizes))
  
  
  CholSave = matrix(0,nrow = sum(SGlevelMatSizes),ncol = d)
  for (lcv2 in 1:SG$d) {
    # Loop over each possible needed correlation matrix
    for (lcv1 in 1:maxSGlevel) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]] # xb are the possible points
      Xbrn = Xbrn[order(Xbrn)] # Sort them low to high, is this necessary? Probably just needs to be consistent.
      S = SG$CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      diag(S) = diag(S) + SG$nugget
      
      KeepMat = length(S)
      # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
      solvetry <- try({
        #CiS[[(lcv2-1)*Q+lcv1]] = solve(S)
        CholSave[cumsumCholSave[lcv1]+(1:SGlevelMatSizes[lcv1]),lcv2] = as.vector((chol(S)))
      })
      if (inherits(solvetry, "try-error")) {return(Inf)}
    }
  }
  
  
  pw = rep(0, length(y)) # Predictive weight for each measured point
  # Loop over blocks selected
  
  pw = pw_c(CholSave,SG$uo,SG$dit,SG$gridsizet,SG$gridsizest,y,SG$w,cumsumCholSave,length(y),SG$uoCOUNT,d,maxSGlevel,max(SG$gridsizet))
  
  if (return_lS) {
    return(list(pw=pw, lS=lS))
  }
  pw
}

#' Calculate derivative of pw
#'
#' @inheritParams calculate_pw
#' @param return_dlS Should dlS be returned?
#'
#' @return derivative matrix of pw with respect to logtheta
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' calculate_pw_and_dpw(SG=SG, y=y, logtheta=c(-.1,.1,.3))
calculate_pw_and_dpw <- function(SG, y, logtheta, return_lS=FALSE, return_dlS=FALSE) {
  
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max level of all blocks
  # Now storing choleskys instead of inverses
  CCS = list(matrix(1,1,1),Q*SG$d) # To store choleskys
  dCS = list(matrix(1,1,1),Q*SG$d)
  #CiS = list(matrix(1,1,1),Q*SG$d) # To store correlation matrices
  #dCiS = list(matrix(1,1,1),Q*SG$d) # To store derivatives of corr mats
  lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store log det S
  dlS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # Store deriv of log det S
  dlS2 = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d) # ???? added for chols
  
  # Loop over each dimension
  for (lcv2 in 1:SG$d) {
    # Loop over depth of each dim
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      S = SG$CorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      diag(S) = diag(S) + SG$nugget
      dS = SG$dCorrMat(Xbrn, Xbrn , logtheta=logtheta[lcv2])
      
      CCS[[(lcv2-1)*Q+lcv1]] = chol(S)#-CiS[[(lcv2-1)*Q+lcv1]]  %*% 
      dCS[[(lcv2-1)*Q+lcv1]] = -dS
      lS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
      V = solve(S,dS);
      dlS[lcv1, lcv2] = sum(diag(V))
    }
  }
  
  pw = rep(0, length(y)) # predictive weights
  
  dpw = matrix(0, nrow = length(y), ncol = SG$d) # derivative of predictive weights
  
  for (lcv1 in 1:SG$uoCOUNT) {
    B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],] = dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],] +
        SG$w[lcv1] * gkronDBS(unlist(CCS[((1:d-1)*Q+SG$uo[lcv1,1:d])]), unlist(dCS[((1:d-1)*Q+SG$uo[lcv1,1:d])]), B, SG$gridsizest[lcv1,], length(B2),  d)
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * B
  }
  
 # print(dpw[1:10,])
  
  out <- list(pw=pw,
              dpw=dpw)
  if (return_lS) {
    out$lS <- lS
  }
  if (return_dlS) {
    out$dlS <- dlS
  }
  out
}
