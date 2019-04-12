CGGP_internal_imputesomegrid <- function(CGGP, y) {
  if(!is.matrix(y)){
    numoutputs = 1
    yimputed <- y
  }else{
    numoutputs = dim(y)[2]
    yimputed <- y
  }
  
  for(oplcv in 1:numoutputs){
    if(!is.matrix(y)){
      y.thisloop = y
      cholS.thisloop = CGGP$cholSs
      thetaMAP.thisloop = CGGP$thetaMAP
    }else{
      y.thisloop = as.vector(y[,oplcv])
      
      if(is.matrix(CGGP$thetaMAP)){
        thetaMAP.thisloop = CGGP$thetaMAP[,oplcv]
        cholS.thisloop = CGGP$cholSs[[oplcv]]
      }else{
        thetaMAP.thisloop = CGGP$thetaMAP
        cholS.thisloop = CGGP$cholSs
      }
    }
    if(any(is.na(y.thisloop))){
      Is = sort(which(is.na(y.thisloop)))
      
      xp = as.matrix(CGGP$design[Is,])
      if(dim(xp)[2]<CGGP$d){
        xp = t(xp)
      }
      n2pred = dim(xp)[1]
      
      brokenblocks = unique(CGGP$blockassign[Is])
      possblocks = 1:max(CGGP$uoCOUNT,20)
      
      
      for (lcv in 1:n2pred){
        Bs = CGGP$blockassign[Is[lcv]]
        possblocks = unique(c(possblocks,CGGP$uala[Bs,1:sum(CGGP$uo[Bs,]>1.5)]))
      }
      possblocks = sort(possblocks)
      keepthisone = rep(1,length(possblocks))
      
      for (lcv in 1:length(possblocks)){
        if (any(Is %in% CGGP$dit[possblocks[lcv], 1:CGGP$gridsizet[possblocks[lcv]]])){
          keepthisone[lcv] = 0;
        }
      }
      possblocks = possblocks[which(keepthisone>0.5)]
      
      # Cp is sigma(x_0) in paper, correlation vector between design points and xp
      Cp = matrix(0,n2pred,CGGP$ss)
      GGGG = list(matrix(1,n2pred,length(CGGP$xb)),CGGP$d)
      for (dimlcv in 1:CGGP$d) { # Loop over dimensions
        V = CGGP$CorrMat(xp[,dimlcv], CGGP$xb[1:CGGP$sizest[max(CGGP$uo[,dimlcv])]],
                         thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                         returnlogs=TRUE)
        GGGG[[dimlcv]] = exp(V)
        Cp = Cp+V[,CGGP$designindex[,dimlcv]]
      }
      Cp = exp(Cp)
      
      ME_t = matrix(1,n2pred,1)
      MSE_v = list(matrix(0,n2pred,2),(CGGP$d+1)*(CGGP$maxlevel+1)) 
      Q  = max(CGGP$uo[1:CGGP$uoCOUNT,])
      for (dimlcv in 1:CGGP$d) {
        for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
          Q  = max(CGGP$uo[1:CGGP$uoCOUNT,])
          gg = (dimlcv-1)*Q
          INDSN = 1:CGGP$sizest[levellcv]
          INDSN = INDSN[sort(CGGP$xb[1:CGGP$sizest[levellcv]],
                             index.return = TRUE)$ix]
          MSE_v[[(dimlcv)*CGGP$maxlevel+levellcv]] =
            CGGP_internal_postvarmatcalc_fromGMat(GGGG[[dimlcv]],
                                                  c(),
                                                  as.matrix(
                                                    cholS.thisloop[[gg+levellcv]]
                                                  ),
                                                  c(),
                                                  INDSN,
                                                  CGGP$numpara,
                                                  returndiag=TRUE)
        }
      }
      
      ME_t = matrix(1,n2pred,length(possblocks))
      for (blocklcv in 1:length(possblocks)) {
        ME_s = matrix(1,n2pred,1)
        for (dimlcv in 1:CGGP$d) {
          levelnow = CGGP$uo[possblocks[blocklcv],dimlcv]
          ME_s = ME_s*as.matrix(MSE_v[[(dimlcv)*CGGP$maxlevel+levelnow]])
        }
        ME_t[,blocklcv] = as.vector(ME_s)
      }
      
      Jstar = possblocks[apply(ME_t, 1, which.max)]
      w = rep(0,n2pred)
      w = 1-apply(ME_t, 1, max)
      
      #find x0
      yn0 = y.thisloop
      yhat0 = rep(0,n2pred)
      for(lcv in 1:length(Jstar)){
        Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max value of all blocks
        gg = (1:CGGP$d-1)*Q
        IS = CGGP$dit[Jstar[lcv], 1:CGGP$gridsizet[Jstar[lcv]]];
        B = y.thisloop[IS]
        rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[Jstar[lcv],]]), B, CGGP$gridsizest[Jstar[lcv],])
        yhat0[lcv] = Cp[lcv,IS]%*%B
        yn0[Is[lcv]] = yhat0[lcv]
      }
      
      #find g at x0
      pwforg = rep(0,length(y.thisloop))
      for (blocklcv in 1:CGGP$uoCOUNT) {
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
        if(any(Is %in% IS) ){
          B = yn0[IS]
          rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
          pwforg[IS] = pwforg[IS]+CGGP$w[blocklcv] * B
        }
      }
      pw0 = pwforg
      pwforg = pwforg[Is]
      dir0 = pwforg
      
      wgst1 = rep(0,length(y.thisloop))
      wgst2 = rep(0,length(y.thisloop))
      wgst1[Is] = w*dir0
      dotheseblocks = rep(0,CGGP$uoCOUNT)
      for (blocklcv in 1:CGGP$uoCOUNT) {
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
        B = wgst1[IS]
        rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
        wgst2[IS] = wgst2[IS]+CGGP$w[blocklcv] * B
      }
      
      lambdas = pmax(pmin(sum(wgst2*yn0)/sum(wgst1[Is]*wgst2[Is]),1.1),-0.1)
      yhat1 = yhat0-lambdas*w*dir0;
      yn1 = yn0
      yn1[Is] = yhat1
      
      pwforg = rep(0,length(y.thisloop))
      for (blocklcv in 1:CGGP$uoCOUNT) {
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
        if(any(Is %in% IS) ){
          B = yn1[IS]
          rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
          pwforg[IS] = pwforg[IS]+CGGP$w[blocklcv] * B
        }
      }
      pw1 = pwforg
      pwforg = pwforg[Is]
      dir1 = pwforg
      
      M = 12
      s = matrix(0,length(dir1),M)
      ny = matrix(0,length(dir1),M)
      dirsave = matrix(0,length(dir1),M)
      xsave = matrix(0,length(dir1),M)
      rho = rep(0,M)
      gamma = rep(0,M)
      
      #x0 = yhat0, x1 = yhat1
      #g(x0)= dir0, x1 = dir1
      lcv = 1
      s[,lcv] = yhat1-yhat0
      ny[,lcv] = dir1-dir0
      xsave[,lcv] = yhat1
      dirsave[,lcv] = dir1
      
      L = rep(0,40)
      for(lcv in 1:40){
        
        q = dir1
        if(lcv > 1.5){
          for(k in 2:min(lcv,M)){#q(min(lcv,M),1,by=-1)
            rho[k] = 1/sum(s[,k]*ny[,k])
            gamma[k] = rho[k]*sum(s[,k]*q)
            q = q-gamma[k]*ny[,k]
          }
        }
        
        r = q*mean(s[,1]*ny[,1])/mean(ny[,1]*ny[,1])
        if(lcv > 1.5){
          for(k in min(lcv,M):2){#q(min(lcv,M),1,by=-1)
            beta = rho[k]*sum(ny[,k]*r)
            r =r+(gamma[k]-beta)*s[,k]
          }
        }
        dir1u = -r
        
        wgst1 = rep(0,length(y.thisloop))
        wgst2 = rep(0,length(y.thisloop))
        wgst1[Is] = dir1u
        dotheseblocks = rep(0,CGGP$uoCOUNT)
        for (blocklcv in 1:CGGP$uoCOUNT) {
          IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
          if(any(Is %in% IS) ){
            B = wgst1[IS]
            rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
            wgst2[IS] = wgst2[IS]+CGGP$w[blocklcv] * B
          }
        }
        
        lambdas =-sum(wgst2*yn1)/sum(wgst1[Is]*wgst2[Is])
        
        for(k in  (min(lcv,M-1)):1){
          s[,k+1] = s[,k]
          ny[,k+1] = ny[,k]
          xsave[,k+1] = xsave[,k]
          dirsave[,k+1] = dirsave[,k]
        }
        
        yhat2 = yhat1+lambdas*dir1u
        
        yn2 = yn1
        yn2[Is] = yhat2
        
        pwforg = rep(0,length(y.thisloop))
        for (blocklcv in 1:CGGP$uoCOUNT) {
          IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
          if(any(Is %in% IS) ){
            B = yn2[IS]
            rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
            pwforg[IS] = pwforg[IS]+CGGP$w[blocklcv] * B
          }
        }
        pw2 = pwforg
        pwforg = pwforg[Is]
        dir2 = pwforg
        
        s[,1] = yhat2-xsave[,1]
        ny[,1] = dir2-dirsave[,2]
        xsave[,1] = yhat2
        dirsave[,1] = dir2
        
        dir1 = dir2
        yhat1 = yhat2
        yn1 = yn2
        pw1 = pw2
        
        L[lcv] = sum( pw2*yn2)
        if(lcv > (M+1)){
          if(max(abs(L[lcv:(lcv-3)]-L[(lcv-1):(lcv-4)])) < 10^(-3)*L[lcv]){
            break
          }
        }
        if(lcv > 5 && lcv <= M+1){
          if(max(abs(L[lcv:(lcv-3)]-L[(lcv-1):(lcv-4)])) < 10^(-3)*L[lcv]){
            break
          }
        }
      }
      
      
      if(!is.matrix(y)){
        yimputed <- yn1
      }else{
        yimputed[,oplcv] <- yn1
      }
    }
  }
  return(yimputed)
}


#' Calculate MSE prediction along a single dimension
#'
#' @param xp Points at which to calculate MSE
#' @param xl Levels along dimension, vector???
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#'
#' @return MSE predictions
#' @export
#'
#' @examples
#' CGGP_internal_MSEpredcalc(c(.4,.52), c(0,.25,.5,.75,1), theta=c(.1,.2),
#'              CorrMat=CGGP_internal_CorrMatCauchySQ)
CGGP_internal_MSEpredcalc <- function(xp,xl,theta,CorrMat) {
  S = CorrMat(xl, xl, theta)
  n = length(xl)
  cholS = chol(S)
  
  Cp = CorrMat(xp, xl, theta)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  
  MSE_val = 1 - rowSums(t(CiCp)*((Cp)))
  return(MSE_val)
}



#' Calculate posterior variance, faster version
#'
#' @param GMat1 Matrix 1
#' @param GMat2 Matrix 2
#' @param cholS Cholesky factorization of S
#' @param INDSN Indices, maybe
#'
#' @return Variance posterior
## @export
#' @noRd
CGGP_internal_postvarmatcalc_fromGMat_asym <- function(GMat1,GMat2,cholS,INDSN) {
  CoinvC1o = backsolve(cholS,
                       backsolve(cholS,t(GMat1[,INDSN]), transpose = TRUE))
  Sigma_mat = (t(CoinvC1o)%*%(t(GMat2[,INDSN])))
  
  return(Sigma_mat)
  
}

#' Calculate posterior variance, faster version
#'
#' @param GMat Matrix
#' @param dGMat Derivative of matrix
#' @param cholS Cholesky factorization of S
#' @param dSMat Deriv of SMat
#' @param INDSN Indices, maybe
#' @param numpara Number of parameters for correlation function
#' @param returnlogs Should log scale be returned
#' @param returnderiratio Should derivative ratio be returned?
#' @param returndG Should dG be returned
#' @param returndiag Should diag be returned
#' @param ... Placeholder
#'
#' @return Variance posterior
# @export
#' @noRd
CGGP_internal_postvarmatcalc_fromGMat <- function(GMat, dGMat,cholS,dSMat,INDSN,numpara,...,
                                                  returnlogs=FALSE, returnderiratio =FALSE,
                                                  returndG = FALSE,returndiag = FALSE) {
  
  # Next line was giving error with single value, so I changed it
  CoinvC1o = backsolve(cholS,backsolve(cholS, t(GMat[,INDSN, drop=F]), transpose = TRUE))
  # backsolve1 <- backsolve(cholS,if (is.matrix(GMat[,INDSN]) && ncol(GMat[,INDSN])>1) t(GMat[,INDSN]) else GMat[,INDSN], transpose = TRUE)
  # CoinvC1o = backsolve(cholS,backsolve1)
  if(returndiag){
    if(!returnlogs){
      Sigma_mat = rowSums(t((CoinvC1o))*((GMat[,INDSN])))
    }else{
      nlSm = rowSums(t((CoinvC1o))*((GMat[,INDSN])))
      Sigma_mat =log(nlSm)
    }
    np = length(Sigma_mat)
  }else{
    if(!returnlogs){
      Sigma_mat = t(CoinvC1o)%*%(t(GMat[,INDSN]))
    }else{
      nlSm =(t(CoinvC1o))%*%(t(GMat[,INDSN]))
      Sigma_mat =log(nlSm)
    }
    np = dim(Sigma_mat)[1]
  }
  
  
  if(returndG){
    nb = dim(GMat)[2]
    nc = dim(dSMat)[1]
    if(returndiag){
      dSigma_mat = matrix(0,np,numpara)
    }else{
      dSigma_mat = matrix(0,np,np*numpara)
    }
    
    for(k in 1:numpara){
      dS = dSMat[,(k-1)*nc+(1:nc)]
      CoinvC1oE = ((as.matrix(dS))%*%t(GMat[,INDSN]))
      if(returndiag){
        dCoinvC1o = backsolve(cholS,backsolve(cholS,CoinvC1oE, transpose = TRUE))
        dCoinvC1o = dCoinvC1o + backsolve(cholS,backsolve(cholS,t(dGMat[,(k-1)*nb+INDSN]), transpose = TRUE))
        
        if(!returnlogs && !returnderiratio){
          dSigma_mat[,k] = rowSums(t(dCoinvC1o)*(GMat[,INDSN]))+rowSums(t(CoinvC1o)*(dGMat[,(k-1)*nb+INDSN]))
        }else if(returnderiratio){
          dSigma_mat[,k] = rowSums(t(dCoinvC1o)*(GMat[,INDSN]))+rowSums(t(CoinvC1o)*(dGMat[,(k-1)*nb+INDSN]))/Sigma_mat
        }else{
          dSigma_mat[,k] = (rowSums(t(dCoinvC1o)*(GMat[,INDSN]))+rowSums(t(CoinvC1o)*(dGMat[,(k-1)*nb+INDSN])))/nlSm
        }
      }else{
        dCoinvC1_part1 = t(CoinvC1o)%*%(CoinvC1oE)
        dCoinvC1_part2 = t(CoinvC1o)%*%t(dGMat[,(k-1)*nb+INDSN])
        dCoinvC1_part2 = dCoinvC1_part2+t(dCoinvC1_part2)
        if(!returnlogs && !returnderiratio){
          dSigma_mat[,(k-1)*np + 1:np] =dCoinvC1_part1+dCoinvC1_part2
        }else if(returnderiratio){
          dSigma_mat[,(k-1)*np + 1:np] =( dCoinvC1_part1+dCoinvC1_part2)/Sigma_mat
        }else{
          dSigma_mat[,(k-1)*np + 1:np] =( dCoinvC1_part1+dCoinvC1_part2)/nlSm
        }
      }
    }
    return(list("Sigma_mat"= Sigma_mat,"dSigma_mat" = dSigma_mat))
  }else{
    return(Sigma_mat)
  }
}
