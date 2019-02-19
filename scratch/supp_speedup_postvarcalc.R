
#' Calculate posterior variance
#'
#' @param x1 Points at which to calculate MSE
#' @param x2 Levels along dimension, vector???
#' @param xo No clue what this is
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#' for vector of 1D points.
#' @param ... Placeholder
#' @param returndPVMC Should dPVMC be returned?
#' @param returndiagonly Should only the diagonal be returned?
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




#' Calculate posterior variance
#'
#' @param x1 Points at which to calculate MSE
#' @param x2 Levels along dimension, vector???
#' @param xo No clue what this is
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#' for vector of 1D points.
#' @param ... Placeholder
#' @param returndPVMC Should dPVMC be returned?
#' @param returndiagonly Should only the diagonal be returned?
#'
#' @return Variance posterior
#' @export
#'
#' @examples
#' SGGP_internal_postvarmatcalc(c(.4,.52), c(0,.25,.5,.75,1),
#'              xo=c(.11), theta=c(.1,.2,.3),
#'              CorrMat=SGGP_internal_CorrMatCauchySQT)
SGGP_internal_postvarmatcalcfaster <- function(GMat, dGMat,Cmat,INDSN,numpara ) {
  cholS = chol(Cmat$C[INDSN,INDSN])
  
  CoinvC1o = backsolve(cholS,backsolve(cholS,t(GMat[,INDSN]), transpose = TRUE))
  Sigma_mat = -((t(CoinvC1o))%*%(t(GMat[,INDSN])))
  dSigma_mat = matrix(0,dim(GMat)[1],dim(GMat)[1]*numpara)
  np = dim(GMat)[1]
  nb = dim((Cmat$C))[1]
  for(k in 1:numpara){
    dS = Cmat$dCdtheta[INDSN,(k-1)*nb+INDSN]
    CoinvC1oE = as.matrix(dS)%*%CoinvC1o
    dCoinvC1o = -backsolve(cholS,backsolve(cholS,CoinvC1oE, transpose = TRUE))
    dCoinvC1o = dCoinvC1o + backsolve(cholS,backsolve(cholS,t(dGMat[,(k-1)*nb+INDSN]), transpose = TRUE))
    dSigma_mat[,((k-1)*np+1):(k*np)] =-t(dCoinvC1o)%*%t(GMat[,INDSN])-t(CoinvC1o)%*%t(dGMat[,(k-1)*nb+INDSN])
  }
  return(list("Sigma_mat"= Sigma_mat,"dSigma_mat" = dSigma_mat))
}