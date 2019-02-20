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
SGGP_internal_postvarmatcalcfaster <- function(GMat, dGMat,cholS,dSMat,INDSN,numpara,returndG = FALSE,returndiag = FALSE) {
  
  CoinvC1o = backsolve(cholS,backsolve(cholS,t(GMat[,INDSN]), transpose = TRUE))
  if(returndiag){
    Sigma_mat = rowSums(t((CoinvC1o))*((GMat[,INDSN])))
  }else{
    Sigma_mat = (t(CoinvC1o))%*%(t(GMat[,INDSN]))
  }
  
  if(returndG){
    if(returndiag){
      dSigma_mat = matrix(0,dim(GMat)[1],numpara)
    }else{
      dSigma_mat = matrix(0,dim(GMat)[1],dim(GMat)[1]*numpara)
    }
    np = dim(GMat)[1]
    nb = dim((cholS))[1]
    for(k in 1:numpara){
      dS = dSMat[,(k-1)*nb+(1:nb)]
      
      CoinvC1oE = ((as.matrix(dS))%*%t(GMat[,INDSN]))
      dCoinvC1o = backsolve(cholS,backsolve(cholS,CoinvC1oE, transpose = TRUE))
      dCoinvC1o = dCoinvC1o + backsolve(cholS,backsolve(cholS,t(dGMat[,(k-1)*nb+INDSN]), transpose = TRUE))
      if(returndiag){
        dSigma_mat[,k] = rowSums(t(dCoinvC1o)*(GMat[,INDSN]))+rowSums(t(CoinvC1o)*(dGMat[,(k-1)*nb+INDSN]))
      }else{
        dSigma_mat[,((k-1)*np+1):(k*np)] =t(dCoinvC1o)%*%t(GMat[,INDSN])+t(CoinvC1o)%*%t(dGMat[,(k-1)*nb+INDSN])
      }
    }
    return(list("Sigma_mat"= Sigma_mat,"dSigma_mat" = dSigma_mat))
  }else{
    return(Sigma_mat)
  }
}