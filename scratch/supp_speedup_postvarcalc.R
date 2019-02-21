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
SGGP_internal_postvarmatcalcfaster <- function(GMat, dGMat,cholS,dSMat,INDSN,numpara,...,returnlogs=FALSE, returnderiratio =FALSE,returndG = FALSE,returndiag = FALSE) {
  
  CoinvC1o = backsolve(cholS,backsolve(cholS,t(GMat[,INDSN]), transpose = TRUE))
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