#' ????????????
#'
#' @param xp 
#' @param xl Levels along dimension, vector???
#' @param theta Correlation parameters
#'
#' @return
#' @export
#'
#' @examples
MSEpred_calc <- function(xp,xl, theta) {
  S = CorrMat(xl, xl, theta)
  t = exp(theta)
  n = length(xl)
  Ci = solve(S)
  
  Cp = CorrMat(xp, xl, theta)
  
  MSE_val = 1-rowSums((Cp%*%Ci)*((Cp)))
  return(MSE_val)
  
}



#' ????????????
#' 
#' Predict using SG with y values at xp?
#' Shouldn't y values already be stored in SG?
#'
#' @param xp 
#' @param SG 
#' @param y 
#' @param theta 
#'
#' @return
#' @export
#'
#' @examples
SGGPpred <- function(xp,SG, y,theta) {
  # Center outputs
  my = mean(y)
  y = y-my
  
  Q  = max(SG$uo[1:SG$uoCOUNT,])
  CiS = list(matrix(1,1,1),Q*SG$d)
  for (lcv2 in 1:SG$d) {
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      S = CorrMat(Xbrn, Xbrn , theta[lcv2])
      CiS[[(lcv2-1)*Q+lcv1]] = solve(S)
    }
  }
  
  
  pw = rep(0, length(y))
  for (lcv1 in 1:SG$uoCOUNT) {
    narrowd = which(SG$uo[lcv1,] > 1.5)
    Ci = 1
    if (lcv1 == 1) {
      narrowd = 1
    } else{
      narrowd = which(SG$uo[lcv1,] > 1.5)
    }
    for (e in narrowd) {
      Ci = kronecker(Ci, CiS[[((e-1)*Q+SG$uo[lcv1,e])]])
    }
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * Ci %*% y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
  }
  sigma_hat = t(y) %*% pw / length(y)
  
  Cp = matrix(1,dim(xp)[1],SG$ss)
  for (e in 1:SG$d) {
    Cp = Cp*CorrMat(xp[,e], SG$design[,e], theta[e])
  }
  
  #    Cact = matrix(1,SG$ss,SG$ss)
  #    for (e in 1:SG$d) {
  #      Cact = Cact*CorrMat(SG$design[,e], SG$design[,e], theta[e])
  #    }
  
  
  MSE_v = array(0, c(SG$d, 8,dim(xp)[1]))
  for (lcv1 in 1:SG$d) {
    for (lcv2 in 1:8) {
      MSE_v[lcv1, lcv2,] = abs(MSEpred_calc(xp[,lcv1],SG$xb[1:SG$sizest[lcv2]], theta[lcv1]))
      if (lcv2 > 1.5) {
        MSE_v[lcv1, lcv2,] = pmin(MSE_v[lcv1, lcv2,], MSE_v[lcv1, lcv2 - 1,])
      }
    }
  }
  
  
  ME_t = rep(1,dim(xp)[1])
  for (lcv1 in 1:SG$uoCOUNT) {
    ME_v = rep(1,dim(xp)[1])
    for (e in 1:SG$d) {
      levelnow = SG$uo[lcv1,e]
      if(levelnow > 1.5){
        ME_v = ME_v*(1-MSE_v[e,levelnow,])
      } else {
        ME_v = ME_v*(1-MSE_v[e,levelnow,])
      }
    }
    ME_t = ME_t-SG$w[lcv1]*ME_v
  }
  
  GP = list("mean" = (my+Cp %*% pw), "var"=sigma_hat[1]*ME_t)
  
  
  return(GP)
}



