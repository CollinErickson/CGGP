MSEpred_calc <- function(xp,xl, theta) {
  S = CorrMat(xl, xl, theta)
  Schol = chol(S)
  n = length(xl)
  Cp = CorrMat(xp, xl, theta)
  
  Cv = backsolve(Schol,backsolve(Schol,t(Cp), transpose = TRUE))
  
  MSE_val = diag_corrMat(xp, theta)-rowSums(t(Cv)*Cp)
  return(MSE_val)
}

SGGPpred <- function(xp,SG, y,theta) {
    my = mean(y)
    y = y-my
    
    Q  = max(SG$uo[1:SG$uoCOUNT,])
    CiS = list(matrix(1,1,1),Q*SG$d)
    CS = list(matrix(1,1,1),Q*SG$d)
    CCS = list(matrix(1,1,1),Q*SG$d)
    for (lcv2 in 1:SG$d) {
      for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
        Xbrn = SG$xb[1:SG$sizest[lcv1]]
        Xbrn = Xbrn[order(Xbrn)]
        S = CorrMat(Xbrn, Xbrn , theta[lcv2])
        CS[[(lcv2-1)*Q+lcv1]] = S
        CCS[[(lcv2-1)*Q+lcv1]] = chol(S)
      }
    }
    
    pw = rep(0, length(y))
    
    for (lcv1 in 1:SG$uoCOUNT) {
      
      B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
      for (e in SG$d:1) {
        if(SG$gridsizest[lcv1,e] > 1.5){
          B <- matrix(as.vector(B),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
          B <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B, transpose = TRUE))
          #B <-  solve(CS[[((e-1)*Q+SG$uo[lcv1,e])]],B)
          B <- t(B)
        }
        else{
          B = as.vector(B)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
        }
      }
      
      pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
        SG$w[lcv1] * B
    }
    sigma_hat = t(y) %*% pw / length(y)
    
    Cp = matrix(1,dim(xp)[1],SG$ss)
    for (e in 1:SG$d) {
      V = CorrMat(xp[,e], SG$xb, theta[e])
      Cp = Cp*V[,SG$designindex[,e]]
      #Cp = Cp*CorrMat(xp[,e], SG$design[,e], theta[e])
    }
#    Cact = matrix(1,SG$ss,SG$ss)
#    for (e in 1:SG$d) {
#      Cact = Cact*CorrMat(SG$design[,e], SG$design[,e], theta[e])
#    }
    
    
    
    MSE_v = array(0, c(SG$d, 9,dim(xp)[1]))
    for (lcv1 in 1:SG$d) {
      MSE_v[lcv1, 1,] = diag_corrMat(xp[,lcv1], theta[lcv1])
    }
    for (lcv1 in 1:SG$d) {
      for (lcv2 in 1:8) {
        MSE_v[lcv1, lcv2+1,] = MSEpred_calc(xp[,lcv1],SG$xb[1:SG$sizest[lcv2]], theta[lcv1])
        MSE_v[lcv1, lcv2+1,] = pmin(MSE_v[lcv1, lcv2+1,], MSE_v[lcv1, lcv2,])
      }
    }
    
    ME_t = prod(MSE_v[,1,],1)
    for (lcv1 in 1:SG$uoCOUNT) {
      ME_v = rep(1,dim(xp)[1])
      for (e in 1:SG$d) {
        levelnow = SG$uo[lcv1,e]
          ME_v = ME_v*(MSE_v[e,1,]-MSE_v[e,levelnow+1,])
      }
      ME_t = ME_t-SG$w[lcv1]*ME_v
    }
    GP = list("mean" = (my+Cp %*% pw), "var"=sigma_hat[1]*ME_t)
    
    
    return(GP)
}



