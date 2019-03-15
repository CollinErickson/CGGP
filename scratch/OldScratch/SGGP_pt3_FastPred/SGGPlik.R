
lik <- function(x, SG, y) {
  theta = x
  
  
  if (max(theta)<3) {
    
    Q  = max(SG$uo[1:SG$uoCOUNT,])
    CS = list(matrix(1,1,1),Q*SG$d)
    CCS = list(matrix(1,1,1),Q*SG$d)
    lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d)
    for (lcv2 in 1:SG$d) {
      for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
        Xbrn = SG$xb[1:SG$sizest[lcv1]]
        Xbrn = Xbrn[order(Xbrn)]
        S= CorrMat(Xbrn, Xbrn , theta[lcv2])
        CCS[[(lcv2-1)*Q+lcv1]] = chol(S)
        lS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
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
    
    
    lDet = 0
    
    for (lcv1 in 1:SG$uoCOUNT) {
      for (lcv2 in 1:SG$d) {
        levelnow = SG$uo[lcv1, lcv2]
        if (levelnow > 1.5) {
          lDet = lDet + (lS[levelnow, lcv2] - lS[levelnow - 1, lcv2]) * (SG$gridsize[lcv1]) /
            (SG$gridsizes[lcv1, lcv2])
        }
      }
    }
    return(log(sigma_hat)+sum(theta^2)/length(y)+1 / length(y) * lDet )
  }else{
    return(Inf)
  }
  
}
glik <- function(x, SG, y) {
  theta = x
  
  
  Q  = max(SG$uo[1:SG$uoCOUNT,])
  CCS = list(matrix(1,1,1),Q*SG$d)
  CiS = list(matrix(1,1,1),Q*SG$d)
  dCS = list(matrix(1,1,1),Q*SG$d)
  lS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d)
  dlS = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d)
  dlS2 = matrix(0, nrow = max(SG$uo[1:SG$uoCOUNT,]), ncol = SG$d)
  for (lcv2 in 1:SG$d) {
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]]
      Xbrn = Xbrn[order(Xbrn)]
      S = CorrMat(Xbrn, Xbrn , theta[lcv2])
      dS = dCorrMat(Xbrn, Xbrn , theta[lcv2])
      CCS[[(lcv2-1)*Q+lcv1]] = chol(S)#-CiS[[(lcv2-1)*Q+lcv1]]  %*% 
      dCS[[(lcv2-1)*Q+lcv1]] = dS
      lS[lcv1, lcv2] = 2*sum(log(diag(CCS[[(lcv2-1)*Q+lcv1]])))
      V = solve(S,dS);
      dlS[lcv1, lcv2] = sum(diag(V))
    }
  }
  
  pw = rep(0, length(y))
  
  dpw = matrix(0, nrow = length(y), ncol = SG$d)
  
  
  for (lcv1 in 1:SG$uoCOUNT) {
    
    B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    for (e in SG$d:1) {
      if(SG$gridsizest[lcv1,e] > 1.5){
        B <- matrix(as.vector(B),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
        B <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B, transpose = TRUE))
        B <- t(B)
      }
      else{
        B = as.vector(B)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
      }
    }
    
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * B
    
    
    B3 = B
    for (e in  SG$d:1) {
      if(SG$gridsizest[lcv1,e] > 1.5){
        B3 <- matrix(as.vector(B3),SG$gridsizest[lcv1,e],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e])
      }
      
      B2 = B3
      
      if(SG$gridsizest[lcv1,e] > 1.5){
        B3 <- t(B3)
      } else{
        B3 = as.vector(B3)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
      }
      
      
      if(SG$gridsizest[lcv1,e] > 1.5){
        B2 <- -dCS[[((e-1)*Q+SG$uo[lcv1,e])]]%*%B2
        B2 <-  backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],backsolve(CCS[[((e-1)*Q+SG$uo[lcv1,e])]],B2, transpose = TRUE))
        B2 = t(B2)
      }else{
        B2 =-as.vector(dCS[[((e-1)*Q+SG$uo[lcv1,e])]])*as.vector(B2)/(as.vector(CCS[[((e-1)*Q+SG$uo[lcv1,e])]])^2)
      }
      
      if(e>1.5){
        for (e2 in  (e-1):1) {
          if(SG$gridsizest[lcv1,e2] > 1.5){
            B2 <- matrix(as.vector(B2),SG$gridsizest[lcv1,e2],SG$gridsizet[lcv1]/SG$gridsizest[lcv1,e2])
            B2 <- t(B2)
          }
          else{
            B2= as.vector(B2)
          }
        }
      }
      
      dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],e] = dpw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],e] +
        SG$w[lcv1] * B2
    }
    
  }
  sigma_hat = t(y) %*% pw / length(y)
  
  
  
  dsigma_hat = t(y) %*% dpw / length(y)
  
  
  
  
  lDet = 0
  
  dlDet = rep(0, SG$d)
  
  for (lcv1 in 1:SG$uoCOUNT) {
    for (lcv2 in 1:SG$d) {
      levelnow = SG$uo[lcv1, lcv2]
      if (levelnow > 1.5) {
        dlDet[lcv2] = dlDet[lcv2] + (dlS[levelnow, lcv2] - dlS[levelnow - 1, lcv2]) * (SG$gridsize[lcv1]) /
          (SG$gridsizes[lcv1, lcv2])
      }
    }
  }
  ddL = dsigma_hat / sigma_hat[1] + 2 / length(y) *theta + dlDet / length(y) 
  return(ddL)
  
}



thetaMLE <- function(SG, y,theta0 = rep(0,SG$d)) {
  x2 = optim(
    theta0,
    fn = lik,
    gr = glik,
    y <- y - mean(y),
    SG = SG,
    method = "L-BFGS-B",
    hessian = FALSE,
    lower =rep(-1,8), upper = rep(2,8)
  )
  return(x2$par)
}

