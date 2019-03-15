SGGPpred2 <- function(xp,SG, y, ..., logtheta, theta) {
  if (missing(theta)) {theta <- exp(logtheta)}
  # Center outputs
  my = mean(y)
  y = y-my
  
  Q  = max(SG$uo[1:SG$uoCOUNT,]) # Max level of evaluated blocks
  # Now store Choleskys instead of inverses
  #CiS = list(matrix(1,1,1),Q*SG$d) # Store correlation matrices
  CS = list(matrix(1,1,1),Q*SG$d)
  CCS = list(matrix(1,1,1),Q*SG$d)
  # Loop over dimensions and possible levels
  for (lcv2 in 1:SG$d) {
    for (lcv1 in 1:max(SG$uo[1:SG$uoCOUNT,lcv2])) {
      Xbrn = SG$xb[1:SG$sizest[lcv1]] # Get x's
      Xbrn = Xbrn[order(Xbrn)] # Sort them
      S = SG$CorrMat(Xbrn, Xbrn , theta=theta[lcv2]) # Calculate corr mat
      diag(S) = diag(S) + SG$nugget
      #CiS[[(lcv2-1)*Q+lcv1]] = solve(S) # Store inversion
      CS[[(lcv2-1)*Q+lcv1]] = S
      CCS[[(lcv2-1)*Q+lcv1]] = chol(S)
      #print(eigen(S)$val)
      #print(det(S))
    }
  }
  
  
  pw = rep(0, length(y)) # ????????????
  pred_v = rep(0, dim(xp)[1]) # ????????????
  
  Corrstore = array(1,c(dim(xp)[1],length(SG$xb),d))
  for (e in 1:SG$d) { # Loop over dimensions
    Corrstore[,,e] = SG$CorrMat(xp[,e], SG$xb, theta=theta[e])
  }
  
  
  Cstor = matrix(1,nrow = dim(xp)[1], ncol = dim(SG$dit)[2])
  
  # Loop over blocks
  for (lcv1 in 1:SG$uoCOUNT) {
    
    
    B = y[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]]
    Av = unlist(CCS[((1:d-1)*Q+SG$uo[lcv1,1:d])])
    B2=B
    B = kronDBS(Av, B2, SG$gridsizest[lcv1,], length(Av), length(B2),  d)
    
    pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] = pw[SG$dit[lcv1, 1:SG$gridsizet[lcv1]]] +
      SG$w[lcv1] * B
    
    Cstor  = matrix(1,nrow = dim(xp)[1], ncol =length(as.vector(B)))
    for(e in 1:SG$d){
      if(SG$gridsizest[lcv1,e] > 1.5){
        Cstor=  Cstor*as.matrix(Corrstore[,SG$designindex[SG$dit[lcv1, 1:SG$gridsizet[lcv1]],e],e])
      }
    }
   # print(dim(as.matrix(Cstor)))
    #print(dim(as.matrix(as.vector(B))))
   pred_adj = as.matrix(Cstor)%*%(as.vector(B))
    for(e in 1:SG$d){
      if(SG$gridsizest[lcv1,e] < 1.5){
        pred_adj=  pred_adj*as.vector(Corrstore[,1,e])
      }
    }
    
    pred_v = pred_v+SG$w[lcv1]*pred_adj
  }
  
  
  sigma_hat = t(y) %*% pw / length(y)
  
  Rv = matrix(1,dim(xp)[1],length(SG$xb))
  Cp = matrix(1,dim(xp)[1],SG$ss)
  
  MSE_v = array(0, c(SG$d, 9,dim(xp)[1]))
  for (lcv1 in 1:SG$d) {
    MSE_v[lcv1, 1,] = SG$diag_corrMat(xp[,lcv1], theta=theta[lcv1], nugget=SG$nugget)
  }
  for (lcv1 in 1:SG$d) {
    for (lcv2 in 1:8) {
      MSE_v[lcv1, lcv2+1,] = abs(MSEpred_calc(xp[,lcv1],SG$xb[1:SG$sizest[lcv2]],
                                              theta=theta[lcv1], nugget=SG$nugget,
                                              CorrMat=SG$CorrMat,
                                              diag_corrMat=SG$diag_corrMat))
      MSE_v[lcv1, lcv2+1,] = pmin(MSE_v[lcv1, lcv2+1,], MSE_v[lcv1, lcv2,])
    }
  }
  
  Fmat <- matrix(MSE_v, nrow = d*9)
  ind_mat2 =SG$uo[1:SG$uoCOUNT,]*d + matrix(seq(1,d),nrow = SG$uoCOUNT, ncol = d,byrow=TRUE)
  
  MEt = rep(1,dim(Fmat)[2])
  
  ME_tt = MEEAGA(Fmat,ind_mat2,SG$w,MEt,as.integer(SG$uoCOUNT),as.integer(dim(Fmat)[2]),as.integer(SG$d))
  
  
#  Cp = matrix(1,dim(xp)[1],SG$ss)
#  for (e in 1:SG$d) { # Loop over dimensions
    #Cp = Cp*CorrMat(xp[,e], SG$design[,e], theta=theta[e]) # Multiply correlation from each dimension
 #   V = SG$CorrMat(xp[,e], SG$xb, theta=theta[e])
#    Cp = Cp*V[,SG$designindex[,e]]
 # }
  
  # Return list with mean and var predictions
  GP = list("mean" = (my+pred_adj), "var"=sigma_hat[1]*ME_tt)
  
  return(GP)
}


#' Multivariate SGGP prediction
#' 
#' Predict output for an SGGP with multivariate output.
#'
#' @param yMV Output matrix, each row is for a design point
#' @inheritParams SGGPpred
#'
#' @return Matrix of predictions
#' @export
#'
#' @examples
#' 
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' y1 <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' y2 <- apply(SG$design, 1, function(x){x[1]^1.3+.4*sin(6*x[2])+rnorm(1,0,.01)})
#' y <- cbind(y1, y2)
#' tx <- cbind(SGGPpredMV(SG$design, SG=SG, yMV=y, theta=c(.1,.1,.1))$mean, y)
#' tx # Columns 1 and 3, and 2 and 4 should be near equal
#' cor(tx)
SGGPpredMV <- function(xp,SG, yMV, ..., logtheta, theta) {
  
  p = dim(yMV)[2]
  
  meanMV = matrix(1,dim(xp)[1],p)
  varMV = matrix(1,dim(xp)[1],p)
  for(i in 1:p){
    GP1d = SGGPpred(xp,SG,as.vector(yMV[,i]),logtheta=logtheta, theta=theta)
    
    meanMV[,i] = GP1d$mean
    varMV[,i] = GP1d$var
  }
  
  # Return list with mean and var predictions
  GPMV = list("mean" = meanMV+ pred_v, "var"= varMV)
  
  return(GPMV)
}
