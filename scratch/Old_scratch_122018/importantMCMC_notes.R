H = rep(0,21)
G = lik(SG$theta,SG,Y)
for(i in 1:21){
  rsad = rep(0,21)
  rsad[i] =10^(-4)
  H[i] = (lik(SG$theta+rsad,SG,Y)-G)*10^(4)
}

glik(SG$theta,SG,Y)


H = matrix(0,nrow=21,ncol=21)
library(tictoc)
tic()
G = glik(SG$theta,SG,Y)
for(c in 1:21){
  rsad = rep(0,21)
  rsad[c] =10^(-4)
  H[c,] = (glik(SG$theta+rsad,SG,Y)-G)*10^(4)
}
H = H/2+t(H)/2
A = eigen(H,TRUE)
Happrox = solve((A$vectors)%*%diag(abs(A$values))%*%t(A$vectors))
cHa = t(chol(Happrox))


source("../R/SGGPlik.R")
theta00 = SG$theta
U <- function(re){
  return(lik(theta00-cHa%*%as.vector(re),SG,Y)*length(Y)/2)
}
grad_U <- function(re){
  return(-((as.vector(glik(theta00-cHa%*%as.vector(re),SG,Y))%*%cHa))*length(Y)/2)
}
q = rep(0,21)

epsilon = 10^(-2.5)
q =rep(0,21)
for(c in 1:100){
  current_q = q
  p = rnorm(length(q),0,1)
  current_p = p
  
  p = p - epsilon * grad_U(q) / 2
  for (i in 1:5)
  {
    q = q + epsilon * p
    p = p - epsilon * grad_U(q)
  }
  
  p = p + epsilon * grad_U(q) / 2
  p = -p
  current_U = U(current_q)
  print(U(current_q))
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  print(U(q))
  proposed_K = sum(p^2) / 2
  if(runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){ print(theta00-cHa%*%as.vector(q)); print("new")}else{q=current_q; print(theta00-cHa%*%as.vector(q)); print("old")}
}
