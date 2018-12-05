rm(list = ls())
library(Rcpp)
source("../R/SGGP_fit_fs.R")
source("../R/SGGP_corr_fs.R")
source("../R/SGGP_create_fs.R")
source("../R/SGGP_append_fs.R")
source("../R/SGGP_pred_fs.R")
source("../R/SGGP_fastcalcassist_fs.R")
sourceCpp("../src/specialkronfunctions.cpp")

borehole <- function(x) {
  rw <- x[, 1] * (0.15 - 0.05) + 0.05
  r <-  x[, 2] * (50000 - 100) + 100
  Tu <- x[, 3] * (115600 - 63070) + 63070
  Hu <- x[, 4] * (1110 - 990) + 990
  Tl <- x[, 5] * (116 - 63.1) + 63.1
  Hl <- x[, 6] * (820 - 700) + 700
  L <-  x[, 7] * (1680 - 1120) + 1120
  Kw <- x[, 8] * (12045 - 9855) + 9855
  
  m1 <- 2 * pi * Tu * (Hu - Hl)
  m2 <- log(r / rw)
  m3 <- 1 + 2 * L * Tu / (m2 * rw ^ 2 * Kw) + Tu / Tl
  return(m1 / m2 / m3)
}


d = 8
testf<-function (x) {  return(borehole(x))} 

Npred <- 1000
library("lhs")
Xp = randomLHS(Npred, d)
Yp = testf(Xp)

SGGP = SGGPcreate(d,801) #create the design.  it has so many entries because i am sloppy
Y = testf(SGGP$design) #the design is $design, simple enough, right?
SGGP = SGGPfit(SGGP,cbind(Y,Y*3,2*Y))
Pred2 = SGGPpred(Xp,SGGP)
SGGP=SGGPappend(SGGP,800)
Y = testf(SGGP$design) #the design is $design, simple enough, right?
SGGP = SGGPfit(SGGP,cbind(Y,Y*3,2*Y))
Pred3 = SGGPpred(Xp,SGGP)
cbind(Pred2$mean[,2],Pred3$mean[,2],Yp^1.001)
cbind(Pred2$mean[,1],Pred3$mean[,1],Yp)


sum(abs(Yp-Pred2$mean)^2)  #prediction should be much better
sum(abs(Yp-Pred2$mean)^2/Pred2$var+log(Pred2$var)) #score should be much better
sum((Yp<= Pred2$mean+1.96*sqrt(Pred2$var))&(Yp>= Pred2$mean-1.96*sqrt(Pred2$var)))  #coverage should be closer to 95 %

# SG=SGappend(SG,800) #add 200 points to the design based on thetahat
# Y = testf(SG$design) #the design is $design, simple enough, right?
# SG = thetaMLE(SG,Y)
# GP = SGGPpred(Xp,SG) #build a full emulator
# sum(abs(Yp-GP$mean)^2)  #prediction should be much better
# sum(abs(Yp-GP$mean)^2/GP$var+log(GP$var)) #score should be much better
# sum((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))  #coverage should be closer to 95 %


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
cHa = (A$vectors)%*%diag(abs(A$values)^(-1/2))%*%t(A$vectors)


source("../R/SGGPlik.R")
theta00 = SG$theta
U <- function(re){
  return(lik(theta00-cHa%*%as.vector(re),SG,Y))
}
grad_U <- function(re){
  return(-((as.vector(glik(theta00-cHa%*%as.vector(re),SG,Y))%*%cHa)))
}
q = rep(0,21)

Ustar = U(rep(0,21))
qstar = rep(0,21)

Uo = U(q)
epsilon = 1/sqrt(sum(grad_U(q)^2))
scalev = 0.5
for(i in 1:200){
  p = rnorm(length(q),0,1)*scalev
  K = sum(p^2)/2/scalev^2
  p = p - epsilon * grad_U(q) 
  qp = q + p
  prev = -p + epsilon * grad_U(qp) 
  
  
  Up = U(qp)
  Kp = sum(prev^2)/2/scalev^2
  # print(theta00-cHa%*%as.vector(q));
  #print(theta00-cHa%*%as.vector(q));
  if(runif(1)<0.75){
  if(runif(1) < exp(Uo-Up+K-Kp)){q=qp;Uo=Up;scalev=exp(log(scalev)+1/sqrt(i+4))}else{scalev=exp(log(scalev)-1/sqrt(i+4));scalev = max(scalev,1/sqrt(length(q))/10)}
  }else{
  if(runif(1) < exp(Uo-Up+K-Kp)){epsilon = exp(log(epsilon)+1/sqrt(i+4)); q = qp; Uo=Up}else{epsilon = exp(log(epsilon)-1/sqrt(i+4))}
  }
  print(scalev)
  print(epsilon)
}
Uo = U(q)
Bs = matrix(0,nrow=100,ncol=21)
for(i in 1:10000){
  p = rnorm(length(q),0,1)*scalev
  K = sum(p^2)/2/scalev^2
  p = p - epsilon * grad_U(q) 
  qp = q + p
  
  prev = -p + epsilon * grad_U(qp) 
  
  
  Up = U(qp)
  Kp = sum(prev^2)/2/scalev^2
#  print(Uo)
#  print(Up)
  # print(theta00-cHa%*%as.vector(q));
  #print(theta00-cHa%*%as.vector(q));
  if(runif(1) < exp(Uo-Up+K-Kp)){q = qp; Uo=Up}
  
  # 
  # 
  # Uo = U(q)
  # for (k in 1:10)
  # {
  #   qp = rnorm(length(q),0,1/sqrt(2))
  #   Up = U(qp)
  #   if(runif(1) < exp(Uo-Up-sum(q^2)+sum(qp^2))){q=qp;Uo=Up}
  # }
  
  
  # print(theta00-cHa%*%as.vector(q));
  #print(theta00-cHa%*%as.vector(q)); 
  if((i%%100)==0){Bs[i/100,]=q; print("saving...")}
}
E = theta00-cHa%*%t(Bs)
plot(Bs[1,])
plot(E[1,])

SG$theta
sqrt(apply(Bs[1:floor(i/100),], 2, var))
sqrt(apply(E[,1:floor(i/100)], 1, var))
sqrt(diag(H^(-1)))


mean((E[1,2:100]-mean(E[1,]))*(E[1,1:99]-mean(E[1,])))/var(E[1,])

1/(1-2*0.2/(1-0.2))