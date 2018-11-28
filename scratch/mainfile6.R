rm(list = ls())
library(Rcpp)
source("../R/SGcreate.R")
source("../R/CorrMatern32.R")
source("../R/SGGPlik.R")
source("../R/SGGPappendstuff.R")
source("../R/SGGPpredstuff.R")
source("../R/calculate_pw.R")
source("../R/calculate_pw_C.R")
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

# d = 10
# testf<-function (x) {  return(wingweight(x))} 

 d = 8
testf<-function (x) {  return(borehole(x))} 

# d = 7
# testf<-function (x) {  return(piston(x))} 

# d = 6
# testf<-function (x) {  return(otlcircuit(x))} 
  

N <- 5001
Npred <- 1000
#install.packages(c("lhs"))
library("lhs")

Xp = randomLHS(Npred, d)
Yp = testf(Xp)


x =randomLHS(N, d)
y= testf(x)

SG = SGcreate(d,2001) #create the design.  it has so many entries because i am sloppy
Y = testf(SG$design) #the design is $design, simple enough, right?
#SG = logthetaMLE(SG,Y)

#for(c in 1:round((N-201)/200)){
#  print(c)
#  SG=SGappend(SG,200,theta=exp(SG$logtheta)) #add 200 points to the design based on thetahat
#  Y = testf(SG$design)
#  if( c< 10 ){  #eventually we stop estimating theta because it takes awhile and the estimates dont change that much
#    SG = logthetaMLE(SG,Y) #estimate the parameter (SG structure is important)
#  }
#}
#Y = testf(SG$design)
source("../R/SGGPlik.R")
SG = logthetaMLE(SG,Y,tol = 1e-3) #do one final parameter estimation,  this should be speed up, but I was lazy
SG$logtheta
SG = logthetaMLE(SG,Y,tol = 1e-3,logtheta0=rep(1,8)) #do one final parameter estimation,  this should be speed up, but I was lazy
SG$logtheta
SG = logthetaMLE(SG,Y,tol = 1e-3,logtheta0=c(-1,2,2,2,1,2,-1,-1)) #do one final parameter estimation,  this should be speed up, but I was lazy
SG$logtheta


source("../R/calculate_pw.R")
P1= posterior(rep(0,8),SG=SG,y=Y) 
P1= posterior(rep(0,8),SG=SG,y=Y) 
P2= posterior(rep(0,8)+c(10^(-5),0,0,0,0,0,0,0),SG=SG,y=Y)  
sum((P2-P1)*100000)
dP1= gposterior(SG$logtheta,SG=SG,y=Y)  
dP1[1]
dP2 = matrix(0,nrow=SG$d,ncol=SG$d)
for(c in 1:d){
  dt = rep(0,SG$d)
  dt[c]=10^(-5)
dP2[c,]= (gposterior(SG$logtheta+dt,SG=SG,y=Y)-dP1 )/10^(-5)
}

Lim = qchisq(0.95,SG$d+2)
Pstar = posterior(SG$logtheta,SG=SG,y=Y) 

x0 = SG$logtheta
x0
(posterior(rep(0,8),SG=SG,y=Y)-Pstar+Lim)
for(t in 1:10){
  print(t)
dphi = -c(0,1,0,rep(0,5))-1/t*gposterior(x0,SG=SG,y=Y)/(posterior(x0,SG=SG,y=Y)-Pstar+Lim)-1/t*rep(1,8)/(4-x0)
#linesearch
eps = 0.01/sum(abs(dphi))
xp = x0 + eps*dphi
Pp = posterior(xp,SG=SG,y=Y)-Pstar+Lim
while(Pp > 0){
  eps = eps*2
  xp = x0 + eps*dphi
  Pp = posterior(xp,SG=SG,y=Y)-Pstar+Lim
}
while(Pp < 0){
  eps = eps/2
  xp = x0 + eps*dphi
  Pp = posterior(xp,SG=SG,y=Y)-Pstar+Lim
}
x0 = x0 + eps*dphi
}
x0



source("../R/CorrFunctions.R")
S1 = CorrMatCauchy(c(1,2,3),c(1.5,2.5),LS = 1)
S2 = CorrMatCauchy(c(1,2,3),c(1.5,2.5),LS = 1+10^(-3))
(S2-S1)*10^(3)
Sstuff = CorrMatCauchy(c(1,2,3),c(1.5,2.5),LS = 1,return_dCdLS = TRUE, return_dCdFD = TRUE, return_dCdHE = TRUE)
Sstuff$dCdLS
S1 = CorrMatCauchy(c(1,2,3),c(1.5,2.5),HE = 1)
S2 = CorrMatCauchy(c(1,2,3),c(1.5,2.5),HE = 1+10^(-3))
(S2-S1)*10^(3)
Sstuff = CorrMatCauchy(c(1,2,3),c(1.5,2.5),HE = 1,return_dCdLS = TRUE, return_dCdFD = TRUE, return_dCdHE = TRUE)
Sstuff$dCdHE
S1 = CorrMatCauchy(c(1,2,3),c(1.5,2.5),HE = 0.5,FD = 1)
S2 = CorrMatCauchy(c(1,2,3),c(1.5,2.5),HE = 0.5,FD = 1+10^(-3))
(S2-S1)*10^(3)
Sstuff = CorrMatCauchy(c(1,2,3),c(1.5,2.5),HE = 0.5,FD = 1,return_dCdLS = TRUE, return_dCdFD = TRUE, return_dCdHE = TRUE)
Sstuff$dCdFD