rm(list = ls())

source("../R/SGcreate.R")
source("../R/CorrMatern32.R")
source("../R/SGGPlik.R")
source("../R/SGGPappendstuff.R")
source("../R/SGGPpredstuff.R")
source("../R/calculate_pw.R")

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

piston <- function(xx)
{
  M  <- xx[,1]*30 + 30
  S  <- xx[,2]*0.015 + 0.005
  V0 <- xx[,3]*0.008 + 0.002
  k  <- xx[,4]*4000 + 1000
  P0 <- xx[,5]*20000+90000
  Ta <- xx[,6]*6 + 290
  T0 <- xx[,7]*20 + 340
  
  Aterm1 <- P0 * S
  Aterm2 <- 19.62 * M
  Aterm3 <- -k*V0 / S
  A <- Aterm1 + Aterm2 + Aterm3
  
  Vfact1 <- S / (2*k)
  Vfact2 <- sqrt(A^2 + 4*k*(P0*V0/T0)*Ta)
  V <- Vfact1 * (Vfact2 - A)
  
  fact1 <- M
  fact2 <- k + (S^2)*(P0*V0/T0)*(Ta/(V^2))
  
  C <- 2 * pi * sqrt(fact1/fact2)
  return(C)
}


wingweight <- function(xx)
{
  Sw      <- xx[,1]*50+150
  Wfw     <- xx[,2]*80+220
  A       <- xx[,3]*4 + 6
  LamCaps <- (xx[,4]*20-10)*pi/180
  q       <- xx[,5]*(45-16)+16
  lam     <- xx[,6]*0.5+0.5
  tc      <- xx[,7]*0.1+0.08
  Nz      <- xx[,8]*3.5 + 2.5
  Wdg     <- xx[,9]*800+1700
  Wp      <- xx[,10]*(0.08-0.025)+0.025
  
  fact1 <- 0.036 * Sw^0.758 * Wfw^0.0035
  fact2 <- (A / ((cos(LamCaps))^2))^0.6
  fact3 <- q^0.006 * lam^0.04
  fact4 <- (100*tc / cos(LamCaps))^(-0.3)
  fact5 <- (Nz*Wdg)^0.49
  
  term1 <- Sw * Wp
  
  y <- fact1*fact2*fact3*fact4*fact5 + term1
  return(y)
}

otlcircuit <- function(xx)
{
  Rb1  <- xx[,1]*100+50
  Rb2  <- xx[,2]*45 + 25
  Rf   <- xx[,3]*2.5 + 0.5
  Rc1  <- xx[,4]*1.3 + 1.2
  Rc2  <- xx[,5]*.95 + .25
  beta <- xx[,6]*250+50
  
  Vb1 <- 12*Rb2 / (Rb1+Rb2)
  term1a <- (Vb1+0.74) * beta * (Rc2+9)
  term1b <- beta*(Rc2+9) + Rf
  term1 <- term1a / term1b
  
  term2a <- 11.35 * Rf
  term2b <- beta*(Rc2+9) + Rf
  term2 <- term2a / term2b
  
  term3a <- 0.74 * Rf * beta * (Rc2+9)
  term3b <- (beta*(Rc2+9)+Rf) * Rc1
  term3 <- term3a / term3b
  
  Vm <- term1 + term2 + term3
  return(Vm)
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

SG = SGcreate(d,35001) #create the design.  it has so many entries because i am sloppy
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
SG = logthetaMLE(SG,Y,tol = 1e-3) #do one final parameter estimation,  this should be speed up, but I was lazy

library(Rcpp)
#library(RcppArmadillo)
library(tictoc)
source("../R/calculate_pw.R")
A = tic()
P1= lik(rep(0,8),SG=SG,y=Y) 
P2= lik(rep(0,8)+c(10^(-5),0,0,0,0,0,0,0),SG=SG,y=Y)  
sum((P2-P1)*100000)
tic()
dP1= glik(rep(0,8),SG=SG,y=Y)  
toc()
dP1[1]
sourceCpp("../src/specialkronfunctions.cpp")
source("calculate_pw_test.R")
A = tic()
P1= lik(rep(0,8),SG=SG,y=Y)  
P2= lik(rep(0,8)+c(10^(-5),0,0,0,0,0,0,0),SG=SG,y=Y)  
sum((P2-P1)*100000)
tic()
dP1= glik(rep(0,8),SG=SG,y=Y) 
toc()
dP1[1]
