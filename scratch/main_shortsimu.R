rm(list = ls())
library(Rcpp)
source("../R/SGcreate.R")
source("../R/CorrFunctions.R")
source("../R/SGGPlik.R")
source("../R/SGGPappendstuff.R")
source("../R/calculate_pw.R")
source("../R/SGGPpredstuff.R")
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


d = 8
testf<-function (x) {  return(borehole(x))} 

#d = 7
#testf<-function (x) {  return(piston(x))} 

Npred <- 1000
library("lhs")
Xp = randomLHS(Npred, d)
Yp = testf(Xp)

SG = SGcreate(d,801) #create the design.  it has so many entries because i am sloppy
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better
sum(abs(Yp-GP$mean)^2/GP$var+log(GP$var)) #score should be much better
sum((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))  #coverage should be closer to 95 %

SG=SGappend(SG,800) #add 200 points to the design based on thetahat
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better
sum(abs(Yp-GP$mean)^2/GP$var+log(GP$var)) #score should be much better
sum((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))  #coverage should be closer to 95 %

SG=SGappend(SG,800) #add 200 points to the design based on thetahat
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better
sum(abs(Yp-GP$mean)^2/GP$var+log(GP$var)) #score should be much better
sum((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))  #coverage should be closer to 95 %


SG=SGappend(SG,800) #add 200 points to the design based on thetahat
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better
sum(abs(Yp-GP$mean)^2/GP$var+log(GP$var)) #score should be much better
sum((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))  #coverage should be closer to 95 %

SG=SGappend(SG,800) #add 200 points to the design based on thetahat
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better
sum(abs(Yp-GP$mean)^2/GP$var+log(GP$var)) #score should be much better
sum((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))  #coverage should be closer to 95 %


# 
# 
library(laGP)
 
x =randomLHS(dim(SG$design)[1], d)
y= testf(x)
formals(aGP)[c("X", "Z", "XX")] <- list(x, y, Xp)
out3 <- aGP(d = list(max = 20))
sum(abs(Yp-out3$mean)^2)  #prediction should be much better
sum(abs(Yp-out3$mean)^2/out3$var+log(out3$var)) #score should be much better
sum((Yp<= out3$mean+1.96*sqrt(out3$var))&(Yp>= out3$mean-1.96*sqrt(out3$var)))  #coverage should be closer to 95 %
 
