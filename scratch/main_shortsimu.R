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


 d = 8
testf<-function (x) {  return(borehole(x))} 

N <- 5001
Npred <- 1000
library("lhs")

Xp = randomLHS(Npred, d)
Yp = testf(Xp)

SG = SGcreate(8,801) #create the design.  it has so many entries because i am sloppy
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better

SG=SGappend(SG,800) #add 200 points to the design based on thetahat
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better



SG=SGappend(SG,800) #add 200 points to the design based on thetahat
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better


SG=SGappend(SG,800) #add 200 points to the design based on thetahat
Y = testf(SG$design) #the design is $design, simple enough, right?
SG = thetaMLE(SG,Y)
GP = SGGPpred(Xp,SG) #build a full emulator
sum(abs(Yp-GP$mean)^2)  #prediction should be much better

