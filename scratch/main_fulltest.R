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
  
  Yn = m1 / m2 / m3
  return(abs(cbind(Yn,Yn^0.75,Yn^0.5,Yn^1.25)))
 # return(Yn)
}


d = 8
testf<-function (x) {  return(borehole(x))} 

Npred <- 1000
library("lhs")
Xp = randomLHS(Npred, d)
Yp = testf(Xp)

SGGP = SGGPcreate(d,1001) #create the design.  it has so many entries because i am sloppy
Y = testf(SGGP$design) #the design is $design, simple enough, right?
source("../R/SGGP_fit_fs.R")
SGGP = SGGPfit(SGGP,Y)
sqrt(colMeans(t(SGGP$thetaPostSamples^2))-colMeans(t(SGGP$thetaPostSamples))^2)

Pred = SGGPpred(Xp,SGGP)

SGGP=SGGPappend(SGGP,800)
Y = testf(SGGP$design) #the design is $design, simple enough, right?
SGGP = SGGPfit(SGGP,Y)

SGGP=SGGPappend(SGGP,800)
Y = testf(SGGP$design) #the design is $design, simple enough, right?
SGGP = SGGPfit(SGGP,Y)
colMeans(t(SGGP$thetaPostSamples^2))-colMeans(t(SGGP$thetaPostSamples))^2

mean(abs(Yp-Pred4$mean)^2)  #prediction should be much better
mean(abs(Yp-Pred4$mean)^2/Pred4$var+log(Pred4$var)) #score should be much better
mean((Yp<= Pred4$mean+1.96*sqrt(Pred4$var))&(Yp>= Pred4$mean-1.96*sqrt(Pred4$var)))  #coverage should be closer to 95 %
