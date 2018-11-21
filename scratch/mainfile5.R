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

d = 8
testf<-function (x) {  return(borehole(x))} 

N <- 5001
Npred <- 1000
library("lhs")

Xp = randomLHS(Npred, d)
Yp = testf(Xp)


x =randomLHS(N, d)
y= testf(x)

SG = SGcreate(d,20001) #create the design.  it has so many entries because i am sloppy
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
#SG = logthetaMLE(SG,Y,tol = 1e-3) #do one final parameter estimation,  this should be speed up, but I was lazy
library(tictoc)
library(Rcpp)

source("../R/calculate_pw.R")
source("../R/calculate_pw_C.R")
A = tic()
R1 <- calculate_pw_and_dpw(SG, Y, rep(0,8))
R2 <- calculate_pw_and_dpw(SG, Y, c(0,0,10^(-3),rep(0,5))+rep(0,8))
toc()
R1$pw[200:210]
(R1$pw[200:210]-R2$pw[200:210])/10^(-3)
R1$dpw[200:210,3]

sourceCpp("../src/specialkronfunctions.cpp")

source("../R/calculate_pw_C.R")
A = tic()
R1 <- calculate_pw_and_dpw_C(SG, Y, rep(0,8))
R2 <- calculate_pw_and_dpw_C(SG, Y, c(0,0,10^(-3),rep(0,5))+rep(0,8))
toc()
R1$pw[200:210]
(R1$pw[200:210]-R2$pw[200:210])/10^(-3)
R1$dpw[200:210,3]

