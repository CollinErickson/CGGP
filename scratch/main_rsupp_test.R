rm(list = ls())
library(SGGP)

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
 # return(abs(cbind(Yn,Yn^0.75,Yn^0.5,Yn^1.1)))
  return(Yn)
}


d = 8
testf<-function (x) {  return(borehole(x))} 

Npred <- 1000
library("lhs")
Xp = randomLHS(Npred, d)
Yp = testf(Xp)

Xs1 = randomLHS(100, d)
Ys1 = testf(Xs1)
Xs2 = randomLHS(200, d)
Ys2 = testf(Xs2)
Xs3 = randomLHS(300, d)
Ys3 = testf(Xs3)

SGGP = SGGPcreate(d,200) #create the design.  it has so many entries because i am sloppy
Y = testf(SGGP$design) 
for(MethodSupp in c( "Ignore", "MarginalValidation","FullValidation")){
SGGP1 = SGGPfit(SGGP,Y,Xs=Xs1,Ys=Ys1,SuppData=MethodSupp)
PredGreedy = SGGPpred(SGGP1,Xp)
SSE1 = mean(abs(Yp-PredGreedy$mean)^2) 
Score1=mean(abs(Yp-PredGreedy$mean)^2/PredGreedy$var+log(PredGreedy$var)) #score should be much better
SGGP1 = SGGPappend(SGGP1,200) #create the design.  it has so many entries because i am sloppy
Y1 = testf(SGGP1$design) #the design is $design, simple enough, right?
SGGP1 = SGGPfit(SGGP1,Y1,Xs=Xs1,Ys=Ys1,SuppData=MethodSupp)
PredGreedy = SGGPpred(SGGP1,Xp)
SSE2 = mean(abs(Yp-PredGreedy$mean)^2) 
Score2=mean(abs(Yp-PredGreedy$mean)^2/PredGreedy$var+log(PredGreedy$var)) #score should be much better
SGGP1 = SGGPappend(SGGP1,200) #create the design.  it has so many entries because i am sloppy
Y1 = testf(SGGP1$design) #the design is $design, simple enough, right?
SGGP1 = SGGPfit(SGGP1,Y1,Xs=Xs1,Ys=Ys1,SuppData=MethodSupp)
PredGreedy = SGGPpred(SGGP1,Xp)
SSE3 =mean(abs(Yp-PredGreedy$mean)^2) 
Score3=mean(abs(Yp-PredGreedy$mean)^2/PredGreedy$var+log(PredGreedy$var)) #score should be much betterv
SGGP1 = SGGPfit(SGGP1,Y1,Xs=Xs2,Ys=Ys2,SuppData=MethodSupp)
PredGreedy = SGGPpred(SGGP1,Xp)
SSE4 =mean(abs(Yp-PredGreedy$mean)^2) 
Score4=mean(abs(Yp-PredGreedy$mean)^2/PredGreedy$var+log(PredGreedy$var)) #score should be much betterr
SGGP1 = SGGPappend(SGGP1,200) #create the design.  it has so many entries because i am sloppy
Y1 = testf(SGGP1$design) #the design is $design, simple enough, right?
SGGP1 = SGGPfit(SGGP1,Y1,Xs=Xs2,Ys=Ys2,SuppData=MethodSupp)
PredGreedy = SGGPpred(SGGP1,Xp)
SSE5 =mean(abs(Yp-PredGreedy$mean)^2) 
Score5= mean(abs(Yp-PredGreedy$mean)^2/PredGreedy$var+log(PredGreedy$var)) #score should be much betterr
SGGP1 = SGGPfit(SGGP1,Y1,Xs=Xs3,Ys=Ys3,SuppData=MethodSupp)
PredGreedy = SGGPpred(SGGP1,Xp)
SSE6 =mean(abs(Yp-PredGreedy$mean)^2) 
Score6= mean(abs(Yp-PredGreedy$mean)^2/PredGreedy$var+log(PredGreedy$var)) #score should be much betterr
print(MethodSupp)
print(cbind(c(300,500,700,800,1000,1100),c(SSE1,SSE2,SSE3,SSE4,SSE5,SSE6),c(Score1,Score2,Score3,Score4,Score5,Score6)))
}