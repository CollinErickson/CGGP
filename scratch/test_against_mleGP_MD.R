rm(list = ls())

boreholeMV <- function(x) {
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
  
  Tot_v = abs(m1 / m2 / m3)
  G1=1/100*(Tot_v)^matrix(seq(0.5,1.5,0.02),nrow=dim(x)[1],ncol=51,byrow=TRUE)
  G = cos(0.05*G1)*G1*10
  return(G)
  
}

d = 8
testf<-function (x) {  return(boreholeMV(x))} 


Npred <- 1000
library("lhs")

Xp = randomLHS(Npred, d)
Yp = testf(Xp)


SGGP = SGGPcreate(d,100) #create the design.  it has so many entries because i am sloppy
Y = testf(SGGP$design) #the design is $design, simple enough, right?
SGGP = SGGPfit(SGGP,Y)
SGGPTS=SGGPappend(SGGP,200,selectionmethod="TS")
YTS = testf(SGGPTS$design) #the design is $design, simple enough, right?
SGGPTS = SGGPfit(SGGPTS,YTS)
SGGPTS=SGGPappend(SGGPTS,300,selectionmethod="TS")
YTS = testf(SGGPTS$design) #the design is $design, simple enough, right?
SGGPTS = SGGPfit(SGGPTS,YTS)
PredTS = SGGPpred(Xp,SGGPTS)
mean(abs(Yp-PredTS$mean)^2)  #prediction should be much better

library(mlegp)
Xmlegp = randomLHS(200, d)
Ymlegp = testf(Xmlegp)
fitMulti_PCA = mlegp(Xmlegp,t(Ymlegp),PC.num=5)
predmlegp_PCA = fitMulti_PCA$UD%*%t(cbind(predict(fitMulti_PCA[[1]],Xp),predict(fitMulti_PCA[[2]],Xp),predict(fitMulti_PCA[[3]],Xp),predict(fitMulti_PCA[[4]],Xp),predict(fitMulti_PCA[[5]],Xp)))


outputindex = 51
SGGPTS2 = SGGPfit(SGGPTS,YTS[,outputindex])
PredTS2 = SGGPpred(Xp,SGGPTS2)
fitMulti_NOPCA = mlegp(Xmlegp,Ymlegp[,outputindex])
predmlegp_NOPCA =predict(fitMulti_NOPCA,Xp)
mean(abs(Yp[,outputindex]-PredTS$mean[,outputindex])^2) 
mean(abs(Yp[,outputindex]-PredTS2$mean)^2) 
mean(abs(Yp[,outputindex]-t(predmlegp_PCA)[,outputindex])^2) 
mean(abs(Yp[,outputindex]-(predmlegp_NOPCA))^2)  #prediction should be much better

