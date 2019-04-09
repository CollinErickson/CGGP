
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

Npred <- 2000
library("lhs")
Xp = randomLHS(Npred, d)
Yp = testf(Xp)

Xs = randomLHS(40, d)
Ys = testf(Xs)
CGGP = CGGPcreate(d,300) 
Y = testf(CGGP$design) 
CGGP = CGGPfit(CGGP,Y,Xs=Xs,Ys=Ys)
CGGP = CGGPappend(CGGP,600) 
library(tictoc)
tic('here')
Pred = CGGPpred(CGGP,Xp)
toc()

Y = testf(CGGP$design) 
CGGP = CGGPfit(CGGP,Y,Xs=Xs,Ys=Ys)

source('./CGGP_impute_fs.R')
Y = testf(CGGP$design) 
CGGP = CGGPfit(CGGP,Y,Xs=Xs,Ys=Ys)
Y = Y - mean(Y)
Yb= Y
Yb[261] = NA
Yb[14] = NA
Yb[56] = NA
Yb[660] = NA
Yb[760] = NA
CGGPimpute(CGGP,Yb,Y)

