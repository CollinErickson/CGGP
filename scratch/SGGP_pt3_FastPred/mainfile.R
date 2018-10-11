#rm(list = ls())

#source("SGcreate.R")
#source("GPstuff.R")
#source("SGGPlik.R")
#source("SGGPappendstuff.R")
#source("SGGPpredstuff.R")

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

#d = 10
#testf<-function (x) {  return(wingweight(x))} 

d = 8
testf<-function (x) {  return(borehole(x))} 

# d = 7
# testf<-function (x) {  return(piston(x))} 

 ##d = 6
 testf<-function (x) {  return(otlcircuit(x))} 
  

N <- 2001
Npred <- 1000
#install.packages(c("lhs"))
library("lhs")

Xp = randomLHS(Npred, d)
Yp = testf(Xp)


#install.packages(c("laGP"))
library(laGP)

#x =randomLHS(N, d)
#y= testf(x)
#formals(aGP)[c("X", "Z", "XX")] <- list(x, y, Xp)
#out3 <- aGP(d = list(max = 20))

SG = SGcreate(rep(0, d), rep(1, d),200)
Y = testf(SG$design)
thetaest = thetaMLE(SG,Y)
for(c in 1:round((N-201)/200)){
  print(c)
  SG=SGappend(SG,200,thetaest)
   Y = testf(SG$design)
  thetaest = thetaMLE(SG,Y)
}
Y = testf(SG$design)
thetaest = thetaMLE(SG,Y,tol = 1e-3)

GP =   SGGPpred(Xp,SG,Y,thetaest)


sum(abs(Yp-out3$mean)^2)
sum(abs(Yp-GP$mean)^2)

sum(abs(Yp-out3$mean)^2/out3$var+log(out3$var))
sum(abs(Yp-GP$mean)^2/GP$var+log(GP$var))

sum((Yp<= out3$mean+1.96*sqrt(out3$var))&(Yp>= out3$mean-1.96*sqrt(out3$var)))
sum((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))

# install.packages("tictoc")
