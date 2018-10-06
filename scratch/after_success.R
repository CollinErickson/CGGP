timestamp()
print("RUNNING after_success.R ...")

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

N <- 10001
Npred <- 1000
#install.packages(c("lhs"))
#library("lhs")
if (!require('lhs', quietly = TRUE)) {
  install.packages('lhs')
  require('lhs')
}

Xp = randomLHS(Npred, d)
Yp = testf(Xp)

goodlogthetaest_old <- c(-0.01932437,  0.82517131,  0.88499983,  0.73263796,  0.86971878,  0.70425694,  0.65443469,  0.80910334)
goodlogthetaest <- log(exp(goodlogthetaest_old)/sqrt(3))
use_goodtheta <- FALSE


SG = SGcreate(rep(0, d), rep(1, d),201) #create the design.  it has so many entries because i am sloppy
Y = testf(SG$design) #the design is $design, simple enough, right?
logthetaest = logthetaMLE(SG,Y)
if (use_goodtheta) logthetaest <- goodlogthetaest
thetaest <- exp(logthetaest)
print(logthetaest)

for(c in 1:round((N-201)/200)){
  print(c)
  SG=SGappend(SG,200,theta=thetaest) #add 200 points to the design based on thetahat
  Y = testf(SG$design)
  if( c< 10 && !use_goodtheta){  #eventually we stop estimating theta because it takes awhile and the estimates dont change that much
    logthetaest = logthetaMLE(SG,Y) #estimate the parameter (SG structure is important)
    thetaest <- exp(logthetaest)
    print(logthetaest)
  }
}
Y = testf(SG$design)
if (!use_goodtheta) {
  logthetaest = logthetaMLE(SG,Y,tol = 1e-3) #do one final parameter estimation,  this should be speed up, but I was lazy
  print(logthetaest)
}
GP = SGGPpred(Xp,SG,Y,logtheta=pmin(logthetaest,2)) #build a full emulator

SumSE <- sum((Yp-GP$mean)^2)  #prediction should be much better
score <- sum((Yp-GP$mean)^2/GP$var+log(GP$var)) #score should be much better
coverage <- sum((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))  #coverage should be closer to 95 %

print(paste("SumSE is     ", SumSE))
print(paste("Score is   ", score))
print(paste("coverage is", coverage))


print("... FINISHED after_success.R")
timestamp()
