timestamp()
cat("RUNNING after_success.R ...\n")
timestart <- Sys.time()

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
testf <- borehole

N <- 10001
Npred <- 1000
#install.packages(c("lhs"))
#library("lhs")
#if (!require('lhs', quietly = TRUE)) {
#  install.packages('lhs')
#  require('lhs')
#}

#Xp = randomLHS(Npred, d)
Xp <- matrix(runif(Npred*d), Npred, d)
Yp = testf(Xp)

# goodlogthetaest_old <- c(-0.01932437,  0.82517131,  0.88499983,  0.73263796,  0.86971878,  0.70425694,  0.65443469,  0.80910334)
# goodlogthetaest <- log(exp(goodlogthetaest_old)/sqrt(3))
# use_goodtheta <- FALSE

require("SGGP")
SG = SGGPcreate(d=d, batchsize=201)
Y = testf(SG$design) #the design is $design, simple enough, right?
# logthetaest = logthetaMLE(SG,Y)
SG = SGGPfit(SG, Y)
# logthetaest <- SG$logtheta
# if (use_goodtheta) logthetaest <- goodlogthetaest
# thetaest <- exp(logthetaest)
# cat(logthetaest, "\n")
cat("Now doing Bayesian\n")

for(c in 1:round((N-201)/200)){
  cat(c, " ")
  SG=SGGPappend(SG,200) #add 200 points to the design based on thetahat
  Y = testf(SG$design)
  if( c< 10){  #eventually we stop estimating theta because it takes awhile and the estimates dont change that much
    SG = SGGPfit(SG,Y) #estimate the parameter (SG structure is important)
    # logthetaest <- SG$logtheta
    # thetaest <- exp(logthetaest)
    # cat(logthetaest,"\n", sep="\t")
  }
}
cat("\n")
Y = testf(SG$design)
timelastlogthetaMLEstart <- Sys.time()
# if (!use_goodtheta) {
  SG = SGGPfit(SG, Y) #do one final parameter estimation,  this should be speed up, but I was lazy
  # logthetaest <- SG$logtheta
  # cat(logthetaest, "\n")
# }
timelastlogthetaMLEend <- Sys.time()

timepredstart <- Sys.time()
GP = SGGPpred(Xp,SG) #build a full emulator
timepredend <- Sys.time()

RMSE <- sqrt(mean(((Yp-GP$mean)^2)))  #prediction should be much better
meanscore <- mean((Yp-GP$mean)^2/GP$var+log(GP$var)) #score should be much better
meancoverage <- mean((Yp<= GP$mean+1.96*sqrt(GP$var))&(Yp>= GP$mean-1.96*sqrt(GP$var)))  #coverage should be closer to 95 %

cat("RMSE is         ", RMSE, "\n")
cat("Mean score is   ", meanscore, "\n")
cat("coverage is     ", meancoverage, "\n")

# Don't count plotting in run time
timeend <- Sys.time()
cat("Total run   time is:", capture.output(timeend - timestart), '\n')
cat("Prediction  time is:", capture.output(timepredend - timepredstart), '\n')
cat("logthetaMLE time is:", capture.output(timelastlogthetaMLEend - timelastlogthetaMLEstart), '\n')

if (T) { # Can Travis just skip this?
  di <- sample(1:nrow(SG$design), 100)
  Y0pred <- SGGPpred(SG$design[di,],SG) #,Y,logtheta=logthetaest)
  plot(Yp, GP$mean, ylim=c(min(GP$mean, Y0pred$m),max(GP$mean, Y0pred$m))); points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
  # Now plot with bars
  #plot(Yp, GP$mean , ylim=c(min(GP$mean, Y0pred$m),max(GP$mean, Y0pred$m)),pch=19)#; points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
  plot(Yp, Yp-GP$mean , ylim=max(sqrt(GP$var))*c(-2,2))#c(min(-2GP$mean, Y0pred$m),max(GP$mean, Y0pred$m)),pch=19,col='white')#; points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
  points(Yp, 0*GP$mean + 2*sqrt(GP$var), col=4, pch=19)#; points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
  points(Yp, 0*GP$mean - 2*sqrt(GP$var), col=5)#; points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
  errmax <- max(sqrt(GP$var), abs(GP$mean - Yp))
  plot(GP$mean-Yp, sqrt(GP$var), xlim=errmax*c(-1,1), ylim=c(0,errmax))#;abline(a=0,b=1,col=2)
  polygon(1.1*errmax*c(0,-2,2),1.1*errmax*c(0,1,1), col=3, density=10, angle=135)
  polygon(1.1*errmax*c(0,-1,1),1.1*errmax*c(0,1,1), col=2, density=30)
  points(GP$mean-Yp, sqrt(GP$var), xlim=errmax*c(-1,1), ylim=c(0,errmax))
}

# cat(capture.output(Sys.time() - timestart), '\n')
cat("... FINISHED after_success.R\n")
timestamp()


cat("Running second example: wingweight\n")
print(getwd())
source("./scratch/after_success_run_one.R")
run_one_SGGP_example(TestFunctions::wingweight, 10, 1000, 6000, 1000, 1000)
