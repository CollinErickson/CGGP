run_one_SGGP_example <- function(testf, d, N0, Nfinal, batchsize, Npred) {
  Xp <- matrix(runif(Npred*d), Npred, d)
  Yp = testf(Xp)
  
  timestart <- Sys.time()
  
  require("SGGP")
  SG = SGGPcreate(d=d, batchsize=201)
  Y = testf(SG$design)
  SG = SGGPfit(SG,Y)
  
  
  
  for(c in 1:ceiling((Nfinal-N0)/batchsize)){
    cat(c, " ")
    SG=SGGPappend(SG,batchsize) #add 200 points to the design based on thetahat
    Y = testf(SG$design)
    if(c < 10){  #eventually we stop estimating theta because it takes awhile and the estimates dont change that much
      SG = SGGPfit(SG,Y) #estimate the parameter (SG structure is important)
    }
  }
  cat("\n")
  Y = testf(SG$design)
  timelastlogthetaMLEstart <- Sys.time()
  SG = SGGPfit(SG, Y) #do one final parameter estimation

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
  
}