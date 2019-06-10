run_one_CGGP_example <- function(testf, d, N0, Nfinal, batchsize, Npred, plotit=T, plotwith="base") {
  Xp <- matrix(runif(Npred*d), Npred, d)
  Yp = testf(Xp)
  
  timestart <- Sys.time()
  
  require("CGGP")
  SG = CGGPcreate(d=d, batchsize=201)
  Y = testf(SG$design)
  SG = CGGPfit(SG,Y)
  
  
  
  for(c in 1:ceiling((Nfinal-N0)/batchsize)){
    cat(c, " ")
    SG=CGGPappend(SG,batchsize, "MAP") #add 200 points to the design based on thetahat
    Y = testf(SG$design)
    if(c < 10){  #eventually we stop estimating theta because it takes awhile and the estimates dont change that much
      SG = CGGPfit(SG,Y) #estimate the parameter (SG structure is important)
    }
  }
  cat("\n")
  Y = testf(SG$design)
  timelastlogthetaMLEstart <- Sys.time()
  SG = CGGPfit(SG, Y) #do one final parameter estimation

  timelastlogthetaMLEend <- Sys.time()
  
  timepredstart <- Sys.time()
  GP = CGGPpred(xp=Xp,CGGP=SG) #build a full emulator
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
  
  if (plotit) {
    if (plotwith=="base") {
      
      di <- sample(1:nrow(SG$design), 100)
      Y0pred <- CGGPpred(xp=SG$design[di,],CGGP=SG) #,Y,logtheta=logthetaest)
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
    } else if (plotwith == "ggplot2") {
      # library(ggplot2)
      tdf <- data.frame(err=GP$mean-Yp, psd=sqrt(GP$var))
      # ggplot(tdf, aes(x=err, y=psd)) + geom_point()
      values <- data.frame(id=factor(c(1, 2)), value=factor(c(1,2)))
      positions <- data.frame(id=rep(values$id, each=3),
                              x=1.1*c(0,errmax*2,-errmax*2, 0,errmax,-errmax),
                              y=1.1*c(0,errmax,errmax,0,errmax,errmax))
      # Currently we need to manually merge the two together
      datapoly <- merge(values, positions, by = c("id"))
      
      # ggplot(datapoly, aes(x = x, y = y)) +
        # geom_polygon(aes(fill = value, group = id))
      # ggplot(tdf, aes(x=err, y=psd)) + geom_polygon(aes(fill = value, group = id, x=x, y=y), datapoly, alpha=.2) + geom_point() +
        # xlab("Predicted - Actual") + ylab("Predicted error") + coord_cartesian(xlim=c(-errmax,errmax), ylim=c(0,errmax))
      ggplot2::ggplot(tdf, ggplot2::aes_string(x='err', y='psd')) + 
        ggplot2::geom_polygon(ggplot2::aes_string(fill = 'value', group = 'id', x='x', y='y'), datapoly, alpha=.2) + 
        ggplot2::geom_point() +
        ggplot2::xlab("Predicted - Actual") + ggplot2::ylab("Predicted error") + 
        ggplot2::coord_cartesian(xlim=c(-errmax,errmax), ylim=c(0,errmax))
      
    }
  }
}