fff1 <- function(nn) {
  
  SG = SGcreate(rep(0, d), rep(1, d),1001) #create the design.  it has so many entries because i am sloppy
  Y = testf(SG$design) #the design is $design, simple enough, right?
  logtheta <- rep(0,8)
  logtheta <- logthetaMLE(SG=SG,y=Y)
  
  n <- 50
  xp <- matrix(runif(d*10),n,8)
  SGpred <- SGGPpred(xp=xp, SG=SG, y=Y, logtheta = logtheta)
  
  my <- mean(Y)
  dy <- Y - my
  Sig <- CorrMat32Full(SG$design, theta=exp(logtheta))
  s <- CorrMat32Full2(xp, SG$design, theta = exp(logtheta))
  
  s2 <- c(t(Y) %*%solve(Sig, Y) / length(Y))
  exvar <- s2 * (1 - colSums(t(s) * solve(Sig, t(s))))
  sigmas <- (1/SGpred$var* exvar)
  if (any(abs(sigmas -sigmas[1])>sigmas[1]/100)) {
    plot(sigmas)
    print(sigmas)
    stop("Not all equal")
  }
  return(sigmas[1])
}
fff1(10)
ns <- c(10,10,10,10,20,20,20,20,30,40,50,60,70,80,90,100,120,150,200,300,400,500)
#ns <- c(10,10,10,10,20,20,20,20)
fs <- sapply(ns, fff1)
plot(ns, fs)
