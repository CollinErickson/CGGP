set.seed(0)
d <- 3
SG <- SGGPcreate(d=d, batchsize=500)
f <- function(x){x[1]+x[2]^2+sin(x[3])^3}
y <- apply(SG$design, 1, f)

nsup <- 400
xsup <- matrix(runif(nsup*3), nsup, 3)
ysup <- apply(xsup, 1, f)

SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup)

SG <- SGGPfit(SG, Y=y)

# Check neglogpost grad matches gneglogpost
theta <- SG$thetaMAP / 2 # Don't want values near -1 or +1
epsval <- 1e-4
thetagrad <- SGGP_internal_gneglogpost(theta, SG, SG$y)
for (i in 1:length(SG$thetaMAP)) {
  eps <- rep(0, length(SG$thetaMAP))
  eps[i] = eps[i] + epsval
  numgrad <- (SGGP_internal_neglogpost(theta + eps, SG, SG$y) - 
                SGGP_internal_neglogpost(theta - eps, SG, SG$y)) / (2*epsval)
  print(cbind(thetagrad[i], numgrad))
}

# Works with supplementary data
nsup <- 130
xsup <- matrix(runif(nsup*3), nsup, 3)
ysup <- apply(xsup, 1, f)
SG <- SGGPfit(SG, Y=SG$Y, Xs=xsup, Ys=ysup)

# Check neglogpost grad matches gneglogpost for all HandlingSuppData options
theta <- SG$thetaMAP / 2 # Don't want values near -1 or +1
epsval <- 1e-4
for (handling in c("Correct", "Only", "Ignore", "Mixture", "MarginalValidation", "FullValidation")) {
  thetagrad <- SGGP_internal_gneglogpost(theta, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)
  numgrad <- rep(0, length(SG$thetaMAP))
  for (i in 1:length(SG$thetaMAP)) {
    eps <- rep(0, length(SG$thetaMAP))
    eps[i] = eps[i] + epsval
    # numgrad <- (SGGP_internal_neglogpost(theta + eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) - 
    #               SGGP_internal_neglogpost(theta - eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)) / (2*epsval)
    numgrad[i] <- (-SGGP_internal_neglogpost(theta + 2*eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) + 
                     8*SGGP_internal_neglogpost(theta + eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) - 
                     8*SGGP_internal_neglogpost(theta - eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) + 
                     SGGP_internal_neglogpost(theta - 2*eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)) / (12*epsval)
  }
  # expect_equal(c(thetagrad), numgrad, tol=1e-4, info = handling)
  print(handling)
  print(rbind(gneglogpost=c(thetagrad), numgrad=numgrad))
}