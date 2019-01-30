f <- TestFunctions::wingweight
d <- 10
n <- 10000

system.time(
  sgwwc <- SGGPcreate(d, n)
)
# ^ 6 sec
system.time(
  sgwwf <- SGGPfit(sgwwc, apply(sgwwc$design, 1, f))
)
# ^ 35 seconds
system.time(
  sgwwfMCMC <- SGGPfit(sgwwc, apply(sgwwc$design, 1, f), laplaceapprox = F)
)
# ^ 1573 sec = 26.2 minutes

# Do profvis on only 4 MCMC samples, should take about a minute
sgwwc2 <- SGGPcreate(d, n)
sgwwc2$numPostSamples <- 4
pvmcmc <- profvis::profvis({
  sgwwf <- SGGPfit(sgwwc2, apply(sgwwc2$design, 1, f), laplaceapprox = F)
}, interval = .05)


profvis::profvis(
  {
    sgwwc <- SGGPcreate(d, n)
    sgwwf <- SGGPfit(sgwwc, apply(sgwwc$design, 1, f))
  }
)



nt <- 1e3
xt <- matrix(runif(nt*d), ncol=d)
system.time(SGGPpred(sgwwf, xt))
system.time(SGGPpred(sgwwf, rbind(xt,xt+.0001)))
