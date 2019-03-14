f <- TestFunctions::wingweight
d <- 10
n <- 10000

system.time(
  sgwwc <- CGGPcreate(d, n)
)
# ^ 6 sec
system.time(
  sgwwf <- CGGPfit(sgwwc, apply(sgwwc$design, 1, f))
)
# ^ 35 seconds
# system.time(
#   sgwwfMCMC <- CGGPfit(sgwwc, apply(sgwwc$design, 1, f), laplaceapprox = F)
# )
# # ^ 1573 sec = 26.2 minutes

# # Do profvis on only 4 MCMC samples, should take about a minute
# sgwwc2 <- CGGPcreate(d, n)
# sgwwc2$numPostSamples <- 4
# pvmcmc <- profvis::profvis({
#   sgwwf <- CGGPfit(sgwwc2, apply(sgwwc2$design, 1, f), laplaceapprox = F)
# }, interval = .05)


profvis::profvis(
  {
    sgwwc <- CGGPcreate(d, n)
    sgwwf <- CGGPfit(sgwwc, apply(sgwwc$design, 1, f))
  }
)


profvis::profvis(
  {
    n <- 500
    sgwwc <- CGGPcreate(d, n)
    sgwwf <- CGGPfit(sgwwc, apply(sgwwc$design, 1, f))
    sgwwf <- CGGPappend(sgwwf, 500)
    sgwwf <- CGGPfit(sgwwf, apply(sgwwf$design, 1, f))
    sgwwf <- CGGPappend(sgwwf, 500)
    sgwwf <- CGGPfit(sgwwf, apply(sgwwf$design, 1, f))
    sgwwf <- CGGPappend(sgwwf, 500)
    sgwwf <- CGGPfit(sgwwf, apply(sgwwf$design, 1, f))
    sgwwf <- CGGPappend(sgwwf, 500)
    sgwwf <- CGGPfit(sgwwf, apply(sgwwf$design, 1, f))
  }
)



nt <- 1e3
xt <- matrix(runif(nt*d), ncol=d)
system.time(CGGPpred(sgwwf, xt))
system.time(CGGPpred(sgwwf, rbind(xt,xt+.0001)))

