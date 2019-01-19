# Trying to check LaPlace approx vs MCMC

SG <- SGGPcreate(d=3, batchsize=1000)
y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
SG1 <- SGGPfit(SG=SG, Y=y, laplaceapprox = T)
SG2 <- SGGPfit(SG=SG, Y=y, laplaceapprox = F)


stripchart(data.frame(t(SG1$thetaPostSamples)), main="SG1 1k LaPlace")
stripchart(data.frame(t(SG2$thetaPostSamples)), main="SG1 1k MCMC")


SG10k <- SGGPcreate(d=4, batchsize=5000)
y <- apply(SG10k$design, 1, function(x){x[1]^1.1+x[2]^2 + .3*sin(2*pi*x[3]*4) + .6*cos(2*pi*x[4]^2*5)})
SG10k1 <- SGGPfit(SG=SG10k, Y=y, laplaceapprox = T)
SG10k2 <- SGGPfit(SG=SG10k, Y=y, laplaceapprox = F)


stripchart(data.frame(t(SG10k1$thetaPostSamples)), main="SG1 5k LaPlace")
stripchart(data.frame(t(SG10k2$thetaPostSamples)), main="SG1 5k MCMC")

SG10k1$thetaMAP
