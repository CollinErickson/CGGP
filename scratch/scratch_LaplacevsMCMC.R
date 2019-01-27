# Trying to check LaPlace approx vs MCMC
f <- function(x){x[1]*x[2]^2 + cos(2*pi*2*x[3]^2)}
d <- 4
nval <- 1e3
Xval <- matrix(runif(d*nval), ncol=d)
Yval <- apply(Xval, 1, f)
SG <- SGGPcreate(d=d, batchsize=1000)
y <- apply(SG$design, 1, f)
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


nlpPS1 <- apply(SG1$thetaPostSamples, 2, SGGP_internal_neglogpost, SG1, SG1$y)
nlpM1 <- SGGP_internal_neglogpost(SG1$thetaMAP, SG1, SG1$y)
nlpPS2 <- apply(SG2$thetaPostSamples, 2, SGGP_internal_neglogpost, SG2, SG2$y)
nlpM2 <- SGGP_internal_neglogpost(SG2$thetaMAP, SG2, SG2$y)
hist(nlpPS1, xlim=c(min(nlpPS1,nlpM1), max(nlpPS1[is.finite(nlpPS1)],nlpM1)),
     main="Hist of neglogpost from Laplace samples"); abline(v=nlpM1, col=2)
nlpPS1 / nlpM1
nlpPS1 / min(nlpPS1, nlpM1)
hist(nlpPS2, xlim=c(min(nlpPS2,nlpM2), max(nlpPS2,nlpM2)),
     main="Hist of neglogpost from MCMC samples"); abline(v=nlpM2, col=2)
nlpPS2 / nlpM2
nlpPS2 / min(nlpPS2, nlpM2)

table(SGGP_internal_neglogpost(SG1$thetaMAP, SG1, SG1$y) < apply(SG1$thetaPostSamples, 2, SGGP_internal_neglogpost, SG1, SG1$y))
table(SGGP_internal_neglogpost(SG2$thetaMAP, SG2, SG2$y) < apply(SG2$thetaPostSamples, 2, SGGP_internal_neglogpost, SG2, SG2$y))

stripchart(as.data.frame(t(SG1$thetaPostSamples)))
stripchart(as.data.frame(t(SG1$thetaMAP)), add=T, col=2, pch=19)

which.min(nlpPS2)
SG3 <- SGGPfit(SG2, Y=SG2$Y, theta0 = SG2$thetaPostSamples[,74])
SGGPvalstats(SG1, Xval, Yval)
SGGPvalstats(SG2, Xval, Yval)
SGGPvalstats(SG3, Xval, Yval)
