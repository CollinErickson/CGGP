
f <- TestFunctions::borehole
d <- 8
set.seed(0)
nval <- 1e4
Xval <- matrix(runif(nval*d), nval, d)
Yval <- apply(Xval, 1, f)

# n0 <- 64
corr <- SGGP_internal_CorrMatCauchySQT
sel.method <- "TS"
batchsize <- 64
pred.fullBayes <- TRUE
append.rimseperpoint <- TRUE
use_laplaceapprox <- TRUE
grid_size <- c(1, 2, 2, 2, 4, 4, 4, 4, 4, 6, 32)




is_power_of_2 <- function(x) {
  abs(log(x, 2) - floor(log(x, 2)+.5)) < 1e-8
}

start.time <- Sys.time()

sg <- SGGPcreate(d, batchsize, corr=corr)
y <- apply(sg$design, 1, f)
sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox)
nallotted <- batchsize

sg.stats <- NULL

while(nallotted < 1024) {
  sg <- SGGPappend(sg, batchsize=batchsize, selectionmethod=sel.method, RIMSEperpoint = append.rimseperpoint)
  y <- apply(sg$design, 1, f)
  sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox)
  nallotted <- nallotted + batchsize
  if (nallotted >= 512 && is_power_of_2(nallotted)) {
    newstats <- SGGPvalstats(sg, Xval = Xval, Yval = Yval)
    newstats$nallotted <- nallotted
    newstats$nused <- nrow(sg$design)
    newstats$elapsedtime <- as.numeric(Sys.time() - start.time, units='secs')
    sg.stats <- if (is.null(sg.stats)) newstats else rbind(sg.stats, newstats)
  }
}
while(nallotted < 16384) {
  sg <- SGGPappend(sg, batchsize=512, selectionmethod=sel.method, RIMSEperpoint = append.rimseperpoint)
  y <- apply(sg$design, 1, f)
  sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox)
  nallotted <- nallotted + 512
  if (nallotted >= 512 && is_power_of_2(nallotted)) {
    newstats <- SGGPvalstats(sg, Xval = Xval, Yval = Yval)
    newstats$nallotted <- nallotted
    newstats$nused <- nrow(sg$design)
    newstats$elapsedtime <- as.numeric(Sys.time() - start.time, units='secs')
    sg.stats <- if (is.null(sg.stats)) newstats else rbind(sg.stats, newstats)
    
  }
}
