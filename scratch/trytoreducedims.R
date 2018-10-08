SGGPprednull <- function(xp,SG, y, ..., logtheta, theta) {
  nulldims <- which(logtheta>=1.999)
  print(paste('nulldims are', nulldims))
  SGnull <- list()
  for (name in names(SG)) {
    # print(SG[name])
    # print(SGnull)
    SGnull[[name]] <- SG[[name]]
  }
  #return(SGnull)
  SGnull$xmin=SG$xmin[-nulldims]
  SGnull$xmax=SG$xmax[-nulldims]
  SGnull$uo=SG$uo[,-nulldims]
  SGnull$po=SG$po[,-nulldims]
  SGnull$gridsizes=SG$gridsizes[,-nulldims]
  SGnull$gridsizest=SG$gridsizest[,-nulldims]
  SGnull$design=SG$design[,-nulldims]
  SGnull$d=SG$d - length(nulldims)
  #   
  #   =SG$[,-nulldims],
  #   =SG$[,-nulldims],
  #   =SG$[,-nulldims],
  #   =SG$[,-nulldims],
  #   =SG$[,-nulldims],
  #   =SG$[,-nulldims],
  #SGnull
  print(str(SGnull))
  SGGPpred(xp=xp, SG=SGnull, y=y, logtheta=logtheta)
}
# This doesn't work. If you remove dimensions, then you have duplicated rows.
# So you have multiple (1,1,..,1,1) rows, but the code is set up so only the first row is that.
# Maybe can fix by removing those, or averagin them.
SGGPprednull(SG$design[di,],SG,Y,logtheta=pmin(logthetaest,2))
