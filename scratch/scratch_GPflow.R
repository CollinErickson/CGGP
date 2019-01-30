library(reticulate)
env1 <- use_condaenv("myenv")
# gpflow <- import("gpflow")
# attributes(gpflow)
# repl_python()
# 
# 
# py$xx
# py$flowpred %>% str
# plot(py$xx, py$flowpred[[1]], type='l')
# points(py$xx, py$flowpred[[1]] + 2*sqrt(py$flowpred[[2]]), col=2, type='l')
# points(py$xx, py$flowpred[[1]] - 2*sqrt(py$flowpred[[2]]), col=2, type='l')
# points(py$xx, py$flowpred[[1]], type='l')
# points(py$X, py$Y, pch=19)
# py_eval("X")
# py$testvar <- 1243
# py_eval('testvar')


f <- function(x){x[1]^1.1+x[1]*x[2]^2 + log(x[1]+1)*sin(2*pi*x[3]^.9)^2}
SG <- SGGPcreate(d=3, batchsize=100)
y <- apply(SG$design, 1, f)
y <- cbind(y, y^1.1)
SG <- SGGPfit(SG=SG, Y=y)

X <- SG$design
# py$X <- X
# py$Y <- y#matrix(y, ncol=1)
xx <- matrix(runif(100*3), ncol=3)
# py$xx <- xx
yy <- apply(xx, 1, f)
yy <- cbind(yy, yy)
# py$yy <- yy
py_run_string("import gpflow
import numpy as np
import matplotlib

k = gpflow.kernels.Matern52(3, lengthscales=0.3)
m = gpflow.models.GPR(r.X, r.y, kern=k)
m
")
py_run_string("m.likelihood.variance = 0.0001
r.gpflowmean, r.gpflowvar = m.predict_y(r.xx)
")
# py$mean
# valstats(py$mean, py$var, yy)
valstats(gpflowmean, gpflowvar, yy)
SGGPvalstats(SG, xx, yy, bydim = F)
SGGPvalplot(SG, xx, yy)
SGGPprojectionplot(SG)
