set.seed(0)

f1 <- function(x){x[1]+x[2]^2}
f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}

# First check MV with PCA
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
y <- cbind(y1, y2)
SG <- SGGPfit(SG, Y=y)
yMVpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(yMVpred[,1], y1, 1e-4)
expect_equal(yMVpred[,2], y2, 1e-4)

# 2. MV without PCA
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
SG <- SGGPfit(SG, Y=y, use_PCA = F)
yMVpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(yMVpred[,1], y1, 1e-4)
expect_equal(yMVpred[,2], y2, 1e-4)

# 3. MV with PCA, separate output par dims
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
y <- cbind(y1, y2)
SG <- SGGPfit(SG, Y=y, separateoutputparameterdimensions = T)
yMVpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(yMVpred[,1], y1, 1e-4)
expect_equal(yMVpred[,2], y2, 1e-4)

set.seed(0)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
SG <- SGGPfit(SG, Y=y1)
SG$thetaMAP
SG$pw
y1Vpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(y1Vpred[,1], y1, 1e-4)
