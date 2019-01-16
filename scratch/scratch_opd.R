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
set.seed(1)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
y <- cbind(y1, y2)
SG <- SGGPfit(SG, Y=y, use_PCA = F)
yMVpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(yMVpred[,1], y1, 1e-4)
expect_equal(yMVpred[,2], y2, 1e-4)

set.seed(1)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
SG <- SGGPfit(SG, Y=y1, use_PCA = F)
y1Vpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(y1Vpred[,1], y1, 1e-4)


# 3. MV with PCA, separate output par dims
set.seed(0)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
y <- cbind(y1, y2)
set.seed(1)
SG <- SGGPfit(SG, Y=y, separateoutputparameterdimensions = T)
yMVpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(yMVpred[,1], y1, 1e-4)
expect_equal(yMVpred[,2], y2, 1e-4)
SGGPappend(SG, 30, selectionmethod = "Greedy")
SGGPappend(SG, 30, selectionmethod = "UCB")
SGGPappend(SG, 30, selectionmethod = "TS")

set.seed(0)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
set.seed(1)
SG <- SGGPfit(SG, Y=y1)
SG$thetaMAP
SG$pw
y1Vpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(y1Vpred[,1], y1, 1e-4)


# 4. MV without PCA, separate output par dims
set.seed(0)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
y <- cbind(y1, y2)
set.seed(1)
SG <- SGGPfit(SG, Y=y, separateoutputparameterdimensions = T, use_PCA = F)
yMVpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(yMVpred[,1], y1, 1e-4)
expect_equal(yMVpred[,2], y2, 1e-4)
SGGPappend(SG, 30, selectionmethod = "Greedy")
SGGPappend(SG, 30, selectionmethod = "UCB")
SGGPappend(SG, 30, selectionmethod = "TS")

set.seed(0)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
set.seed(1)
SG <- SGGPfit(SG, Y=y1)
SG$thetaMAP
SG$pw
y1Vpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(y1Vpred[,1], y1, 1e-4)






# ---------------------------
#    With supplementary data
# -----------------------------



# 3 w/ SUPP. MV with PCA, separate output par dims
set.seed(0)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
y <- cbind(y1, y2)
set.seed(0)
nsup <- 20
xsup <- matrix(runif(nsup*3), nsup, 3)
ysup1 <- apply(xsup, 1, f1)
ysup2 <- apply(xsup, 1, f2)
ysup <- cbind(ysup1, ysup2)
set.seed(1)
SG <- SGGPfit(SG, Y=y, separateoutputparameterdimensions = T, Xs=xsup, Ys=ysup)
yMVpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(yMVpred[,1], y1, 1e-4)
expect_equal(yMVpred[,2], y2, 1e-4)
yMVpredsup <- SGGPpred(xsup, SG=SG)$mean
expect_equal(yMVpredsup[,1], ysup1, 1e-4)
expect_equal(yMVpredsup[,2], ysup2, 1e-4)
SGGPappend(SG, 30, selectionmethod = "Greedy")
SGGPappend(SG, 30, selectionmethod = "UCB")
SGGPappend(SG, 30, selectionmethod = "TS")



# 4. MV without PCA, separate output par dims
set.seed(0)
SG <- SGGPcreate(d=3, batchsize=100)
y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
y <- cbind(y1, y2)
set.seed(0)
nsup <- 20
xsup <- matrix(runif(nsup*3), nsup, 3)
ysup1 <- apply(xsup, 1, f1)
ysup2 <- apply(xsup, 1, f2)
ysup <- cbind(ysup1, ysup2)
set.seed(1)
SG <- SGGPfit(SG, Y=y, Xs=xsup, Ys=ysup,
              separateoutputparameterdimensions = T, use_PCA = F)
yMVpred <- SGGPpred(SG$design, SG=SG)$mean
expect_equal(yMVpred[,1], y1, 1e-4)
expect_equal(yMVpred[,2], y2, 1e-4)
yMVpredsup <- SGGPpred(xsup, SG=SG)$mean
expect_equal(yMVpredsup[,1], ysup1, 1e-4)
expect_equal(yMVpredsup[,2], ysup2, 1e-4)
SGGPappend(SG, 30, selectionmethod = "Greedy")
SGGPappend(SG, 30, selectionmethod = "UCB")
SGGPappend(SG, 30, selectionmethod = "TS")
