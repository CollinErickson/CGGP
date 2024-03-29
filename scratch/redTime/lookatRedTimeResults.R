lo.in <- read.csv("../../../Desktop/redTimeData/LHS1L_n1000_s1225_all_output.csv")
lo <- lo.in[,-1]
lo <- as.matrix(lo)
set.seed(1225)
lx <- lhs::maximinLHS(n=1000, k=9)

sg8 <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_8039.rds")
sg1 <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_1319.rds")
sg3 <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_3119.rds")
SGGPvalplot(SGGP=sg8, Xval=lx, Yval=lo, d=3)


Ylog <- log(sg8$Y, 10)
sg8log <- SGGPfit(sg8, Y=Ylog)
lo.log <- log(lo, 10)
SGGPvalplot(sg8log, lx, lo.log, d=3)

ylhs8039 <- as.matrix(
  read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_all_output.csv")[,-1])
if (F) {
  set.seed(1226)
  xlhs8039 <- lhs::maximinLHS(n=8039, k=9)
  write.csv(xlhs8039, "../../../Desktop/redTimeData/LHS1L_n8039_s1226_matrix.csv")
} else {
  xlhs8039 <- read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_matrix.csv")
}
ylhs8039_100 <- ylhs8039[1:100,]
xlhs8039_100 <- xlhs8039[1:100,]

ylhs1000_100 <- lo[1:100,]
xlhs1000_100 <- lx[1:100,]
if (F) {
  mod.mlegp <- mlegp::mlegp(X=xlhs1000_100, Z=ylhs1000_100)
  saveRDS(mod.mlegp, "../../../Desktop/redTimeData/modmlegp100.rds")
  mod.mlegp.pca <- mlegp::mlegp(X=xlhs1000_100, Z=ylhs1000_100, PC.percent = 99.999)
  saveRDS(mod.mlegp.pca, "../../../Desktop/redTimeData/modmlegppca100.rds")
} else {
  mod.mlegp <- readRDS("../../../Desktop/redTimeData/modmlegp100.rds")
  mod.mlegp.pca <- readRDS("../../../Desktop/redTimeData/modmlegppca100.rds")
}
# mod.mlegp1 <- mlegp::mlegp(X=xlhs1000_100, Z=ylhs1000_100[,1])
# Try 300 pts
ylhs1000_300 <- lo[1:300,]
xlhs1000_300 <- lx[1:300,]
if (F) {
  mod.mlegp.300 <- mlegp::mlegp(X=xlhs1000_300, Z=ylhs1000_300)
  saveRDS(mod.mlegp.300, "../../../Desktop/redTimeData/modmlegp300.rds")
  mod.mlegp.pca.300 <- mlegp::mlegp(X=xlhs1000_300, Z=t(ylhs1000_300),
                                    PC.percent = 99.999)
  saveRDS(mod.mlegp.pca.300, "../../../Desktop/redTimeData/modmlegppca300.rds")
} else {
  mod.mlegp.300 <- readRDS("../../../Desktop/redTimeData/modmlegp300.rds")
  mod.mlegp.pca.300 <- readRDS("../../../Desktop/redTimeData/modmlegppca300.rds")
}

# pca.p1 <- mlegp::predict.gp(mod.mlegp.pca[[1]], lx)
# plot(lo[,1], pca.p1)
lapply(mod.mlegp.pca, function(mod) {predict(mod, lx)})
pmatpre <- sapply(1:mod.mlegp.pca$numGPs, function(i)predict(mod.mlegp.pca[[i]], lx))
pmat <- mod.mlegp.pca$UD %*% t(pmatpre)
mlegperrors <- as.matrix(lo) - t(pmat)

# ytest <- lo[101:1000,]
# xtest <- lx[101:1000,]
xtest <- xlhs8039
ytest <- ylhs8039
# sg1_d1 <- SGGPfit(SGGP = sg1, Y = sg1$Y[,1])
comp1D <- function(d) {
  pred.mlegp.pca <- mlegp::predict.gp(mod.mlegp.pca, xtest)
  pred.mlegp <- mlegp::predict.gp(mod.mlegp1, xtest)
  pred.sggp1_d1 <- predict(sg1_d1, xtest)
  pred.sggp1 <- predict(sg1, xtest)
  pred.sggp3 <- predict(sg3, xtest)
  pred.sggp8 <- predict(sg8, xtest)
  rmse.mlegp <- sqrt(mean((pred.mlegp - ytest[,1])^2))
  rmse.sggp1_1d <- sqrt(mean((pred.sggp1_d1$me[,1] - ytest[,1])^2))
  rmse.sggp1 <- sqrt(mean((pred.sggp1$me[,1] - ytest[,1])^2))
  rmse.sggp3 <- sqrt(mean((pred.sggp3$me[,1] - ytest[,1])^2))
  rmse.sggp8 <- sqrt(mean((pred.sggp8$me[,1] - ytest[,1])^2))
  c('mlegp'=rmse.mlegp, 'sggp1_d1'=rmse.sggp1_1d,
    'sggp1'=rmse.sggp1, 'sggp3'=rmse.sggp3, 'sggp8'=rmse.sggp8)
}
comp1D()

sapply(1:100, function(i) {forecast::BoxCox.lambda(lo[,i])})

compmods <- function() {
  # pred.mlegp.pca <- mlegp::predict.gp(mod.mlegp.pca, xtest)
  pred.mlegp.pca.pretrans <- sapply(1:mod.mlegp.pca$numGPs,
                                    function(i)predict(mod.mlegp.pca[[i]], xtest))
  pred.mlegp.pca <- t(mod.mlegp.pca$UD %*% t(pred.mlegp.pca.pretrans))
  # browser()
  pred.mlegp <- sapply(1:mod.mlegp$numGPs,
                       function(i)predict(mod.mlegp[[i]], xtest))
  pred.mlegp.300 <- sapply(1:mod.mlegp$numGPs,
                           function(i)predict(mod.mlegp.300[[i]], xtest))
  # pred.sggp1 <- predict(sg1, xtest)
  pred.sggp1 <- predict(sg1, xtest)
  pred.sggp3 <- predict(sg3, xtest)
  pred.sggp8 <- predict(sg8, xtest)
  rmse.mlegp.pca <- sqrt(mean((pred.mlegp.pca - ytest)^2))
  rmse.mlegp <- sqrt(mean((pred.mlegp - ytest)^2))
  rmse.mlegp.300 <- sqrt(mean((pred.mlegp.300 - ytest)^2))
  rmse.sggp1 <- sqrt(mean((pred.sggp1$me - ytest)^2))
  rmse.sggp1 <- sqrt(mean((pred.sggp1$me - ytest)^2))
  rmse.sggp3 <- sqrt(mean((pred.sggp3$me - ytest)^2))
  rmse.sggp8 <- sqrt(mean((pred.sggp8$me - ytest)^2))
  c('mlegp.pca'=rmse.mlegp.pca, 'mlegp'=rmse.mlegp,
    'mlegp.300'=rmse.mlegp.300,
    # 'sggp1_d1'=rmse.sggp1_1d, 
    'sggp1'=rmse.sggp1, 
    'sggp3'=rmse.sggp3, 'sggp8'=rmse.sggp8)
}
compmods()

# # Compare mlegp and SGGP
# mlegp.pca     mlegp mlegp.300     sggp1     sggp3     sggp8 
# 4999.9849  641.2593  144.7101 2108.9625  126.4526  106.1754 

# ^ PCA was done wrong, didn't use t(), was 100x100 so it didn't give error