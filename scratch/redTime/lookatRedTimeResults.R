lo.in <- read.csv("../../../Desktop/redTimeData/LHS1L_n1000_s1225_all_output.csv")
lo <- lo.in[,-1]
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

ylhs8039 <- read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_all_output.csv")
set.seed(1226)
xlhs8039 <- lhs::maximinLHS(n=8039, k=9)
ylhs8039_100 <- ylhs8039[1:100,]
xlhs8039_100 <- xlhs8039[1:100,]

ylhs1000_100 <- lo[1:100,]
xlhs1000_100 <- lx[1:100,]
mod.mlegp <- mlegp::mlegp(X=xlhs1000_100, Z=ylhs1000_100)
mod.mlegp.pca <- mlegp::mlegp(X=xlhs1000_100, Z=ylhs1000_100, PC.percent = 99.999)
mod.mlegp1 <- mlegp::mlegp(X=xlhs1000_100, Z=ylhs1000_100[,1])

ytest <- lo[101:1000,]
xtest <- lx[101:1000,]
sg1_d1 <- SGGPfit(SGGP = sg1, Y = sg1$Y[,1])
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
  c('mlegp'=rmse.mlegp, 'sggp1_d1'=rmse.sggp1_1d, 'sggp1'=rmse.sggp1, 'sggp3'=rmse.sggp3, 'sggp8'=rmse.sggp8)
}
comp1D()

sapply(1:100, function(i) {forecast::BoxCox.lambda(lo[,i])})