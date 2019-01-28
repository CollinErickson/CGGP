# Everything in this one will be on the log scale
redtimefolder <- "../../../Desktop/redTimeData/"

lo.in <- read.csv("../../../Desktop/redTimeData/LHS1L_n1000_s1225_all_output.csv")
lo <- lo.in[,-1]
lo <- log(as.matrix(lo))
set.seed(1225)
lx <- lhs::maximinLHS(n=1000, k=9)

# These were fit on non-log scale
sg8in <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_8039.rds")
sg1in <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_1319.rds")
sg3in <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_3119.rds")
# SGGPvalplot(SGGP=sg8, Xval=lx, Yval=lo, d=3)

# Refit on log scale, all four MV options:
#  with/out PCA, shared/unshared params
if (F) {
  sg1ps <- SGGPfit(sg1in, Y=sg1in$Y, use_PCA=T, separateoutputparameterdimensions=T)
  saveRDS(sg1ps, paste0(redtimefolder, "sg1ps.rds"))
  sg1po <- SGGPfit(sg1in, Y=sg1in$Y, use_PCA=T, separateoutputparameterdimensions=F)
  saveRDS(sg1po, paste0(redtimefolder, "sg1po.rds"))
  sg1ns <- SGGPfit(sg1in, Y=sg1in$Y, use_PCA=F, separateoutputparameterdimensions=T)
  saveRDS(sg1ns, paste0(redtimefolder, "sg1ns.rds"))
  sg1no <- SGGPfit(sg1in, Y=sg1in$Y, use_PCA=F, separateoutputparameterdimensions=F)
  saveRDS(sg1no, paste0(redtimefolder, "sg1no.rds"))
  
  sg3ps <- SGGPfit(sg3in, Y=sg3in$Y, use_PCA=T, separateoutputparameterdimensions=T)
  saveRDS(sg3ps, paste0(redtimefolder, "sg3ps.rds"))
  sg3po <- SGGPfit(sg3in, Y=sg3in$Y, use_PCA=T, separateoutputparameterdimensions=F)
  saveRDS(sg3po, paste0(redtimefolder, "sg3po.rds"))
  sg3ns <- SGGPfit(sg3in, Y=sg3in$Y, use_PCA=F, separateoutputparameterdimensions=T)
  saveRDS(sg3ns, paste0(redtimefolder, "sg3ns.rds"))
  sg3no <- SGGPfit(sg3in, Y=sg3in$Y, use_PCA=F, separateoutputparameterdimensions=F)
  saveRDS(sg3no, paste0(redtimefolder, "sg3no.rds"))
  
  sg8ps <- SGGPfit(sg8in, Y=sg8in$Y, use_PCA=T, separateoutputparameterdimensions=T)
  saveRDS(sg8ps, paste0(redtimefolder, "sg8ps.rds"))
  sg8po <- SGGPfit(sg8in, Y=sg8in$Y, use_PCA=T, separateoutputparameterdimensions=F)
  saveRDS(sg8po, paste0(redtimefolder, "sg8po.rds"))
  sg8ns <- SGGPfit(sg8in, Y=sg8in$Y, use_PCA=F, separateoutputparameterdimensions=T)
  saveRDS(sg8ns, paste0(redtimefolder, "sg8ns.rds"))
  sg8no <- SGGPfit(sg8in, Y=sg8in$Y, use_PCA=F, separateoutputparameterdimensions=F)
  saveRDS(sg8no, paste0(redtimefolder, "sg8no.rds"))
} else { # Just read in if already fit and saved
  sg1ps <- readRDS(paste0(redtimefolder, "sg1ps.rds"))
  sg1po <- readRDS(paste0(redtimefolder, "sg1po.rds"))
  sg1ns <- readRDS(paste0(redtimefolder, "sg1ns.rds"))
  sg1no <- readRDS(paste0(redtimefolder, "sg1no.rds"))
  sg3ps <- readRDS(paste0(redtimefolder, "sg3ps.rds"))
  sg3po <- readRDS(paste0(redtimefolder, "sg3po.rds"))
  sg3ns <- readRDS(paste0(redtimefolder, "sg3ns.rds"))
  sg3no <- readRDS(paste0(redtimefolder, "sg3no.rds"))
  sg8ps <- readRDS(paste0(redtimefolder, "sg8ps.rds"))
  sg8po <- readRDS(paste0(redtimefolder, "sg8po.rds"))
  sg8ns <- readRDS(paste0(redtimefolder, "sg8ns.rds"))
  sg8no <- readRDS(paste0(redtimefolder, "sg8no.rds"))
}

# Ylog <- log(sg8$Y, 10)
# sg8log <- SGGPfit(sg8, Y=Ylog)
# lo.log <- log(lo, 10)
# SGGPvalplot(sg8log, lx, lo.log, d=3)

ylhs8039 <- log(as.matrix(
  read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_all_output.csv")[,-1]))
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
  saveRDS(mod.mlegp, "../../../Desktop/redTimeData/modmlegp100_log.rds")
  mod.mlegp.pca <- mlegp::mlegp(X=xlhs1000_100, Z=t(ylhs1000_100), PC.percent = 99.999)
  saveRDS(mod.mlegp.pca, "../../../Desktop/redTimeData/modmlegppca100_log.rds")
} else {
  mod.mlegp <- readRDS("../../../Desktop/redTimeData/modmlegp100_log.rds")
  mod.mlegp.pca <- readRDS("../../../Desktop/redTimeData/modmlegppca100_log.rds")
}
# mod.mlegp1 <- mlegp::mlegp(X=xlhs1000_100, Z=ylhs1000_100[,1])
# Try 300 pts
ylhs1000_300 <- lo[1:300,]
xlhs1000_300 <- lx[1:300,]
if (F) {
  mod.mlegp.300 <- mlegp::mlegp(X=xlhs1000_300, Z=ylhs1000_300)
  saveRDS(mod.mlegp.300, "../../../Desktop/redTimeData/modmlegp300_log.rds")
  mod.mlegp.pca.300 <- mlegp::mlegp(X=xlhs1000_300, Z=t(ylhs1000_300),
                                    PC.percent = 99.999)
  saveRDS(mod.mlegp.pca.300, "../../../Desktop/redTimeData/modmlegppca300_log.rds")
} else {
  mod.mlegp.300 <- readRDS("../../../Desktop/redTimeData/modmlegp300_log.rds")
  mod.mlegp.pca.300 <- readRDS("../../../Desktop/redTimeData/modmlegppca300_log.rds")
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
# # sg1_d1 <- SGGPfit(SGGP = sg1, Y = sg1$Y[,1])
# comp1D <- function(d) {
#   pred.mlegp.pca <- mlegp::predict.gp(mod.mlegp.pca, xtest)
#   pred.mlegp <- mlegp::predict.gp(mod.mlegp1, xtest)
#   pred.sggp1_d1 <- predict(sg1_d1, xtest)
#   pred.sggp1 <- predict(sg1, xtest)
#   pred.sggp3 <- predict(sg3, xtest)
#   pred.sggp8 <- predict(sg8, xtest)
#   rmse.mlegp <- sqrt(mean((pred.mlegp - ytest[,1])^2))
#   rmse.sggp1_1d <- sqrt(mean((pred.sggp1_d1$me[,1] - ytest[,1])^2))
#   rmse.sggp1 <- sqrt(mean((pred.sggp1$me[,1] - ytest[,1])^2))
#   rmse.sggp3 <- sqrt(mean((pred.sggp3$me[,1] - ytest[,1])^2))
#   rmse.sggp8 <- sqrt(mean((pred.sggp8$me[,1] - ytest[,1])^2))
#   c('mlegp'=rmse.mlegp, 'sggp1_d1'=rmse.sggp1_1d,
#     'sggp1'=rmse.sggp1, 'sggp3'=rmse.sggp3, 'sggp8'=rmse.sggp8)
# }
# comp1D()

# sapply(1:100, function(i) {forecast::BoxCox.lambda(lo[,i])})

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
  # pred.sggp1 <- predict(sg1, xtest)
  # pred.sggp3 <- predict(sg3, xtest)
  # pred.sggp8 <- predict(sg8, xtest)
  # rmse.mlegp.pca <- sqrt(mean((pred.mlegp.pca - ytest)^2))
  # rmse.mlegp <- sqrt(mean((pred.mlegp - ytest)^2))
  # rmse.mlegp.300 <- sqrt(mean((pred.mlegp.300 - ytest)^2))
  # rmse.sggp1 <- sqrt(mean((pred.sggp1$me - ytest)^2))
  # rmse.sggp1 <- sqrt(mean((pred.sggp1$me - ytest)^2))
  # rmse.sggp3 <- sqrt(mean((pred.sggp3$me - ytest)^2))
  # rmse.sggp8 <- sqrt(mean((pred.sggp8$me - ytest)^2))
  browser()
  mlegp.100 <- valstats(pred.mlegp$mean, pred.mlegp$var, ytest)
  sg1ps.stats <- SGGPvalstats(sg1ps, xtest, ytest)
  sg1po.stats <- SGGPvalstats(sg1po, xtest, ytest)
  sg1ns.stats <- SGGPvalstats(sg1ns, xtest, ytest)
  sg1no.stats <- SGGPvalstats(sg1no, xtest, ytest)
  sg3ps.stats <- SGGPvalstats(sg3ps, xtest, ytest)
  sg3po.stats <- SGGPvalstats(sg3po, xtest, ytest)
  sg3ns.stats <- SGGPvalstats(sg3ns, xtest, ytest)
  sg3no.stats <- SGGPvalstats(sg3no, xtest, ytest)
  sg8ps.stats <- SGGPvalstats(sg8ps, xtest, ytest)
  sg8po.stats <- SGGPvalstats(sg8po, xtest, ytest)
  sg8ns.stats <- SGGPvalstats(sg8ns, xtest, ytest)
  sg8no.stats <- SGGPvalstats(sg8no, xtest, ytest)
  # c('mlegp.pca'=rmse.mlegp.pca, 'mlegp'=rmse.mlegp,
  #   'mlegp.300'=rmse.mlegp.300,
    # 'sggp1_d1'=rmse.sggp1_1d, 
    # 'sggp1'=rmse.sggp1, 
    # 'sggp3'=rmse.sggp3, 'sggp8'=rmse.sggp8)
  browser()
  rbind(
    sg1ps.stats, sg1po.stats, sg1ns.stats, sg1no.stats,
    sg3ps.stats, sg3po.stats, sg3ns.stats, sg3no.stats,
    sg8ps.stats, sg8po.stats, sg8ns.stats, sg8no.stats
  )
}
compmods()

# # Compare mlegp and SGGP
# mlegp.pca     mlegp mlegp.300     sggp1     sggp3     sggp8 
# 4999.9849  641.2593  144.7101 2108.9625  126.4526  106.1754 

# ^ PCA was done wrong, didn't use t(), was 100x100 so it didn't give error