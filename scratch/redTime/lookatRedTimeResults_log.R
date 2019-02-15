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
  sg1ps <- SGGPfit(sg1in, Y=log(sg1in$Y), use_PCA=T, separateoutputparameterdimensions=T)
  saveRDS(sg1ps, paste0(redtimefolder, "sg1ps.rds"))
  sg1po <- SGGPfit(sg1in, Y=log(sg1in$Y), use_PCA=T, separateoutputparameterdimensions=F)
  saveRDS(sg1po, paste0(redtimefolder, "sg1po.rds"))
  sg1ns <- SGGPfit(sg1in, Y=log(sg1in$Y), use_PCA=F, separateoutputparameterdimensions=T)
  saveRDS(sg1ns, paste0(redtimefolder, "sg1ns.rds"))
  sg1no <- SGGPfit(sg1in, Y=log(sg1in$Y), use_PCA=F, separateoutputparameterdimensions=F)
  saveRDS(sg1no, paste0(redtimefolder, "sg1no.rds"))
  
  sg3ps <- SGGPfit(sg3in, Y=log(sg3in$Y), use_PCA=T, separateoutputparameterdimensions=T)
  saveRDS(sg3ps, paste0(redtimefolder, "sg3ps.rds"))
  sg3po <- SGGPfit(sg3in, Y=log(sg3in$Y), use_PCA=T, separateoutputparameterdimensions=F)
  saveRDS(sg3po, paste0(redtimefolder, "sg3po.rds"))
  sg3ns <- SGGPfit(sg3in, Y=log(sg3in$Y), use_PCA=F, separateoutputparameterdimensions=T)
  saveRDS(sg3ns, paste0(redtimefolder, "sg3ns.rds"))
  sg3no <- SGGPfit(sg3in, Y=log(sg3in$Y), use_PCA=F, separateoutputparameterdimensions=F)
  saveRDS(sg3no, paste0(redtimefolder, "sg3no.rds"))
  
  sg8ps <- SGGPfit(sg8in, Y=log(sg8in$Y), use_PCA=T, separateoutputparameterdimensions=T)
  saveRDS(sg8ps, paste0(redtimefolder, "sg8ps.rds"))
  sg8po <- SGGPfit(sg8in, Y=log(sg8in$Y), use_PCA=T, separateoutputparameterdimensions=F)
  saveRDS(sg8po, paste0(redtimefolder, "sg8po.rds"))
  sg8ns <- SGGPfit(sg8in, Y=log(sg8in$Y), use_PCA=F, separateoutputparameterdimensions=T)
  saveRDS(sg8ns, paste0(redtimefolder, "sg8ns.rds"))
  sg8no <- SGGPfit(sg8in, Y=log(sg8in$Y), use_PCA=F, separateoutputparameterdimensions=F)
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

# SGGP from Test2, was run on log scale. Used PCA and single theta
sgT2_227 <- readRDS(paste0(redtimefolder, "out_T2_SGGP-227.rds"))
sgT2_455 <- readRDS(paste0(redtimefolder, "out_T2_SGGP-455.rds"))
sgT2_1063 <- readRDS(paste0(redtimefolder, "out_T2_SGGP-1063.rds"))
sgT3_227 <- readRDS(paste0(redtimefolder, "out_T3_SGGP-227.rds"))
sgT3_455 <- readRDS(paste0(redtimefolder, "out_T3_SGGP-455.rds"))
sgT3_1061 <- readRDS(paste0(redtimefolder, "out_T3_SGGP-1061.rds"))
sgT3_1669 <- readRDS(paste0(redtimefolder, "out_T3_SGGP-1669.rds"))
sgT4_227 <- readRDS(paste0(redtimefolder, "out_T4_SGGP-227.rds"))
sgT4_455 <- readRDS(paste0(redtimefolder, "out_T4_SGGP-455.rds"))
sgT4_1061 <- readRDS(paste0(redtimefolder, "out_T4_SGGP-1061.rds"))
sgT4_1517 <- readRDS(paste0(redtimefolder, "out_T4_SGGP-1517.rds"))
sgS1_227 <- readRDS(paste0(redtimefolder, "out_Sup1_SGGP-227.rds"))
sgS1_453 <- readRDS(paste0(redtimefolder, "out_Sup1_SGGP-453.rds"))
sgO90_227 <- readRDS(paste0(redtimefolder, "out_O90_SGGP-227.rds"))
sgO90_455 <- readRDS(paste0(redtimefolder, "out_O90_SGGP-455.rds"))
sgO50_227 <- readRDS(paste0(redtimefolder, "out_O50_SGGP-227.rds"))
sgO50_455 <- readRDS(paste0(redtimefolder, "out_O50_SGGP-455.rds"))
sgO50_1061 <- readRDS(paste0(redtimefolder, "out_O50_SGGP-1061.rds"))
sgO50_2051 <- readRDS(paste0(redtimefolder, "out_O50_SGGP-2051.rds"))


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
  xlhs8039 <- unname(as.matrix(read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_matrix.csv")[,-1]))
}
ylhs8039_100 <- ylhs8039[1:100,]
xlhs8039_100 <- xlhs8039[1:100,]

ylhs1000_100 <- lo[1:100,]
xlhs1000_100 <- lx[1:100,]
library(mlegp)
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
# lapply(mod.mlegp.pca, function(mod) {predict(mod, lx)})
# pmatpre <- sapply(1:mod.mlegp.pca$numGPs, function(i)predict(mod.mlegp.pca[[i]], lx))
# pmat <- mod.mlegp.pca$UD %*% t(pmatpre)
# mlegperrors <- as.matrix(lo) - t(pmat)

# ytest <- lo[101:1000,]
# xtest <- lx[101:1000,]
xtest <- xlhs8039[1:1000,]
ytest <- ylhs8039[1:1000,]
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

# compmods <- function() {#browser()
# pred.mlegp.pca <- mlegp::predict.gp(mod.mlegp.pca, xtest)
pred.mlegp.pca.pretrans <- lapply(1:mod.mlegp.pca$numGPs,
                                  function(i)predict(mod.mlegp.pca[[i]], xtest, se.fit=T))
pred.mlegp.pca.pretrans.means <- do.call(cbind, lapply(pred.mlegp.pca.pretrans, function(x) x$fit))
pred.mlegp.pca.pretrans.ses <- do.call(cbind, lapply(pred.mlegp.pca.pretrans, function(x) x$se))
pred.mlegp.pca <- list(
  mean = t(mod.mlegp.pca$UD %*% t(pred.mlegp.pca.pretrans.means)),
  var = (t(mod.mlegp.pca$UD %*% t(pred.mlegp.pca.pretrans.ses)))^2)
# browser()
pred.mlegp.raw <- lapply(1:mod.mlegp$numGPs,
                         function(i)predict(mod.mlegp[[i]], xtest, se.fit=T))
pred.mlegp <- list(
  mean = do.call(cbind, lapply(pred.mlegp.raw, function(x) x$fit)),
  var = do.call(cbind, lapply(pred.mlegp.raw, function(x) x$se^2))
)

pred.mlegp.pca.300.pretrans <- lapply(1:mod.mlegp.pca.300$numGPs,
                                      function(i)predict(mod.mlegp.pca.300[[i]], xtest, se.fit=T))
pred.mlegp.pca.300.pretrans.means <- do.call(cbind, lapply(pred.mlegp.pca.300.pretrans, function(x) x$fit))
pred.mlegp.pca.300.pretrans.ses <- do.call(cbind, lapply(pred.mlegp.pca.300.pretrans, function(x) x$se))
pred.mlegp.pca.300 <- list(
  mean = t(mod.mlegp.pca.300$UD %*% t(pred.mlegp.pca.300.pretrans.means)),
  var = (t(mod.mlegp.pca.300$UD %*% t(pred.mlegp.pca.300.pretrans.ses)))^2)

pred.mlegp.300.raw <- lapply(1:mod.mlegp$numGPs,
                             function(i)predict(mod.mlegp.300[[i]], xtest, se.fit=T))
pred.mlegp.300 <- list(
  mean = do.call(cbind, lapply(pred.mlegp.300.raw, function(x) x$fit)),
  var = do.call(cbind, lapply(pred.mlegp.300.raw, function(x) x$se^2))
)

mlegp.pca.100.stats <- valstats(pred.mlegp.pca$mean, pred.mlegp.pca$var, ytest)
mlegp.100.stats <- valstats(pred.mlegp$mean, pmax(1e-8, pred.mlegp$var), ytest)
mlegp.pca.300.stats <- valstats(pred.mlegp.pca.300$mean, pmax(1e-8, pred.mlegp.pca.300$var), ytest)
mlegp.300.stats <- valstats(pred.mlegp.300$mean, pmax(1e-8, pred.mlegp.300$var), ytest)
sg1ps.stats <- SGGPvalstats(sg1ps, xtest, ytest, bydim=F)
sg1po.stats <- SGGPvalstats(sg1po, xtest, ytest, bydim=F)
sg1ns.stats <- SGGPvalstats(sg1ns, xtest, ytest, bydim=F)
sg1no.stats <- SGGPvalstats(sg1no, xtest, ytest, bydim=F)
sg3ps.stats <- SGGPvalstats(sg3ps, xtest, ytest, bydim=F)
sg3po.stats <- SGGPvalstats(sg3po, xtest, ytest, bydim=F)
sg3ns.stats <- SGGPvalstats(sg3ns, xtest, ytest, bydim=F)
sg3no.stats <- SGGPvalstats(sg3no, xtest, ytest, bydim=F)
sg8ps.stats <- SGGPvalstats(sg8ps, xtest, ytest, bydim=F)
sg8po.stats <- SGGPvalstats(sg8po, xtest, ytest, bydim=F)
sg8ns.stats <- SGGPvalstats(sg8ns, xtest, ytest, bydim=F)
sg8no.stats <- SGGPvalstats(sg8no, xtest, ytest, bydim=F)
sgT2_227.stats <- SGGPvalstats(sgT2_227, xtest, ytest, bydim=F)
sgT2_455.stats <- SGGPvalstats(sgT2_455, xtest, ytest, bydim=F)
sgT2_1063.stats <- SGGPvalstats(sgT2_1063, xtest, ytest, bydim=F)
sgT3_227.stats <- SGGPvalstats(sgT3_227, xtest, ytest, bydim=F)
sgT3_455.stats <- SGGPvalstats(sgT3_455, xtest, ytest, bydim=F)
sgT3_1061.stats <- SGGPvalstats(sgT3_1061, xtest, ytest, bydim=F)
sgT3_1669.stats <- SGGPvalstats(sgT3_1669, xtest, ytest, bydim=F)
sgT4_227.stats <- SGGPvalstats(sgT4_227, xtest, ytest, bydim=F)
sgT4_455.stats <- SGGPvalstats(sgT4_455, xtest, ytest, bydim=F)
sgT4_1061.stats <- SGGPvalstats(sgT4_1061, xtest, ytest, bydim=F)
sgT4_1517.stats <- SGGPvalstats(sgT4_1517, xtest, ytest, bydim=F)
sgS1_227.stats <- SGGPvalstats(sgS1_227, xtest, ytest, bydim=F)
sgS1_453.stats <- SGGPvalstats(sgS1_453, xtest, ytest, bydim=F)
sgO90_227.stats <- SGGPvalstats(sgO90_227, xtest, ytest[,90])
sgO90_455.stats <- SGGPvalstats(sgO90_455, xtest, ytest[,90])
sgO50_227.stats <- SGGPvalstats(sgO50_227, xtest, ytest[,50])
sgO50_455.stats <- SGGPvalstats(sgO50_455, xtest, ytest[,50])
sgO50_1061.stats <- SGGPvalstats(sgO50_1061, xtest, ytest[,50])
sgO50_2051.stats <- SGGPvalstats(sgO50_2051, xtest, ytest[,50])


rbind(
  mlegp.pca.100=mlegp.pca.100.stats,
  mlegp.100=mlegp.100.stats,
  mlegp.pca.300=mlegp.pca.300.stats,
  mlegp.300=mlegp.300.stats,
  sg1ps=sg1ps.stats,
  sg1po=sg1po.stats,
  sg1ns=sg1ns.stats,
  sg1no=sg1no.stats,
  sg3ps=sg3ps.stats,
  sg3po=sg3po.stats,
  sg3ns=sg3ns.stats,
  sg3no=sg3no.stats,
  sg8ps=sg8ps.stats,
  sg8po=sg8po.stats,
  sg8ns=sg8ns.stats,
  sg8no=sg8no.stats
)
# }
# compmodout <- compmods()

# Compmods on 101 test points
#                   RMSE        score    CRPscore  coverage
# mlegp.pca.100 0.0034514503  -8.383659 0.001930322 0.7356436
# mlegp.100     0.0023564587        NaN         NaN 0.6925743
# mlegp.pca.300 0.0026704462        NaN         NaN 0.3819802
# mlegp.300     0.0007643889        NaN         NaN 0.5571287
# sg1ps         0.0115567585        NaN         NaN        NA
# sg1po         0.2333113866  -1.755597 0.135213044 0.9997030
# sg1ns         0.0111944670  -1.678100 0.146692876 1.0000000
# sg1no         0.0113708078  -5.194722 0.018646834 1.0000000
# sg3ps         0.0091577089        NaN         NaN        NA
# sg3po         0.0064022953  -6.951926 0.007851538 1.0000000
# sg3ns         0.0089211176  -2.405820 0.091034747 0.9997030
# sg3no         0.0050380924  -8.267539 0.004245088 1.0000000
# sg8ps         0.0024128541        NaN         NaN        NA
# sg8po         0.0033162834  -8.231127 0.004110924 1.0000000
# sg8ns         0.0029106681 -10.864965 0.001490032 0.9880198
# sg8no         0.0031820504 -10.700086 0.001578966 0.9898020

# Compmods on 1000 test points

#               RMSE         score     CRPscore     coverage corr   R2 RMSEnorm
# mlegp.pca.100 0.006266007  -8.395314 0.003103312  0.75365
# mlegp.100     0.004634059        NaN         NaN  0.71058
# mlegp.pca.300 0.003455452        NaN         NaN  0.39905
# mlegp.300     0.001589360        NaN         NaN  0.56924
# sg1ps         0.014889951  -2.827439 0.059072836  1.00000
# sg1po         0.276976502  -1.639099 0.153035155  0.99591
# sg1ns         0.014714046  -1.609238 0.152545708  1.00000
# sg1no         0.015580332  -5.166784 0.019306609  1.00000
# sg3ps         0.010514521  -2.984185 0.053966859  1.00000
# sg3po         0.009217553  -6.815351 0.008686614  0.99999
# sg3ns         0.010226738  -2.377687 0.091978759  0.99920
# sg3no         0.006792111  -8.111536 0.004778535  0.99880
# sg8ps         0.003774826 -10.667542 0.001689869  0.96988
# sg8po         0.004603387  -8.071398 0.004567326  0.99896
# sg8ns         0.004269050 -10.467252 0.001873210  0.96632
# sg8no         0.004635364 -10.252414 0.001998025  0.96958
# sgT2_227      0.043509    -4.712862  0.02779118        1
# SGT2_455      0.01959066  -5.633467  0.01592225        1
# SGT2_1063     0.01213821  -6.330312  0.01097829        1
# SGT3_227      0.009422946 -8.570286  0.004714287  0.97654 0.9999836            0.01404303
# SGT3_455      0.004432364 -9.997243  0.002174571  0.99006 0.9999963            0.006551142
# SGT3_1061     0.002942663 -10.90342  0.001363505  0.98783 0.9999984            0.004319948
# SGT3_1669     0.001346887 -11.25519  0.0009854101 0.99982 0.9999997            0.002063127
# SGT4_227      0.007950075 -8.654136  0.004338907  0.98703 0.9999878 0.9999757  0.01223752
# SGT4_455      0.00443485  -10.02579  0.002146854  0.98864 0.9999964 0.9999924  0.006340489
# SGT4_1061     0.002773989 -10.90239  0.001325069   0.9913 0.9999985 0.999997   0.00404263
# SGT4_1517     0.001605968 -11.07487  0.001113512  0.99981 0.9999995 0.999999   0.002448991
# SGS1_255+100  0.005890464 -9.658773  0.003019816  0.92948 0.9999933 0.9999866  0.008979054
# SGS1_453+100  0.004320681 -10.34532  0.002162065  0.91754 0.9999965 0.9999928  0.006389303


# Single output
require(mlegp)
mlegp.O90.100 <- mlegp(xlhs1000_100, ylhs1000_100[,90])
pred.mlegp.O90.100 <- predict(mlegp.O90.100, xtest, se=T)
pred.sgO90_227 <- predict(sgO90_227, xtest)
pred.sgO90_455 <- predict(sgO90_455, xtest)
valstats(c(pred.mlegp.O90.100$f), c(pred.mlegp.O90.100$s), ytest[,90])
plot(ytest[,90], pred.sgO90_227$mean - ytest[,90])
points(ytest[,90], pred.sgO90_455$mean - ytest[,90], col=2)
points(ytest[,90], pred.mlegp.O90.100$fi-ytest[,90], col=3, pch=4)
# mlegp.O90.100  0.004768259   NaN      NaN    0.731 0.9999808 0.999961
# sgO90_227      0.008338798 -8.663197 0.004349465    0.993 0.9999416 0.9998807
# sgO90_455      0.01031961 -8.749555 0.004792865    0.961 0.9999106 0.9998173


mlegp.O50.100 <- mlegp(xlhs1000_100, ylhs1000_100[,50])
pred.mlegp.O50.100 <- predict(mlegp.O50.100, xtest, se=T)
valstats(c(pred.mlegp.O50.100$f), c(pred.mlegp.O50.100$s), ytest[,50])
mlegp.O50.300 <- mlegp(xlhs1000_300, ylhs1000_300[,50])
pred.mlegp.O50.300 <- predict(mlegp.O50.300, xtest, se=T)
valstats(c(pred.mlegp.O50.300$f), c(pred.mlegp.O50.300$s), ytest[,50])
pred.sgO50_227 <- predict(sgO50_227, xtest)
pred.sgO50_455 <- predict(sgO50_455, xtest)
plot(ytest[,50], pred.sgO50_227$mean - ytest[,50])
points(ytest[,50], pred.sgO50_455$mean - ytest[,50], col=2)
points(ytest[,50], pred.mlegp.O50.100$fi-ytest[,50], col=3, pch=4)
# For output 50
# mlegp.O50.100  0.00293514   NaN      NaN    0.783 0.9999872 0.999974
# mlegp.O50.300  0.001315041   NaN      NaN    0.342 0.9999974 0.9999948
# sgO50_227      0.007884243 -8.801122 0.004248075    0.966 0.9999095 0.9998126
# sgO50_455      0.003161804 -10.30824 0.001793238    0.993 0.999985 0.9999699
# sgO50_1061     0.002270208 -11.16324 0.001183156    0.989 0.9999924 0.9999845
# sgO50_2051     0.001029277 -11.77373 0.0007522737    0.998 0.9999984 0.9999968


