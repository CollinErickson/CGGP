# Using expandedRanges2 on redTime

x100 <- unname(as.matrix(read.csv("../../../Desktop/redTimeData/ExpandedRanges2_LHS1L_n100_s0228_all_input.csv")[,-1]))
y100 <- log(unname(as.matrix(read.csv("../../../Desktop/redTimeData/ExpandedRanges2_LHS1L_n100_s0228_all_output.csv")[,-1])))
x1000 <- unname(as.matrix(read.csv("../../../Desktop/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
y1000 <- log(unname(as.matrix(read.csv("../../../Desktop/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))


# Ran SGGP on output dimension 50
rt.sggp.199 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-199.rds")
rt.sggp.299 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-299.rds")
rt.sggp.399 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-399.rds")
rt.sggp.499 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-499.rds")
rt.sggp.599 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-599.rds")
rt.sggp.699 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-699.rds")
rt.sggp.799 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-799.rds")
rt.sggp.899 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-899.rds")
rt.sggp.999 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-999.rds")
rt.sggp.1099 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-1099.rds")
rt.sggp.1299 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-1299.rds")
rt.sggp.1699 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-1699.rds")
rt.sggp.2099 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-2099.rds")
rt.sggp.2499 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-2499.rds")
rt.sggp.3099 <- readRDS("../../../Desktop/redTimeData/out_S2o50_SGGP-3099.rds")
stats.rt.sggp.199 <- SGGPvalstats(rt.sggp.199, x1000, y1000[,50])
stats.rt.sggp.299 <- SGGPvalstats(rt.sggp.299, x1000, y1000[,50])
stats.rt.sggp.399 <- SGGPvalstats(rt.sggp.399, x1000, y1000[,50])
stats.rt.sggp.499 <- SGGPvalstats(rt.sggp.499, x1000, y1000[,50])
stats.rt.sggp.599 <- SGGPvalstats(rt.sggp.599, x1000, y1000[,50])
stats.rt.sggp.699 <- SGGPvalstats(rt.sggp.699, x1000, y1000[,50])
stats.rt.sggp.799 <- SGGPvalstats(rt.sggp.799, x1000, y1000[,50])
stats.rt.sggp.899 <- SGGPvalstats(rt.sggp.899, x1000, y1000[,50])
stats.rt.sggp.999 <- SGGPvalstats(rt.sggp.999, x1000, y1000[,50])
stats.rt.sggp.1099 <- SGGPvalstats(rt.sggp.1099, x1000, y1000[,50])
stats.rt.sggp.1299 <- SGGPvalstats(rt.sggp.1299, x1000, y1000[,50])
stats.rt.sggp.1699 <- SGGPvalstats(rt.sggp.1699, x1000, y1000[,50])
stats.rt.sggp.2099 <- SGGPvalstats(rt.sggp.2099, x1000, y1000[,50])
stats.rt.sggp.2499 <- SGGPvalstats(rt.sggp.2499, x1000, y1000[,50])
stats.rt.sggp.3099 <- SGGPvalstats(rt.sggp.3099, x1000, y1000[,50])

# SGGP with other correlation functions
stats.rt.sggp.199.m3 <- SGGPvalstats(SGGPfit(rt.sggp.199, rt.sggp.199$Y, Xs=rt.sggp.199$Xs, Ys=rt.sggp.199$Ys, corr="m32"), x1000, y1000[,50])
stats.rt.sggp.199.cauchy <- SGGPvalstats(SGGPfit(rt.sggp.199, rt.sggp.199$Y, Xs=rt.sggp.199$Xs, Ys=rt.sggp.199$Ys, corr="cauchy"), x1000, y1000[,50])
stats.rt.sggp.199.pe <- SGGPvalstats(SGGPfit(rt.sggp.199, rt.sggp.199$Y, Xs=rt.sggp.199$Xs, Ys=rt.sggp.199$Ys, corr="pe"), x1000, y1000[,50])


# Run with mlegp
mod.mlegp.50 <- mlegp::mlegp(x100[1:50,], y100[1:50,50])
pred.mlegp.50 <- predict(mod.mlegp.50, x1000, se=T)
stats.mlegp.50 <- valstats(pred.mlegp.50$fit, pred.mlegp.50$se, y1000[,50])
mod.mlegp.75 <- mlegp::mlegp(x100[1:75,], y100[1:75,50])
pred.mlegp.75 <- predict(mod.mlegp.75, x1000, se=T)
stats.mlegp.75 <- valstats(pred.mlegp.75$fit, pred.mlegp.75$se, y1000[,50])
mod.mlegp.100 <- mlegp::mlegp(x100, y100[,50])
pred.mlegp.100 <- predict(mod.mlegp.100, x1000, se=T)
stats.mlegp.100 <- valstats(pred.mlegp.100$fit, pred.mlegp.100$se, y1000[,50])

# Run with DK
mod.DK.100 <- DiceKriging::km(design=x100, response=y100[,50])
pred.DK.100 <- DiceKriging::predict.km(mod.DK.100, x1000, se=T, type = "SK")
stats.DK.100 <- valstats(pred.DK.100$mean, pred.DK.100$sd, y1000[,50])

# 100+199 0.03076381 -6.262226 0.01522756    0.927 0.9996693 0.9993362
stats50 <- list(
  # SGGP as fit while running redTime
   data.frame("SGGP", 100, 199, 0.03076381, -6.262226, 0.01522756,    0.927, 0.9996693, 0.9993362)
  ,data.frame("SGGP", 100, 299, 0.03451267, -5.969878, 0.01650027,    0.904, 0.9995829, 0.9991646)
  ,data.frame("SGGP", 100, 399, 0.03300229, -6.088676, 0.01463064,    0.913, 0.9996187, 0.9992361)
  ,data.frame("SGGP", 100, 499, 0.02625919, -6.672529, 0.0116401,    0.945, 0.9997592, 0.9995164)
  ,data.frame("SGGP", 100, 599, 0.02180873, -6.890593, 0.01057118 ,    0.951, 0.9998343, 0.9996664)
  ,data.frame("SGGP", 100, 699, 0.01777888, -7.284645, 0.008713988,    0.945, 0.9998899, 0.9997783)
  ,data.frame("SGGP", 100, 799, 0.02116162, -7.101041, 0.009290412,    0.946, 0.999843 , 0.9996859)
  ,data.frame("SGGP", 100, 899, 0.02314205, -6.917019, 0.009970185,    0.959, 0.9998133, 0.9996244)
  ,data.frame("SGGP", 100, 999, 0.02157915, -7.034552, 0.009963331,     0.947, 0.9998376, 0.9996734)
  ,data.frame("SGGP", 100, 1099, 0.01534432, -7.531176, 0.007727928,     0.969, 0.999918 , 0.9998349)
  ,data.frame("SGGP", 100, 1299, 0.01276961, -7.843283, 0.00664756 ,     0.97 , 0.9999435, 0.9998856)
  ,data.frame("SGGP", 100, 1699, 0.01244242, -7.857461, 0.006426012,     0.97 , 0.9999461, 0.9998914)
  ,data.frame("SGGP", 100, 2099, 0.01097685, -8.078817, 0.005756391,    0.964, 0.9999579, 0.9999155)
  ,data.frame("SGGP", 100, 2499, 0.0100807,  -8.237784, 0.005289297,    0.966, 0.9999644, 0.9999287)
  ,data.frame("SGGP", 100, 3099, 0.008719812, -8.298772, 0.004743257,   0.981, 0.9999734, 0.9999467)
  
  
  ,data.frame("SGGP", 100, 199, 0.03216928, -6.083981, 0.01645107,    0.981, 0.9996371, 0.9992742, "m32")
  # ,data.frame("SGGP", 100, 199, )
  # ,data.frame("SGGP", 100, 199, )
  ,data.frame("mlegp", 50, 0, 0.2270544, -1.913619, 0.1250952,    0.999, 0.9831718, 0.9638425)
  ,data.frame("mlegp", 75, 0, 0.10766, -2.792874, 0.07135557,        1, 0.996209, 0.9918708)
  ,data.frame("mlegp", 100, 0, 0.06409675, -3.36158, 0.05061328,        1, 0.9985716, 0.9971186)
  # ,data.frame("mlegp", 0, 100, )
  # ,data.frame("mlegp", 0, 100, )
  ,data.frame("DK", 100, 0, 0.1124202, -2.070725, 0.09362545,        1, 0.9958498, 0.9911361)
)
stats50 <- lapply(stats50, function(x){colnames(x) <- c("Package", 'Nsup',"Ngrid","RMSE","score","CRPscore","coverage","corr","R2");x})
stats50 <- do.call(rbind, stats50)
stats50$Ntotal <- stats50$Nsup + stats50$Ngrid
library(ggplot2)
ggplot(data=stats50, mapping=aes(Ntotal, RMSE, color=Package, shape=as.factor(Nsup))) + geom_point(size=3)
ggplot(data=stats50, mapping=aes(Ntotal, RMSE, color=Package, shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10() + scale_y_log10()
