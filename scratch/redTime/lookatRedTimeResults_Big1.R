x100 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n100_s0228_all_input.csv")[,-1]))
y100 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n100_s0228_all_output.csv")[,-1])))
x1000 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
y1000 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
x1000_2 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0304_all_input.csv")[,-1]))
y1000_2 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0304_all_output.csv")[,-1])))


# Ran SGGP on all outputs, w/ 90 supp pts. Didn't save until 1699 b/c of error
rt.sggp.1699 <- readRDS("./scratch/redTime/redTimeData/out_Big1_SGGP-1699.rds")
rt.sggp.2195 <- readRDS("./scratch/redTime/redTimeData/out_Big1_SGGP-2195.rds")
rt.sggp.2695 <- readRDS("./scratch/redTime/redTimeData/out_Big1_SGGP-2695.rds")


# Get stats
stats.rt.sggp.1699 <- CGGPvalstats(rt.sggp.1699, x1000, y1000, bydim=T)
stats.rt.sggp.2195 <- CGGPvalstats(rt.sggp.2195, x1000, y1000, bydim=T)
stats.rt.sggp.2695 <- CGGPvalstats(rt.sggp.2695, x1000, y1000, bydim=T)

# Check stats on 50th dim
CGGPvalstats(CGGPfit(rt.sggp.1699, rt.sggp.1699$Y[,50], Xs=rt.sggp.1699$Xs, Ys=rt.sggp.1699$Ys[,50]), x1000, y1000[,50], bydim=F)



# Run with mlegp
mod.mlegp.50 <- mlegp::mlegp(x100[1:50,], y100[1:50,])
pred.mlegp.50 <- lapply(1:100, function(i) predict(mod.mlegp.50[[i]], x1000, se=T)) %>% {list(fit={do.call(cbind, lapply(., function(i) i$fit))}, se.fit={do.call(cbind, lapply(., function(i) i$se.fit))})}
stats.mlegp.50 <- valstats(pred.mlegp.50$fit, pred.mlegp.50$se, y1000, bydim=F)
mod.mlegp.75 <- mlegp::mlegp(x100[1:75,], y100[1:75,])
pred.mlegp.75 <- lapply(1:100, function(i) predict(mod.mlegp.75[[i]], x1000, se=T)) %>% {list(fit={do.call(cbind, lapply(., function(i) i$fit))}, se.fit={do.call(cbind, lapply(., function(i) i$se.fit))})}
stats.mlegp.75 <- valstats(pred.mlegp.75$fit, pred.mlegp.75$se, y1000, bydim=F)
mod.mlegp.100 <- mlegp::mlegp(x100, y100)
pred.mlegp.100 <- lapply(1:100, function(i) predict(mod.mlegp.100[[i]], x1000, se=T)) %>% {list(fit={do.call(cbind, lapply(., function(i) i$fit))}, se.fit={do.call(cbind, lapply(., function(i) i$se.fit))})}
stats.mlegp.100 <- valstats(pred.mlegp.100$fit, pred.mlegp.100$se, y1000, bydim=F)
mod.mlegp.200 <- mlegp::mlegp(x1000_2, y1000_2)
pred.mlegp.200 <- predict(mod.mlegp.200, x1000, se=T)
stats.mlegp.200 <- valstats(pred.mlegp.200$fit, pred.mlegp.200$se, y1000, bydim=F)
mod.mlegp.300 <- mlegp::mlegp(x1000_2, y1000_2)
pred.mlegp.300 <- predict(mod.mlegp.300, x1000, se=T)
stats.mlegp.300 <- valstats(pred.mlegp.300$fit, pred.mlegp.300$se, y1000, bydim=F)
mod.mlegp.400 <- mlegp::mlegp(x1000_2, y1000_2)
pred.mlegp.400 <- predict(mod.mlegp.400, x1000, se=T)
stats.mlegp.400 <- valstats(pred.mlegp.400$fit, pred.mlegp.400$se, y1000, bydim=F)
mod.mlegp.500 <- mlegp::mlegp(x1000_2, y1000_2)
pred.mlegp.500 <- predict(mod.mlegp.500, x1000, se=T)
stats.mlegp.500 <- valstats(pred.mlegp.500$fit, pred.mlegp.500$se, y1000, bydim=F)




allstats <- list(
  # CGGP
  data.frame("CGGP", 100, 1699, 0.03876258, -5.951558, 0.01723906,  0.92504, 0.9998292, 0.9996583, 0.02941151),
  data.frame("CGGP", 100, 2195, 0.03973626, -5.939358, 0.01805724,  0.94663, 0.9998205, 0.9996409, 0.0301442),
  data.frame("CGGP", 100, 2695, 0.03895927, -5.709819, 0.01752450,  0.89132, 0.9998274, 0.9996548, 0.02955428),
  # mlegp
  data.frame("mlegp", 0, 50, 0.2121086, -2.136166, 0.1122491,   0.9941, 0.9949418, 0.9897676, 0.1614444),
  data.frame("mlegp", 0, 75, 0.1266577, -2.708455, 0.0773198,  0.99936, 0.9982067, 0.9963514, 0.09525388),
  data.frame("mlegp", 0,100, 0.1005995, -2.896368, 0.06750803,   0.9998, 0.9988543, 0.9976983, 0.07583394)
)
allstats <- lapply(allstats, function(x){colnames(x) <- c("Package", 'Nsup',"Ngrid","RMSE","score","CRPscore","coverage","corr","R2","RMSEnorm");x})
allstats <- do.call(rbind, allstats)
allstats$Ntotal <- allstats$Nsup + allstats$Ngrid
library(ggplot2)
ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=interaction(Package), shape=as.factor(Nsup))) + geom_point(size=3)
ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=interaction(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10() + scale_y_log10()
ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=Nsup)) + geom_point(size=3) + facet_grid(. ~ Package) + scale_x_log10() + scale_y_log10()
ggplot(data=allstats %>% filter(score<1e5), mapping=aes(Ntotal, score, color=interaction(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10()




tdf <- rbind(cbind(Ngrid=1699, stats.rt.sggp.1699),
      cbind(Ngrid=2195, stats.rt.sggp.2195),
      cbind(Ngrid=2695, stats.rt.sggp.2695))
ggplot(data=tdf, mapping=aes(Ngrid, RMSE)) + geom_point()
