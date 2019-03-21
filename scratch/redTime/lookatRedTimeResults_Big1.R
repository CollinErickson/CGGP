# This has the results from Big1 and Big2.
# Big1 used UCB and gave bad results.
# Big2 used Greedy and seemed to work well.
# We fixed UCB while it was running, so UCB should work from now on.

x100 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n100_s0228_all_input.csv")[,-1]))
y100 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n100_s0228_all_output.csv")[,-1])))
x1000 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
y1000 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
x1000_2 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0304_all_input.csv")[,-1]))
y1000_2 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0304_all_output.csv")[,-1])))


# Big1
# Ran SGGP on all outputs, w/ 90 supp pts. Didn't save until 1699 b/c of error
rt.sggp.1699 <- readRDS("./scratch/redTime/redTimeData/out_Big1_SGGP-1699.rds")
rt.sggp.2195 <- readRDS("./scratch/redTime/redTimeData/out_Big1_SGGP-2195.rds")
rt.sggp.2695 <- readRDS("./scratch/redTime/redTimeData/out_Big1_SGGP-2695.rds")
# Big2
r2.sggp.399 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-399.rds")
r2.sggp.599 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-599.rds")
r2.sggp.799 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-799.rds")
r2.sggp.999 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-999.rds")
r2.sggp.1199 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-1199.rds")
r2.sggp.1699 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-1699.rds")
r2.sggp.2199 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-2199.rds")
r2.sggp.2699 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-2699.rds")
r2.sggp.3199 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-3199.rds")
r2.sggp.4199 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-4199.rds")
r2.sggp.5199 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-5199.rds")
r2.sggp.6195 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-6195.rds")
r2.sggp.7187 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-7187.rds")
r2.sggp.8179 <- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-8179.rds")
r2.sggp.10163<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-10163.rds")
r2.sggp.12159<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-12159.rds")
r2.sggp.16159<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-16159.rds")
r2.sggp.18159<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-18159.rds")
r2.sggp.20155<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-20155.rds")
r2.sggp.22155<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-22155.rds")
r2.sggp.24155<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-24155.rds")
r2.sggp.26155<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-26155.rds")
# r2.sggp.<- readRDS("./scratch/redTime/redTimeData/out_Big2_SGGP-.rds")
r2.sggp.26155.Ignore<- CGGPfit(r2.sggp.26155, Y=r2.sggp.26155$Y,Xs=r2.sggp.26155$Xs,Ys=r2.sggp.26155$Ys,HandlingSuppData="Ignore")
r2.sggp.26155.Mat32<- CGGPfit(r2.sggp.26155, Y=r2.sggp.26155$Y,Xs=r2.sggp.26155$Xs,Ys=r2.sggp.26155$Ys,corr="m3")
r2.sggp.26155.18ktheta <- CGGPfit(r2.sggp.26155, Y=r2.sggp.26155$Y,Xs=r2.sggp.26155$Xs,Ys=r2.sggp.26155$Ys,set_thetaMAP_to = r2.sggp.18159$thetaMAP)
r2.sggp.26155.18kthetanosupp <- CGGPfit(r2.sggp.26155, Y=r2.sggp.26155$Y,set_thetaMAP_to = r2.sggp.18159$thetaMAP)

# Get stats
# Big1
stats.rt.sggp.1699 <- CGGPvalstats(rt.sggp.1699, x1000, y1000, bydim=T)
stats.rt.sggp.2195 <- CGGPvalstats(rt.sggp.2195, x1000, y1000, bydim=T)
stats.rt.sggp.2695 <- CGGPvalstats(rt.sggp.2695, x1000, y1000, bydim=T)
# Big2
stats.r2.sggp.399 <- CGGPvalstats(r2.sggp.399, x1000, y1000, bydim=F)
stats.r2.sggp.599 <- CGGPvalstats(r2.sggp.599, x1000, y1000, bydim=F)
stats.r2.sggp.799 <- CGGPvalstats(r2.sggp.799, x1000, y1000, bydim=F)
stats.r2.sggp.999 <- CGGPvalstats(r2.sggp.999, x1000, y1000, bydim=F)
stats.r2.sggp.1199 <- CGGPvalstats(r2.sggp.1199, x1000, y1000, bydim=F)
stats.r2.sggp.1699 <- CGGPvalstats(r2.sggp.1699, x1000, y1000, bydim=F)
stats.r2.sggp.2199 <- CGGPvalstats(r2.sggp.2199, x1000, y1000, bydim=F)
stats.r2.sggp.2699 <- CGGPvalstats(r2.sggp.2699, x1000, y1000, bydim=F)
stats.r2.sggp.3199 <- CGGPvalstats(r2.sggp.3199, x1000, y1000, bydim=F)
stats.r2.sggp.4199 <- CGGPvalstats(r2.sggp.4199, x1000, y1000, bydim=F)
stats.r2.sggp.5199 <- CGGPvalstats(r2.sggp.5199, x1000, y1000, bydim=F)
stats.r2.sggp.6195 <- CGGPvalstats(r2.sggp.6195, x1000, y1000, bydim=F)
stats.r2.sggp.7187 <- CGGPvalstats(r2.sggp.7187, x1000, y1000, bydim=F)
stats.r2.sggp.8179 <- CGGPvalstats(r2.sggp.8179, x1000, y1000, bydim=F)
stats.r2.sggp.10163<- CGGPvalstats(r2.sggp.10163, x1000, y1000, bydim=F)
stats.r2.sggp.12159<- CGGPvalstats(r2.sggp.12159, x1000, y1000, bydim=F)
stats.r2.sggp.16159<- CGGPvalstats(r2.sggp.18159, x1000, y1000, bydim=F)
stats.r2.sggp.18159<- CGGPvalstats(r2.sggp.18159, x1000, y1000, bydim=F)
stats.r2.sggp.20155<- CGGPvalstats(r2.sggp.20155, x1000, y1000, bydim=F)
stats.r2.sggp.22155<- CGGPvalstats(r2.sggp.22155, x1000, y1000, bydim=F)
stats.r2.sggp.24155<- CGGPvalstats(r2.sggp.24155, x1000, y1000, bydim=F)
stats.r2.sggp.26155<- CGGPvalstats(r2.sggp.26155, x1000, y1000, bydim=F)
# stats.r2.sggp.<- CGGPvalstats(r2.sggp., x1000, y1000, bydim=F)
stats.r2.sggp.26155.Ignore<- CGGPvalstats(r2.sggp.26155.Ignore, x1000, y1000, bydim=F)
stats.r2.sggp.26155.Mat32<- CGGPvalstats(r2.sggp.26155.Mat32, x1000, y1000, bydim=F)
stats.r2.sggp.26155.18ktheta<- CGGPvalstats(r2.sggp.26155.18ktheta, x1000, y1000, bydim=F)
stats.r2.sggp.26155.18kthetanosupp<- CGGPvalstats(r2.sggp.26155.18kthetanosupp, x1000, y1000, bydim=F)


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
mod.mlegp.200 <- sample(1:1000, 200) %>% {mlegp::mlegp(x1000_2[.,], y1000_2[.,])}
pred.mlegp.200 <- lapply(1:100, function(i) predict(mod.mlegp.200[[i]], x1000, se=T)) %>% {list(fit={do.call(cbind, lapply(., function(i) i$fit))}, se.fit={do.call(cbind, lapply(., function(i) i$se.fit))})}
stats.mlegp.200 <- valstats(pred.mlegp.200$fit, pred.mlegp.200$se^2, y1000, bydim=F)
# mod.mlegp.300 <- mlegp::mlegp(x1000_2, y1000_2)
# pred.mlegp.300 <- predict(mod.mlegp.300, x1000, se=T)
# stats.mlegp.300 <- valstats(pred.mlegp.300$fit, pred.mlegp.300$se^2, y1000, bydim=F)
# mod.mlegp.400 <- mlegp::mlegp(x1000_2, y1000_2)
# pred.mlegp.400 <- predict(mod.mlegp.400, x1000, se=T)
# stats.mlegp.400 <- valstats(pred.mlegp.400$fit, pred.mlegp.400$se^2, y1000, bydim=F)
# mod.mlegp.500 <- mlegp::mlegp(x1000_2, y1000_2)
# pred.mlegp.500 <- predict(mod.mlegp.500, x1000, se=T)
# stats.mlegp.500 <- valstats(pred.mlegp.500$fit, pred.mlegp.500$se^2, y1000, bydim=F)




allstats <- list(
  # CGGP
  data.frame("CGGP1", 90, 1699, 0.03876258, -5.951558, 0.01723906,  0.92504, 0.9998292, 0.9996583, 0.02941151),
  data.frame("CGGP1", 90, 2195, 0.03973626, -5.939358, 0.01805724,  0.94663, 0.9998205, 0.9996409, 0.0301442),
  data.frame("CGGP1", 90, 2695, 0.03895927, -5.709819, 0.01752450,  0.89132, 0.9998274, 0.9996548, 0.02955428),
  data.frame("CGGP", 90, 399,  0.03292523, -6.292347, 0.01476177,   0.93929, 0.9998769, 0.9997534, 0.02502954),
  data.frame("CGGP", 90, 599,  0.02445636, -6.717207, 0.01163083,   0.95573, 0.999932,  0.999864,  0.01859961),
  data.frame("CGGP", 90, 799,  0.02426409, -6.831546, 0.01119061,   0.95034, 0.9999334, 0.9998661, 0.01841507),
  data.frame("CGGP", 90, 999,  0.02418435, -6.846158, 0.01062984,   0.95242, 0.9999338, 0.999867,  0.01836773),
  data.frame("CGGP", 90, 1199, 0.02205223, -7.085375, 0.009788555,  0.95352, 0.9999454, 0.9998894, 0.01662363),
  data.frame("CGGP", 90, 1699, 0.01527591, -7.696946, 0.007269341,  0.97937, 0.9999735, 0.9999469, 0.01145086),
  data.frame("CGGP", 90, 2199, 0.01370636, -7.730417, 0.006900497,   0.9812, 0.9999787, 0.9999573, 0.01036107),
  data.frame("CGGP", 90, 2699, 0.014988  , -7.711674, 0.006816212,  0.98524, 0.9999745, 0.9999489, 0.01092403),
  data.frame("CGGP", 90, 3199, 0.01357989, -8.056265, 0.005935432,  0.98193, 0.9999791, 0.9999581, 0.009891518),
  data.frame("CGGP", 90, 4199, 0.01130804, -8.147965, 0.005525203,  0.97767, 0.9999855, 0.9999709, 0.008477671),
  data.frame("CGGP", 90, 5199, 0.009739982,-8.454266, 0.004678459,  0.98001, 0.9999892, 0.9999784, 0.007307837),
  data.frame("CGGP", 90, 6195, 0.007799561, -8.601915, 0.004069435,  0.99209, 0.9999931, 0.9999862, 0.005720459),
  data.frame("CGGP", 90, 7187, 0.007602274, -8.832714, 0.003677643,   0.9891, 0.9999934, 0.9999869, 0.005552559),
  data.frame("CGGP", 90, 8179, 0.006772034, -9.045547, 0.003310222,  0.98978, 0.9999948, 0.9999896, 0.004955794),
  data.frame("CGGP", 90, 10163,0.005923366, -9.315209, 0.002928836,  0.98941, 0.999996,  0.999992,  0.004356485),
  data.frame("CGGP", 90, 12159,0.004279414, -9.339384, 0.002790037,  0.99488, 0.9999979, 0.9999958, 0.003197408),
  data.frame("CGGP", 90, 16159,0.003443654, -9.731616, 0.002267292,  0.99656, 0.9999987, 0.9999973, 0.002545515),
  data.frame("CGGP", 90, 18159,0.003443654, -9.731616, 0.002267292,  0.99656, 0.9999987, 0.9999973, 0.002545515),
  data.frame("CGGP", 90, 20155,0.003984703, -9.738379, 0.002263422,  0.99591, 0.9999982, 0.9999964, 0.002905657),
  data.frame("CGGP", 90, 22155,0.004459119, -9.739725, 0.002274312,  0.99623, 0.9999977, 0.9999955, 0.00319555),
  data.frame("CGGP", 90, 24155,0.004439284, -9.703938, 0.002317367,  0.99618, 0.9999978, 0.9999955, 0.003193208),
  data.frame("CGGP", 90, 26155,0.004961369, -9.659113, 0.002383479,  0.99554, 0.9999972, 0.9999944, 0.00355455),
  data.frame("CGGPIgnore", 90, 26155,0.004959593, -9.656147, 0.002386668,  0.99554, 0.9999972, 0.9999944, 0.003553383),
  data.frame("CGGPMat32", 90, 26155, 0.004838821, -8.49521 , 0.003830926,  0.99942, 0.9999973, 0.9999947, 0.003493466),
  data.frame("CGGP18ktheta", 90, 26155, 0.004768937, -7.572385, 0.005782846,  0.99923, 0.9999974, 0.9999948, 0.003410945),
  data.frame("CGGP18kthetanosupp", 90, 26155, 0.00473138, -7.165235, 0.007330487,  0.99932, 0.9999975, 0.9999949, 0.003383729),
  # data.frame("CGGP", 90, 2199, ),
  # data.frame("CGGP", 90, , ),
  # data.frame("CGGP", 90, , ),
  # mlegp
  data.frame("mlegp", 0, 50, 0.2121086, -2.136166, 0.1122491,   0.9941, 0.9949418, 0.9897676, 0.1614444),
  data.frame("mlegp", 0, 75, 0.1266577, -2.708455, 0.0773198,  0.99936, 0.9982067, 0.9963514, 0.09525388),
  data.frame("mlegp", 0,100, 0.1005995, -2.896368, 0.06750803,   0.9998, 0.9988543, 0.9976983, 0.07583394)
)
allstats <- lapply(allstats, function(x){colnames(x) <- c("Package", 'Nsup',"Ngrid","RMSE","score","CRPscore","coverage","corr","R2","RMSEnorm");x})
allstats <- do.call(rbind, allstats)
allstats$Ntotal <- allstats$Nsup + allstats$Ngrid
library(ggplot2)
ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_y_log10()
ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10() + scale_y_log10()
ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=Nsup)) + geom_point(size=3) + facet_grid(. ~ Package) + scale_x_log10() + scale_y_log10()
ggplot(data=allstats, mapping=aes(Ntotal, score, color=(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10()
ggplot(data=allstats, mapping=aes(Ntotal, CRPscore, color=(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10() + scale_y_log10()
# ggplot(data=allstats %>% filter(score<1e5), mapping=aes(Ntotal, score, color=(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10()




tdf <- rbind(cbind(Ngrid=1699, stats.rt.sggp.1699),
      cbind(Ngrid=2195, stats.rt.sggp.2195),
      cbind(Ngrid=2695, stats.rt.sggp.2695))
ggplot(data=tdf, mapping=aes(Ngrid, RMSE)) + geom_point()
