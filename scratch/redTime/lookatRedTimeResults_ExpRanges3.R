# This has results from ExpandedRanges3.
# Done after thesis, to be used in publication.
# ER3a used no supp, Greedy, power exp.
# Only using first 80 output dimensions.

x100 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges3_LHS1L_n100_s0228_all_input.csv")[,-1]))
y100 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges3_LHS1L_n100_s0228_all_output.csv")[,-1])))[,1:80]
x1000 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges3_LHS1L_n1000_s0429_all_input.csv")[,-1]))
y1000 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges3_LHS1L_n1000_s0429_all_output.csv")[,-1])))[,1:80]
x1000_2 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges3_LHS1L_n1000_s0304_all_input.csv")[,-1]))
y1000_2 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges3_LHS1L_n1000_s0304_all_output.csv")[,-1])))[,1:80]

# ER3a
ra.sggp.399 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-399.rds")
ra.sggp.599 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-599.rds")
ra.sggp.799 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-799.rds")
ra.sggp.999 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-999.rds")
ra.sggp.1199 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-1199.rds")
ra.sggp.1699 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-1699.rds")
ra.sggp.2199 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-2199.rds")
ra.sggp.2699 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-2699.rds")
ra.sggp.3195 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-3195.rds")
ra.sggp.4191 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-4191.rds")
ra.sggp.5187 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-5187.rds")
ra.sggp.6179 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-6179.rds")
ra.sggp.7171 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-7171.rds")
ra.sggp.8171 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-8171.rds")
ra.sggp.10159 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-10159.rds")
ra.sggp.12155 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-12155.rds")
ra.sggp.16155 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-16155.rds")
ra.sggp.18155 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-18155.rds")
# ra.sggp.20159 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-20159.rds")
ra.sggp.30151 <- readRDS("./scratch/redTime/redTimeData/out_ER3a_SGGP-30151.rds")
# ra.sggp.20159.PE<- CGGPfit(ra.sggp.20159, Y=ra.sggp.20159$Y,Xs=ra.sggp.20159$Xs,Ys=ra.sggp.20159$Ys,corr="PowerExp")
ra.sggp.30151.PE<- CGGPfit(ra.sggp.30151, Y=ra.sggp.30151$Y,Xs=ra.sggp.30151$Xs,Ys=ra.sggp.30151$Ys,corr="PowerExp")
# ra.sggp.20159.C <- CGGPfit(ra.sggp.20159, Y=ra.sggp.20159$Y,Xs=ra.sggp.20159$Xs,Ys=ra.sggp.20159$Ys,corr="Cauchy")
ra.sggp.30151.C <- CGGPfit(ra.sggp.30151, Y=ra.sggp.30151$Y,Xs=ra.sggp.30151$Xs,Ys=ra.sggp.30151$Ys,corr="Cauchy")



# Get stats
# Big 3
stats.ra.sggp.399  <- CGGPvalstats(ra.sggp.399,   x1000, y1000, bydim=F)
stats.ra.sggp.599  <- CGGPvalstats(ra.sggp.599,   x1000, y1000, bydim=F)
stats.ra.sggp.799  <- CGGPvalstats(ra.sggp.799,   x1000, y1000, bydim=F)
stats.ra.sggp.999  <- CGGPvalstats(ra.sggp.999,   x1000, y1000, bydim=F)
stats.ra.sggp.1199 <- CGGPvalstats(ra.sggp.1199,  x1000, y1000, bydim=F)
stats.ra.sggp.1699 <- CGGPvalstats(ra.sggp.1699,  x1000, y1000, bydim=F)
stats.ra.sggp.2199 <- CGGPvalstats(ra.sggp.2199,  x1000, y1000, bydim=F)
stats.ra.sggp.2699 <- CGGPvalstats(ra.sggp.2699,  x1000, y1000, bydim=F)
stats.ra.sggp.3195 <- CGGPvalstats(ra.sggp.3195,  x1000, y1000, bydim=F)
stats.ra.sggp.4191 <- CGGPvalstats(ra.sggp.4191,  x1000, y1000, bydim=F)
stats.ra.sggp.5187 <- CGGPvalstats(ra.sggp.5187,  x1000, y1000, bydim=F)
stats.ra.sggp.6179 <- CGGPvalstats(ra.sggp.6179,  x1000, y1000, bydim=F)
stats.ra.sggp.7171 <- CGGPvalstats(ra.sggp.7171,  x1000, y1000, bydim=F)
stats.ra.sggp.8171 <- CGGPvalstats(ra.sggp.8171,  x1000, y1000, bydim=F)
stats.ra.sggp.10159<- CGGPvalstats(ra.sggp.10159, x1000, y1000, bydim=F)
stats.ra.sggp.12155<- CGGPvalstats(ra.sggp.12155, x1000, y1000, bydim=F)
stats.ra.sggp.16155<- CGGPvalstats(ra.sggp.16155, x1000, y1000, bydim=F)
stats.ra.sggp.18155<- CGGPvalstats(ra.sggp.18155, x1000, y1000, bydim=F)
# stats.ra.sggp.20159<- CGGPvalstats(ra.sggp.20159, x1000, y1000, bydim=F)
stats.ra.sggp.30151<- CGGPvalstats(ra.sggp.30151, x1000, y1000, bydim=F)
# stats.ra.sggp.20159.PE<- CGGPvalstats(ra.sggp.20159.PE, x1000, y1000, bydim=F)
stats.ra.sggp.30151.PE<- CGGPvalstats(ra.sggp.30151.PE, x1000, y1000, bydim=F)
# stats.ra.sggp.20159.C <- CGGPvalstats(ra.sggp.20159.C , x1000, y1000, bydim=F)
stats.ra.sggp.30151.C <- CGGPvalstats(ra.sggp.30151.C , x1000, y1000, bydim=F)

allstats <- list(
  # CGGP
  data.frame("CGGP", 0,   399, 0.1270009 , -3.61207 , 0.05555768,  0.94545 , 0.9969352, 0.9935222, 0.09223503),
  data.frame("CGGP", 0,   599, 0.08656919, -4.095473, 0.03655804, 0.954425 , 0.9985132, 0.9969902, 0.06312542),
  data.frame("CGGP", 0,   799, 0.07907009, -4.313236, 0.03075301, 0.9594125, 0.9987519, 0.9974891, 0.05773612),
  data.frame("CGGP", 0,   999, 0.0750431 , -4.276698, 0.02822047,   0.9563 , 0.9988783, 0.9977383, 0.05491472),
  data.frame("CGGP", 0,  1199, 0.06771744, -4.496644, 0.0243313 ,  0.96335 , 0.9990838, 0.9981583, 0.04956006),
  data.frame("CGGP", 0,  1699, 0.04583473, -5.088112, 0.0241948 , 0.9818875, 0.9995893, 0.9991563, 0.03327931),
  data.frame("CGGP", 0,  1199, 0.03068066, -4.947636, 0.02252122, 0.9933375, 0.9998112, 0.999622,  0.02247564),
  data.frame("CGGP", 0,  2699, 0.02704855, -5.738293, 0.01620703, 0.9936125, 0.9998545, 0.9997062, 0.01972948),
  data.frame("CGGP", 0,  3195, 0.02155968, -6.241509, 0.0126585,  0.9906125, 0.9999072, 0.9998133, 0.01573683),
  data.frame("CGGP", 0,  4191, 0.02113381, -6.562632, 0.01087367, 0.985575 , 0.9999103, 0.9998206, 0.01543448),
  data.frame("CGGP", 0,  5187, 0.01550605, -6.789563, 0.00939717, 0.993375 , 0.9999519, 0.9999034, 0.01128335),
  data.frame("CGGP", 0,  6179, 0.01353115, -6.940049, 0.008560179, 0.992275, 0.9999634, 0.9999265, 0.009878027),
  data.frame("CGGP", 0,  7171, 0.01272007, -7.123897, 0.007880539,0.9931625, 0.9999676, 0.999935,  0.009295405),
  data.frame("CGGP", 0,  8171, 0.01134226, -7.209546, 0.007439439,0.9928125, 0.9999742, 0.9999483, 0.008268045),
  data.frame("CGGP", 0, 10159, 0.01038528, -7.561542, 0.00635735, 0.9885125, 0.9999784, 0.9999567, 0.0075761),
  data.frame("CGGP", 0, 12159, 0.009005683,-7.803026, 0.005657783, 0.991825, 0.9999838, 0.9999674, 0.006588723),
  data.frame("CGGP", 0, 16155, 0.007374222,-7.991828, 0.005007521,0.9932375, 0.9999891, 0.9999782, 0.005393246),
  data.frame("CGGP", 0, 18155, 0.006947166,-8.225993, 0.004492074, 0.991475, 0.9999903, 0.9999806, 0.005082594),
  data.frame("CGGP", 0, 30151, 0.005188444,-8.647838, 0.003515817,   0.9932, 0.9999946, 0.9999892, 0.003800047)
)

allstats <- lapply(allstats, function(x){colnames(x) <- c("Package", 'Nsup',"Ngrid","RMSE","score","CRPscore","coverage","corr","R2","RMSEnorm");x})
allstats <- do.call(rbind, allstats)
allstats$Ntotal <- allstats$Nsup + allstats$Ngrid
allstats$Package <- as.character(allstats$Package)
library(ggplot2)
# ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=Package, shape=as.factor(Nsup))) + geom_point(size=3) + scale_y_log10()
ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=Package, shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10() + scale_y_log10()
ggplot(data=allstats, mapping=aes(Ntotal, RMSE, color=Nsup)) + geom_point(size=3) + facet_grid(. ~ Package) + scale_x_log10() + scale_y_log10()
ggplot(data=allstats, mapping=aes(Ntotal, score, color=(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10()
ggplot(data=allstats, mapping=aes(Ntotal, CRPscore, color=(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10() + scale_y_log10()
# ggplot(data=allstats %>% filter(score<1e5), mapping=aes(Ntotal, score, color=(Package), shape=as.factor(Nsup))) + geom_point(size=3) + scale_x_log10()

