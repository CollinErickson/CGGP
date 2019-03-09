# Running an experiment to test redTime with all 100 outputs.
# we have 2x2 options.
# Using SGGP data from run that only fit to dimension 50.

# Exp2 is now running with sup data

evfunc <- function(Ngrid, use_PCA, separateoutputparameterdimensions, Nsup) {
  require("SGGP")
  
  if (version$os=="linux-gnu") {
    x1000 <- unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
    y1000 <- log(unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
    x1000_2 <- unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0304_all_input.csv")[,-1]))
    y1000_2 <- log(unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0304_all_output.csv")[,-1])))
    sg.base <- readRDS(paste0("~/scratch/redTime_v0.1/SGGPruns/redTimeTestSup2o50/output_files/out_S2o50_SGGP-", Ngrid, ".rds"))
    Y.8099 <- unname(log(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/redTimeTestSup2o50_all_SGGP_output-8099.csv")[,-1])))
  } else if (version$os == "mingw32") {
    x1000 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
    y1000 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
    x1000_2 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0304_all_input.csv")[,-1]))
    y1000_2 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0304_all_output.csv")[,-1])))
    sg.base <- readRDS(paste0("./scratch/redTime/redTimeData/out_S2o50_SGGP-", Ngrid, ".rds"))
    Y.8099 <- unname(log(as.matrix(read.csv("./scratch/redTime/redTimeData/redTimeTestSup2o50_all_SGGP_output-8099.csv")[,-1])))
  } else {
    stop("version$os doesn't match any")
  }
  
  if (Nsup > 0) {
    supinds <- sample(1:1000, Nsup, replace=F)
    Xs <- x1000_2[supinds,]
    Ys <- y1000_2[supinds,]
  } else {
    Xs <- NULL
    Ys <- NULL
  }
  
  sg <- SGGPfit(sg.base, Y=Y.8099[1:Ngrid,], use_PCA = use_PCA,
                Xs=Xs, Ys=Ys, HandlingSuppData = "Correct",
                separateoutputparameterdimensions = separateoutputparameterdimensions)
  
  SGGPvalstats(sg, x1000, y1000)
}


e1 <- comparer::ffexp$new(
  Ngrid = c(199, 299, 399, 499, 699, 1299, 1299, 2499, 4099, 6099, 8099),
  Nsup = c(0, 100, 200, 300),
  use_PCA = c(T,F),
  separateoutputparameterdimensions = c(T,F),
  eval_func = evfunc,
  folder_path = if (version$os=="linux-gnu") {"/home/collin/scratch/SGGP/scratch/redTime/redTimeMVExp2/Exp2"}
  else if (version$os=="mingw32") {"./scratch/redTime/redTimeMVExp2/Exp2"}
  else {stop("bad folderpath")},
  # folder_path = "./scratch/InternalComparison/redTimeMVout_from_S2o50",
  parallel=T,
  parallel_cores=30
)
e1$recover_parallel_temp_save(delete_after = F)

e1$run_all(run_order="random", delete_parallel_temp_save_after = F, parallel_temp_save = T, write_start_files = T, write_error_files = T)

e1$save_self()

if (F) {
  e1 <- readRDS("C:/Users/cbe117/Documents/GitHub/SGGP/scratch/redTime/redTimeMVExp2/Exp2_167_of_176.rds")
  edf <- e1$outcleandf
  edf$outdim <- rep(1:100, 176)
  edf <- edf[!is.na(edf$RMSE),]
  edf$use_PCAchar <- c("noPCA","PCA")[edf$use_PCA+1]
  edf$sopdchar <- c("1opd", "sopd")[edf$separateoutputparameterdimensions+1]
  edf$N <- edf$Ngrid + edf$Nsup
  
  ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=outdim,shape=interaction(use_PCAchar, sopdchar))) + geom_point() + facet_grid(. ~ as.factor(Nsup)) + scale_y_log10()
  ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=outdim)) + geom_point(size=2) + facet_grid(Nsup ~ use_PCAchar + sopdchar) + scale_y_log10() + scale_x_log10()
  ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar) + scale_x_log10() + scale_y_log10()
  ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=interaction(use_PCA, separateoutputparameterdimensions))) + geom_point()
  
  ggplot(data=edf, mapping=aes(x=N, y=runtime, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar) + scale_x_log10() + scale_y_log10()
  
  ggplot(data=edf, mapping=aes(x=N, y=score, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar) + scale_x_log10()
  
}

# Compare to mlegp

if (F) {
  mlegp.exp <- readRDS("./scratch/redTime/redTimeMVExp1/Exp1_mlegp/Exp1_mlegp_completed.rds")
  mdf <- mlegp.exp$outcleandf
  mdf$outdim <- rep(1:100, 8)
  mdf$use_PCAchar <- c("noPCA","PCA")[mdf$use_PCA+1]
  mdf$sopdchar <- c("1opd", "sopd")[mdf$separateoutputparameterdimensions+1]
  mdf$Ngrid <- mdf$N
  mdf$Nsup <- 0
  
  # Combine the two
  adf <- rbind(cbind(edf,package="SGGP"), cbind(mdf, package="mlegp"))
  
  library(ggplot2)
  ggplot(data=adf, mapping=aes(x=N, y=RMSE, color=outdim,shape=interaction(package, use_PCAchar, sopdchar))) + geom_point()
  ggplot(data=adf, mapping=aes(x=N, y=RMSE, color=outdim)) + geom_point(size=2) + facet_grid(package ~ use_PCAchar + sopdchar)
  ggplot(data=adf, mapping=aes(x=N, y=RMSE, color=outdim)) + geom_point(size=2) + facet_grid(package ~ use_PCAchar + sopdchar) + scale_x_log10() + scale_y_log10()
  ggplot(data=adf, mapping=aes(x=N, y=CRPscore, color=outdim)) + geom_point(size=2) + facet_grid(package ~ use_PCAchar + sopdchar) + scale_x_log10() + scale_y_log10()
  ggplot(data=adf, mapping=aes(x=N, y=RMSE, color=interaction(use_PCA, separateoutputparameterdimensions, package))) + geom_point()
  
  ggplot(data=adf, mapping=aes(x=N, y=runtime, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar) + scale_x_log10() + scale_y_log10()
  
  ggplot(data=adf, mapping=aes(x=N, y=score, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar) + scale_x_log10()
  
  # Too many variables, cut down 100 outdim to single
  adf2 <- plyr::ddply(adf, c("package", "N", "Ngrid", "Nsup", "use_PCAchar", "sopdchar"),
                      function(x) {colMeans(x[,c("RMSE","score","CRPscore","coverage","corr","R2","runtime")])})
  ggplot(data=adf2, mapping=aes(x=N, y=RMSE, color=as.factor(Nsup),shape=interaction(package, use_PCAchar, sopdchar))) + geom_point(size=3)
  ggplot(data=adf2, mapping=aes(x=N, y=RMSE, color=as.factor(Nsup))) + geom_point(size=3) + facet_grid(package ~ use_PCAchar + sopdchar) + scale_x_log10() + scale_y_log10()
  ggplot(data=adf2, mapping=aes(x=N, y=CRPscore, color=as.factor(Nsup))) + geom_point(size=3) + facet_grid(package ~ use_PCAchar + sopdchar) + scale_x_log10() + scale_y_log10()
  
}