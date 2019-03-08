# Running an experiment to test redTime with all 100 outputs.
# we have 2x2 options.
# Using SGGP data from run that only fit to dimension 50.

evfunc <- function(N, use_PCA, separateoutputparameterdimensions) {
  require("SGGP")
  
  if (version$os=="linux-gnu") {
    x1000 <- unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
    y1000 <- log(unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
    sg.base <- readRDS(paste0("~/scratch/redTime_v0.1/SGGPruns/redTimeTestSup2o50/output_files/out_S2o50_SGGP-", N, ".rds"))
    Y.8099 <- unname(log(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/redTimeTestSup2o50_all_SGGP_output-8099.csv")[,-1])))
  } else if (version$os == "mingw32") {
    x1000 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
    y1000 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
    sg.base <- readRDS(paste0("./scratch/redTime/redTimeData/out_S2o50_SGGP-", N, ".rds"))
    Y.8099 <- unname(log(as.matrix(read.csv("./scratch/redTime/redTimeData/redTimeTestSup2o50_all_SGGP_output-6889.csv")[,-1])))
  } else {
    stop("version$os doesn't match any")
  }
  
  
  sg <- SGGPfit(sg.base, Y=Y.8099[1:N,], use_PCA = use_PCA,
                separateoutputparameterdimensions = separateoutputparameterdimensions)
  
  SGGPvalstats(sg, x1000, y1000)
}


e1 <- comparer::ffexp$new(
  N = c(199, 299, 399, 499, 699, 1299, 1299, 2499, 4099, 6099, 8099),
  use_PCA = c(T,F),
  separateoutputparameterdimensions = c(T,F),
  eval_func = evfunc,
  folder_path = if (version$os=="linux-gnu") {"/home/collin/scratch/SGGP/scratch/redTime/redTimeMVExp1/Exp1"}
  else if (version$os=="mingw32") {"./scratch/redTime/redTimeMVExp1/Exp1"}
  else {stop("bad folderpath")},
  # folder_path = "./scratch/InternalComparison/redTimeMVout_from_S2o50",
  parallel=T,
  parallel_cores=10
)
e1$recover_parallel_temp_save(delete_after = F)

e1$run_all(delete_parallel_temp_save_after = F, parallel_temp_save = T, write_start_files = T, write_error_files = T)

e1$save_self()

if (F) {
  e1 <- readRDS("C:/Users/cbe117/Documents/GitHub/SGGP/scratch/redTime/redTimeMVExp1/redTimeMVExp1_completed.rds")
  edf <- e1$outcleandf
  edf$outdim <- rep(1:100, 44)
  edf$use_PCAchar <- c("noPCA","PCA")[edf$use_PCA+1]
  edf$sopdchar <- c("1opd", "sopd")[edf$separateoutputparameterdimensions+1]
  
  ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=outdim,shape=interaction(use_PCAchar, sopdchar))) + geom_point()
  ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar)
  ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar) + scale_x_log10() + scale_y_log10()
  ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=interaction(use_PCA, separateoutputparameterdimensions))) + geom_point()
  
  ggplot(data=edf, mapping=aes(x=N, y=runtime, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar) + scale_x_log10() + scale_y_log10()
  
  ggplot(data=edf, mapping=aes(x=N, y=score, color=outdim)) + geom_point(size=2) + facet_grid(use_PCAchar ~ sopdchar) + scale_x_log10()
  
}

if (F) {
  # Run mlegp to combine results
  mlegp_evalfunc <- function(N, use_PCA, separateoutputparameterdimensions) {
    
    
    if (version$os=="linux-gnu") {
      x1000 <- unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
      y1000 <- log(unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
      x1000_2 <- unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0304_all_input.csv")[,-1]))
      y1000_2 <- log(unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0304_all_output.csv")[,-1])))
    } else if (version$os == "mingw32") {
      x1000 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
      y1000 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
      x1000_2 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0304_all_input.csv")[,-1]))
      y1000_2 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0304_all_output.csv")[,-1])))
    } else {
      stop("version$os doesn't match any")
    }
    inds <- sample(1:1000, N, replace=F)
    X <- x1000_2[inds,]
    Y <- y1000_2[inds,]
    if (use_PCA) {
      mod.mlegp.pca <- mlegp::mlegp(X=X, Z=t(Y),
                                        PC.percent = 99.999)
      
      pred.mlegp.pca.pretrans <- lapply(1:mod.mlegp.pca$numGPs,
                                        function(i)predict(mod.mlegp.pca[[i]], x1000, se.fit=T))
      pred.mlegp.pca.pretrans.means <- do.call(cbind, lapply(pred.mlegp.pca.pretrans, function(x) x$fit))
      pred.mlegp.pca.pretrans.ses <- do.call(cbind, lapply(pred.mlegp.pca.pretrans, function(x) x$se))
      pred.mlegp <- list(
        mean = t(mod.mlegp.pca$UD %*% t(pred.mlegp.pca.pretrans.means)),
        var = (t(mod.mlegp.pca$UD %*% t(pred.mlegp.pca.pretrans.ses)))^2)
      
    } else {
      # browser()
      mod.mlegp <- mlegp::mlegp(X=X, Z=Y)
      pred.mlegp.raw <- lapply(1:mod.mlegp$numGPs,
                               function(i)predict(mod.mlegp[[i]], x1000, se.fit=T))
      pred.mlegp <- list(
        mean = do.call(cbind, lapply(pred.mlegp.raw, function(x) x$fit)),
        var = do.call(cbind, lapply(pred.mlegp.raw, function(x) x$se^2)))
    }
    # browser()
    SGGP::valstats(pred.mlegp$mean, pred.mlegp$var, y1000)
  }
  mlegp.exp <- comparer::ffexp$new(
    N=c(100, 200, 300, 500),
    use_PCA=c(T, F),
    separateoutputparameterdimensions=c(T),
    eval_func=function(...)mlegp_evalfunc(...),
    parallel=F,
    folder_path = "./scratch/redTime/redTimeMVExp1/Exp1_mlegp"
  )
  mlegp.exp$rungrid
  # mlegp.exp$run_one(5)
  mlegp.exp$run_all()
  if (F) {
    mlegp.exp <- readRDS("./scratch/redTime/redTimeMVExp1/Exp1_mlegp/Exp1_mlegp_completed.rds")
    mdf <- mlegp.exp$outcleandf
    mdf$outdim <- rep(1:100, 8)
    mdf$use_PCAchar <- c("noPCA","PCA")[mdf$use_PCA+1]
    mdf$sopdchar <- c("1opd", "sopd")[mdf$separateoutputparameterdimensions+1]
    
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
    
  }
}