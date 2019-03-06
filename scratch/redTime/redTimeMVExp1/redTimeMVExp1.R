# Running an experiment to test redTime with all 100 outputs.
# we have 2x2 options.

evfunc <- function(N, use_PCA, separateoutputparameterdimensions) {
  if (version$os=="linux-gnu") {
    x1000 <- unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
    y1000 <- log(unname(as.matrix(read.csv("~/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
    sg.base <- readRDS(paste0("~/scratch/redTime_v0.1/SGGPruns/redTimeTestSup2o50/output_files/out_S2o50_SGGP-", N, ".rds"))
    
    
  } else if (version$os == "mingw32") {
    x1000 <- unname(as.matrix(read.csv("./scratch/redTime/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_input.csv")[,-1]))
    y1000 <- log(unname(as.matrix(read.csv("./scratch/redTime/redTimeData/redTimeData/ExpandedRanges2_LHS1L_n1000_s0303_all_output.csv")[,-1])))
    sg.base <- readRDS(paste0("./scratch/redTime/redTimeData/redTimeData/out_S2o50_SGGP-", N, ".rds"))
  } else {
    stop("version$os doesn't match any")
  }
  
  
  sg <- SGGPfit(sg.base, Y=Y.6889[1:N,], use_PCA = use_PCA,
                separateoutputparameterdimensions = separateoutputparameterdimensions)
  
  SGGPvalstats(sg, x1000, y1000)
}


e1 <- comparer::ffexp$new(
  N = c(199, 299, 399, 499, 699, 1299, 1299, 2499, 4099, 6099, 8099),
  use_PCA = c(T,F),
  separateoutputparameterdimensions = c(T,F),
  eval_func = evfunc,
  folder_path = if (version$os=="linux-gnu") {"/home/collin/scratch/SGGP/scratch/redTime/redTimeMVExp1/Exp1"} else if (version$os=="mingw32") {"./scratch/InternalComparison/redTimeMVout_from_S2o50"} else {stop("bad folderpath")}
  # folder_path = "./scratch/InternalComparison/redTimeMVout_from_S2o50",
  parallel=T,
  parallel_cores=10
)
e1$recover_parallel_temp_save(delete_after = F)

e1$run_all(delete_parallel_temp_save_after = F, parallel_temp_save = T, write_start_files = T, write_error_files = T)

e1$save_self()

edf <- e1$outcleandf
edf$outdim <- rep(1:100, 32)

ggplot(data=edf, mapping=aes(x=N, y=RMSE, color=outdim,shape=interaction(use_PCA, separateoutputparameterdimensions))) + geom_point()
