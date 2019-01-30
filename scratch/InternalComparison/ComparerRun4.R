### ComparerRun4

# This will include some slow options

sggpexp_func <- function(corr, sel.method, f, d, batchsize, pred.fullBayes,
                         append.rimseperpoint, use_laplaceapprox, grid_size) {
  require("SGGP")
  f <- eval(parse(text=paste0("TestFunctions::", f)))
  set.seed(0)
  nval <- 1e4
  Xval <- matrix(runif(nval*d), nval, d)
  Yval <- apply(Xval, 1, f)
  # Reset seed so you don't get repeats
  set.seed(Sys.time())
  
  # grid_size <- if (grid_size=="fast") {c(1, 2, 4, 4, 4, 6, 8, 32)}
  grid_size <- if (grid_size=="fast") {c(1,2,4,4,8,12,32)}
  else if (grid_size=="slow") {c(1, 2, 2, 2, 4, 4, 4, 4, 4, 6, 32)}
  else {stop("bad grid size")}
  
  
  
  
  is_power_of_2 <- function(x) {
    abs(log(x, 2) - floor(log(x, 2)+.5)) < 1e-8
  }
  
  start.time <- Sys.time()
  
  sg <- SGGPcreate(d, batchsize, corr=corr, grid_sizes = grid_size)
  y <- apply(sg$design, 1, f)
  sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox)
  nallotted <- batchsize
  
  sg.stats <- NULL
  
  while(nallotted < 1024) {
    sg <- SGGPappend(sg, batchsize=batchsize, selectionmethod=sel.method, RIMSEperpoint = append.rimseperpoint)
    y <- apply(sg$design, 1, f)
    sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox)
    nallotted <- nallotted + batchsize
    if (nallotted >= 512 && is_power_of_2(nallotted)) {
      newstats <- SGGPvalstats(sg, Xval = Xval, Yval = Yval, fullBayesian=pred.fullBayes)
      newstats$nallotted <- nallotted
      newstats$nused <- nrow(sg$design)
      newstats$elapsedtime <- as.numeric(Sys.time() - start.time, units='secs')
      sg.stats <- if (is.null(sg.stats)) newstats else rbind(sg.stats, newstats)
    }
    print(nallotted)
  }
  while(nallotted < 16384) {
    sg <- SGGPappend(sg, batchsize=512, selectionmethod=sel.method, RIMSEperpoint = append.rimseperpoint)
    y <- apply(sg$design, 1, f)
    sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox)
    nallotted <- nallotted + 512
    if (nallotted >= 512 && is_power_of_2(nallotted)) {
      newstats <- SGGPvalstats(sg, Xval = Xval, Yval = Yval, fullBayesian=pred.fullBayes)
      newstats$nallotted <- nallotted
      newstats$nused <- nrow(sg$design)
      newstats$elapsedtime <- as.numeric(Sys.time() - start.time, units='secs')
      sg.stats <- if (is.null(sg.stats)) newstats else rbind(sg.stats, newstats)

    }
    print(nallotted)
  }
  sg.stats
}


require("comparer")

# This was main one using five correlation functions
# e2 <- ffexp$new(
#   eval_func = sggpexp_func,
#   corr = c("cauchysqt", "gaussian", "powerexp", "cauchy", "cauchysq"), #, "m32", "m52", "cauchysq", "cauchy"),
#   sel.method = c("TS", "UCB", "Greedy"),
#   fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"), d=c(3,6,7,8,10),
#                 row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
#   batchsize=c(64, 256),
#   append.rimseperpoint=c(TRUE, FALSE),
#   use_laplaceapprox=c(TRUE), # No MCMC, too slow
#   pred.fullBayes=c(FALSE), # No full Bayes prediction, too slow
#   grid_size=c("fast", "slow"),
#   parallel=TRUE,
#   parallel_cores = 37,
#   folder_path= "/home/collin/scratch/SGGP/scratch/InternalComparison/ComparerRun4a" #"./scratch/sggpout"
# )

# Going to run again with other correlation functions
e2 <- ffexp$new(
  eval_func = sggpexp_func,
  corr = c("m32", "m52"), #"cauchysqt", "gaussian", "powerexp", "cauchy", "cauchysq"), #, "m32", "m52", "cauchysq", "cauchy"),
  sel.method = c("TS", "UCB", "Greedy"),
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"), d=c(3,6,7,8,10),
                row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
  batchsize=c(64, 256),
  append.rimseperpoint=c(TRUE, FALSE),
  use_laplaceapprox=c(TRUE), # No MCMC, too slow
  pred.fullBayes=c(FALSE), # No full Bayes prediction, too slow
  grid_size=c("fast", "slow"),
  parallel=TRUE,
  parallel_cores = 37,
  folder_path= "/home/collin/scratch/SGGP/scratch/InternalComparison/ComparerRun4b" #"./scratch/sggpout"
)

e2$rungrid
# try because it gave delete error before, but shouldn't need it now
try(e2$recover_parallel_temp_save(delete_after = FALSE))
e2$save_self()
# if (F) {
e2$run_all(parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE)
# }
e2$save_self()

print("Completed all runs in ComparerRun4.R")
