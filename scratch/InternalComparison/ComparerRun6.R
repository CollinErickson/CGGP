### ComparerRun6

# Run with supp data, different supp methods

sggpexp_func <- function(corr, sel.method, f, d, batchsize, pred.fullBayes,
                         supppd, sup.method, Nsupppd,
                         replicate, # Not used, just to allow it to do multiples
                         append.rimseperpoint, use_laplaceapprox, grid_size) {browser()
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
  else if (grid_size=="medium") {c(1, 2, 4, 4, 4, 6, 8, 32)}
  else {stop("bad grid size")}
  
  Ns <- d*Nsupppd
  if (Ns>0) {
    Xs <- lhs::maximinLHS(Ns, d)
    Ys <- apply(Xs, 1, f)
  } else {
    Xs <- NULL
    Ys <- NULL
  }
  
  is_power_of_2 <- function(x) {
    abs(log(x, 2) - floor(log(x, 2)+.5)) < 1e-8
  }
  
  start.time <- Sys.time()
  
  # First batch is zero points if supp data available
  batchsize0 <- if (Ns>0) 0 else batchsize
  
  sg <- SGGPcreate(d, batchsize0, corr=corr,
                   Xs=Xs, Ys=Ys, HandlingSuppData = sup.method,
                   grid_sizes = grid_size)
  if (Ns == 0) {
    y <- apply(sg$design, 1, f)
    sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox)
  }
  nallotted <- batchsize0
  
  sg.stats <- NULL
  
  while(nallotted < 1024) {
    sg <- SGGPappend(sg, batchsize=batchsize, selectionmethod=sel.method, RIMSEperpoint = append.rimseperpoint)
    y <- apply(sg$design, 1, f)
    sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox,
                  Xs=Xs, Ys=Ys, HandlingSuppData=sup.method)
    nallotted <- nallotted + batchsize
    if (nallotted >= 256 && is_power_of_2(nallotted)) {
      newstats <- SGGPvalstats(sg, Xval = Xval, Yval = Yval, fullBayesian=pred.fullBayes)
      newstats$nallotted <- nallotted
      newstats$nused <- nrow(sg$design)
      newstats$nsupp <- nrow(sg$Xs)
      newstats$ntotal <- nrow(sg$design) + if (!is.null(sg$Xs)) {nrow(sg$Xs)} else {0}
      newstats$elapsedtime <- as.numeric(Sys.time() - start.time, units='secs')
      sg.stats <- if (is.null(sg.stats)) newstats else rbind(sg.stats, newstats)
    }
    print(nallotted)
  }
  # Aafter 1024 add in batches of larger size
  batchsize2 <- 512 # was 512
  while(nallotted < 16384) {
    sg <- SGGPappend(sg, batchsize=batchsize2, selectionmethod=sel.method, RIMSEperpoint = append.rimseperpoint)
    y <- apply(sg$design, 1, f)
    sg <- SGGPfit(sg, Y=y, laplaceapprox=use_laplaceapprox,
                  Xs=Xs, Ys=Ys, HandlingSuppData=sup.method)
    nallotted <- nallotted + batchsize2
    if (nallotted >= batchsize2 && is_power_of_2(nallotted)) {
      newstats <- SGGPvalstats(sg, Xval = Xval, Yval = Yval, fullBayesian=pred.fullBayes)
      newstats$nallotted <- nallotted
      newstats$nused <- nrow(sg$design)
      newstats$nsupp <- nrow(sg$Xs)
      newstats$ntotal <- nrow(sg$design) + if (!is.null(sg$Xs)) {nrow(sg$Xs)} else {0}
      newstats$elapsedtime <- as.numeric(Sys.time() - start.time, units='secs')
      sg.stats <- if (is.null(sg.stats)) newstats else rbind(sg.stats, newstats)
      
    }
    print(nallotted)
  }
  sg.stats
}


require("comparer")

e2 <- ffexp$new(
  eval_func = sggpexp_func,
  corr = c("cauchysq", "powerexp"), #"cauchysqt", "gaussian", "powerexp", "cauchy", "cauchysq"), #, "m32", "m52", "cauchysq", "cauchy"),
  sel.method = c("UCB"), #c("TS", "UCB", "Greedy"),
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"), d=c(3,6,7,8,10),
                row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
  batchsize=c(64),
  append.rimseperpoint=c(TRUE), #, FALSE),
  sup.method=c("Ignore", "Only", "Correct","Mixture", "MarginalValidation","FullValidation"),
  Nsupppd=c(0,5,10,20),
  use_laplaceapprox=c(TRUE), # No MCMC, too slow
  pred.fullBayes=c(FALSE), # No full Bayes prediction, too slow
  grid_size=c("medium"), #"fast", "slow"),
  parallel=TRUE,
  parallel_cores = 20,
  replicate=1:3,
  folder_path= "/home/collin/scratch/SGGP/scratch/InternalComparison/ComparerRun6" #"./scratch/sggpout"
  # folder_path="./scratch/delthisnow/"
)

e2$rungrid
# try because it gave delete error before, but shouldn't need it now
try(e2$recover_parallel_temp_save(delete_after = FALSE))
e2$save_self()
# if (F) {
e2$run_all(parallel_temp_save=TRUE, delete_parallel_temp_save_after=FALSE,
           write_start_files=TRUE)
# }
e2$save_self()

print("Completed all runs in ComparerRun6.R")