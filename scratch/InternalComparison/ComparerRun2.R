# I'm trying to see if my function comparer::ffexp can be used to run our
# internal comparison experiment. It should make parallelizing and sending
# out trials to be run easy.

# This below is just a small run on a single function.
# I commented out the part that runs it to 16384 so it will stop at 1024.
# So it may take 10x longer, but it fits less often,
# so it shouldn't be that bad.

sggpexp_func <- function(corr, sel.method, f, d, batchsize, pred.fullBayes,
                         append.rimseperpoint, use_laplaceapprox, grid_size) {
  require("SGGP")
  # f <- TestFunctions::borehole
  # d <- 8
  # browser()
  f <- eval(parse(text=paste0("TestFunctions::", f)))
  set.seed(0)
  nval <- 1e4
  Xval <- matrix(runif(nval*d), nval, d)
  Yval <- apply(Xval, 1, f)

  # n0 <- 64
  # corr <- SGGP_internal_CorrMatCauchySQT
  # sel.method <- "TS"
  # batchsize <- 64
  # pred.fullBayes <- TRUE
  # append.rimseperpoint <- TRUE
  # use_laplaceapprox <- TRUE
  # grid_size <- c(1, 2, 2, 2, 4, 4, 4, 4, 4, 6, 32)
  grid_size <- if (grid_size=="fast") {c(1, 2, 4, 4, 4, 6, 8, 32)}
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



e2 <- ffexp$new(
  eval_func = sggpexp_func,
  corr = c("cauchysqt", "gaussian", "powerexp", "m32", "m52", "cauchysq", "cauchy"),
  sel.method = c("TS", "UCB", "Greedy"),
  fd=data.frame(f="borehole", d=8, row.names = c("borehole"), stringsAsFactors = F),
  batchsize=c(64, 256),
  append.rimseperpoint=c(TRUE, FALSE),
  use_laplaceapprox=c(TRUE), #c(TRUE, FALSE),
  pred.fullBayes=c(FALSE), #c(TRUE, FALSE),
  grid_size=c("fast", "slow"),
  parallel=TRUE,
  parallel_cores = 'detect-1',
  folder_path="./scratch/sggpout"
)
e2$rungrid
if (F) {
  e2$run_all(parallel_temp_save = TRUE)
}
