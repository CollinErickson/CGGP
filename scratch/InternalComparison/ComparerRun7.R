### ComparerRun7

# Created 6/4/19
# We found that power exp was the only corr func
# that didn't bounce up at end in ExternalComp 5/6.
# Then we found it had a nug, while the Cauchy funcs
# didn't. But Matern3/2 had a nug and still bounced.
# The purpose of InternalComp7 is to run the corr funcs
# with and without nug to see if nug affects the bounce.
# I made each corr func have no nug, then added a version
# with a nug, same name with Nug at the end.

sggpexp_func <- function(corr, sel.method, f, d, batchsize, 
                         sup.method, Nsupppd,
                         replicate) {#browser()
  require("CGGP")
  f <- eval(parse(text=paste0("TestFunctions::", f)))
  set.seed(0)
  nval <- 1e4
  Xval <- matrix(runif(nval*d), nval, d)
  Yval <- apply(Xval, 1, f)
  nmax <- 10000 * d#8192 #16384
  # Reset seed so you don't get repeats
  set.seed(Sys.time())
  
  # grid_size <- if (grid_size=="fast") {c(1, 2, 4, 4, 4, 6, 8, 32)}
  # grid_size <- if (grid_size=="fast") {c(1,2,4,4,8,12,32)}
  # else if (grid_size=="slow") {c(1, 2, 2, 2, 4, 4, 4, 4, 4, 6, 32)}
  # else if (grid_size=="medium") {c(1, 2, 4, 4, 4, 6, 8, 32)}
  # else {stop("bad grid size")}
  
  Ns <- d*Nsupppd
  if (Ns>0) {
    Xs <- lhs::maximinLHS(Ns, d)
    Ys <- apply(Xs, 1, f)
  } else {
    Xs <- NULL
    Ys <- NULL
  }
  
  # is_power_of_2 <- function(x) {
  #   abs(log(x, 2) - floor(log(x, 2)+.5)) < 1e-8
  # }
  
  start.time <- Sys.time()
  
  # First batch is zero points if supp data available
  batchsize0 <- if (Ns>0) 0 else batchsize
  
  sg <- CGGPcreate(d, batchsize0, corr=corr,
                   Xs=Xs, Ys=Ys, HandlingSuppData = sup.method)
  if (Ns == 0) {
    y <- apply(sg$design, 1, f)
    sg <- CGGPfit(sg, Y=y)
  }
  nallotted <- batchsize0
  
  if (sg$HandlingSuppData != sup.method) {stop(paste("HandlingSuppData is wrong"))}
  
  sg.stats <- NULL
  
  # while(nallotted < 1024) {
  #   sg <- CGGPappend(sg, batchsize=batchsize, selectionmethod=sel.method)
  #   y <- apply(sg$design, 1, f)
  #   sg <- CGGPfit(sg, Y=y,
  #                 Xs=Xs, Ys=Ys, HandlingSuppData=sup.method)
  #   if (sg$HandlingSuppData != sup.method) {stop(paste("HandlingSuppData is wrong"))}
  #   nallotted <- nallotted + batchsize
  #   if (nallotted >= 256 && is_power_of_2(nallotted)) {
  #     newstats <- CGGPvalstats(sg, Xval = Xval, Yval = Yval)
  #     newstats$nallotted <- nallotted
  #     newstats$ngrid <- nrow(sg$design)
  #     newstats$nsupp <- if (is.null(sg$Xs)) {0} else {nrow(sg$Xs)}
  #     newstats$ntotal <- nrow(sg$design) + if (!is.null(sg$Xs)) {nrow(sg$Xs)} else {0}
  #     newstats$elapsedtime <- as.numeric(Sys.time() - start.time, units='secs')
  #     sg.stats <- if (is.null(sg.stats)) newstats else rbind(sg.stats, newstats)
  #   }
  #   print(nallotted)
  # }
  # Aafter 1024 add in batches of larger size
  # batchsize2 <- 512 # was 512
  while(nallotted < nmax) { #16384) {
    
    Nalready <- (if (!is.null(sg$design)) {nrow(sg$design)} else {0}) + nrow(sg$Xs)
    # ni <- if (Nalready<1000) {200} else if (Nalready<10000) {500} else if (Nalready<20000) {2000} else {10000}
    batchsizei <- if (Nalready<1000) {200} else if (Nalready<4000) {500} else if (Nalready<20000) {2000} else {10000}
    sg <- CGGPappend(sg, batchsize=batchsizei, selectionmethod=sel.method)
    y <- apply(sg$design, 1, f)
    sg <- CGGPfit(sg, Y=y,
                  Xs=Xs, Ys=Ys, HandlingSuppData=sup.method)
    if (sg$HandlingSuppData != sup.method) {stop(paste("HandlingSuppData is wrong"))}
    nallotted <- nallotted + batchsizei
    if (TRUE) { # Always save  #nallotted >= batchsize2 && is_power_of_2(nallotted)) {
      newstats <- CGGPvalstats(sg, Xval = Xval, Yval = Yval)
      newstats$nallotted <- nallotted
      newstats$ngrid <- nrow(sg$design)
      newstats$nsupp <- if (is.null(sg$Xs)) {0} else {nrow(sg$Xs)}
      newstats$ntotal <- nrow(sg$design) + if (!is.null(sg$Xs)) {nrow(sg$Xs)} else {0}
      newstats$elapsedtime <- as.numeric(Sys.time() - start.time, units='secs')
      sg.stats <- if (is.null(sg.stats)) newstats else rbind(sg.stats, newstats)
      
    }
    print(nallotted)
  }
  sg.stats
}


require("comparer")

# The supp options when it is 0 are a waste, they are same as ignore

# sup.df <- data.frame(sup.method=rep(
#   c("Ignore", "Only", "Correct","Mixture", "MarginalValidation","FullValidation"), 4),
#   Nsupppd=rep(c(0,5,10,20), each=6), stringsAsFactors=F)
# sup.df <- sup.df[!(sup.df$Nsupppd==0 & sup.df$sup.method!="Ignore"),]

e2 <- ffexp$new(
  eval_func = sggpexp_func,
  # corr = c("cauchysq", "powerexp")[1], #"cauchysqt", "gaussian", "powerexp", "cauchy", "cauchysq"), #, "m32", "m52", "cauchysq", "cauchy"),
  corr = c("cauchysq", "cauchysqt", "cauchy", "powerexp", "m32", "m52", #"gauss",
           "cauchysqnug", "cauchysqtnug", "cauchynug", "powerexpnug", "m32nug", "m52nug", "gaussnug"),
  sel.method = c("Greedy"), #c("TS", "UCB", "Greedy"),
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"), d=c(3,6,7,8,10),
                row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
  batchsize=128,
  # sup.method=c("Ignore", "Only", "Correct","Mixture", "MarginalValidation","FullValidation"),
  # Nsupppd=c(0,5,10,20),
  sup.method="Correct",
  Nsupppd=10,
  # sup.df=sup.df,
  parallel=if (version$os =="linux-gnu") {TRUE} else {FALSE},
  parallel_cores = if (version$os =="linux-gnu") {20} else {3},
  replicate=1:3,
  folder_path= if (version$os =="linux-gnu") {"/home/collin/scratch/SGGP/scratch/InternalComparison/ComparerRun7/"}
               else {"./scratch/InternalComparison/ComparerRun7"}
)

e2$rungrid
# try because it gave delete error before, but shouldn't need it now
try(e2$recover_parallel_temp_save(delete_after = FALSE))
e2$save_self()
# e2$run_one(1)
# if (F) {
# e2$run_all(parallel_temp_save=TRUE, delete_parallel_temp_save_after=FALSE,
#            write_start_files=!TRUE, write_error_files=!T, run_order = "random")
e2$run_all(parallel_temp_save=TRUE, delete_parallel_temp_save_after=FALSE,
           write_start_files=TRUE, write_error_files=T)
# }
e2$recover_parallel_temp_save(delete_after = F)
e2$save_self()

print("Completed all runs in ComparerRun7.R")

if (F) {
  # e2 <- readRDS("./scratch/InternalComparison/ComparerRun6_object_328_of_360.rds")
  e2 <- readRDS("C:/Users/cbe117/Documents/GitHub/CGGP/scratch/redTime/redTimeData/ComparerRun7_completed.rds")
  e2$completed_runs %>% table
  e2dfcomp <- e2$outcleandf[e2$completed_runs,]
  e2df <- e2dfcomp[!is.na(e2dfcomp$RMSE),]
  colnames(e2df)[1] <- "corr.func"
  e2df$sup.method[e2df$nsupp==0] <- "Ignore"
  
  library(ggplot2)
  ggplot(data=e2df, mapping=aes(nallotted, RMSE, color=as.factor(Nsupppd))) + geom_point(size=4) + facet_grid(f ~ sup.method, scales="free_y") + scale_y_log10()
  ggplot(data=e2df, mapping=aes(nallotted, RMSE, color=as.factor(Nsupppd))) + geom_point(size=4) + facet_grid(f ~ sup.method, scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=e2df, mapping=aes(nallotted, score, color=as.factor(Nsupppd))) + geom_point(size=4) + facet_grid(f ~ sup.method, scales="free_y")
  ggplot(data=e2df, mapping=aes(nallotted, CRPscore, color=as.factor(Nsupppd))) + geom_point(size=4) + facet_grid(f ~ sup.method, scales="free_y") + scale_y_log10()
  ggplot(data=e2df, mapping=aes(nallotted, elapsedtime, color=as.factor(Nsupppd))) + geom_point(size=4) + facet_grid(f ~ sup.method, scales="free_y") + scale_y_log10()
  
}

if (F) {
  e7 <- readRDS("./scratch/InternalComparison/ComparerRun7_completed.rds")
  require(ggplot2); require(dplyr); require(magrittr)
  e7$completed_runs %>% table
  e7df <- e7$outcleandf[e7$completed_runs, ]
  colnames(e7df)[1] <- "corr.func"
  e7df$RMSE %>% summary
  e7df$RMSE %>% is.na %>% summary
  e7df$nug <- sapply(e7df$corr.func, function(xx) {grepl("nug", xx)})
  e7df$corr.funcnn <- sapply(e7df$corr.func, function(xx) {stringr::str_remove(xx, "nug")})
  ggplot(data=e7df, mapping=aes(nallotted, RMSE, color=corr.funcnn)) + geom_point() + scale_y_log10() + facet_grid(. ~ nug)
  ggplot(data=e7df, mapping=aes(nallotted, RMSE, color=corr.funcnn)) + geom_point() + scale_y_log10() + facet_grid(nug ~ corr.funcnn)
  ggplot(data=e7df %>% mutate(RMSE=pmin(10, RMSE)), mapping=aes(nallotted, RMSE, color=corr.funcnn)) + geom_point() + scale_y_log10() + facet_grid(nug ~ corr.funcnn)
}