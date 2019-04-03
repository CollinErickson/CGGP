# Comparison for paper

# decentLHS <- sFFLHD::decentLHS
decentLHS <- function(n, d, ndes, max.time) {
  if (missing(ndes) && missing(max.time)) {
    stop("Must give ndes or max.time to decentLHS")
  }
  
  # MaxProMeasure in 10d with 3e4 pts takes 3 min, crashes with 1e5 pts.
  if (n > 1e4) {
    return(lhs::randomLHS(n, d))
  }
  
  start.time <- Sys.time()
  bestdes <- NULL
  bestcrit <- Inf
  i <- 1
  while (TRUE) {
    # Make new LHS
    x <- lhs::randomLHS(n, d)
    crit <- MaxPro::MaxProMeasure(x)
    
    # Check if best yet
    if (crit < bestcrit) {
      bestcrit <- crit
      bestdes <- x
    }
    
    # Check if done
    i <- i+1
    if (!missing(ndes) && (i > ndes)) {
      break
    }
    if (!missing(max.time) && (as.numeric(Sys.time() - start.time, units="secs") > max.time)) {
      break
    }
  }
  # Return bestdes
  bestdes
}


run_BASS <- function(Ntotal, Nappend, f, d, x, y, xtest, ytest, seed) {
  if (!missing(seed)) {set.seed(seed)}
  if (missing(x) && missing(y)) {
    if (Ntotal<=2000) {x <- lhs::maximinLHS(Ntotal, d)}
    else {x <- decentLHS(Ntotal, d, max.time=15)}
    y <- apply(x, 1, f)
  }
  fit.time.start <- Sys.time()
  mod <- BASS::bass(xx=x, y=y)
  fit.time.end <- Sys.time()
  pred.time.start <- Sys.time()
  pred <- predict(mod, xtest)
  pred.time.end <- Sys.time()
  list(mean=colMeans(pred), var=rep(NaN, nrow(xtest)), n=nrow(x),
       pred.time=as.numeric(pred.time.end - pred.time.start, units="secs"),
       fit.time =as.numeric(fit.time.end  - fit.time.start , units="secs"))
}




# Need a generic function that passes to specific ones
run_one <- function(package, selection.method, correlation, HandlingSuppData,
                    f, d, npd, replicate) {#browser()
  # package <- psch$package
  n <- npd * d
  # if (n!= 500) {stop('bad n')}
  f <- eval(parse(text=paste0("TestFunctions::", f)))
  
  ntest <- 1e4
  # xtest <- matrix(runif(ntest*d), ncol=d)
  xtest <- decentLHS(ntest, d, max.time = 10)
  ytest <- apply(xtest, 1, f)
  
  # if (package == "CGGP") {
  #   out <- run_CGGP(Nappend=floor(n * (1:5)/5), f=f, d=d, xtest=xtest)
  if (package == "CGGPsupp") {
    # out <- run_CGGP(Nappend=floor(n*(2:5)/5), Nlhs=floor(.2*n), f=f, d=d, xtest=xtest)
    out <- run_CGGPsupp(Ntotal=n, Nlhs=10*d, f=f, d=d, xtest=xtest, HandlingSuppData=HandlingSuppData, selection.method=selection.method, correlation=correlation)
    # } else if (package == "CGGPsupponly") {
    #   out <- run_CGGP(Nappend=c(), Nlhs=n, f=f, d=d, xtest=xtest)
  } else if (package == "CGGPsupponly") {
    out <- run_CGGPsupponly(Ntotal=n, f=f, d=d, xtest=xtest, correlation=correlation)
  } else if (package == "CGGPoneshot") {
    out <- run_CGGPoneshot(Ntotal=n, f=f, d=d, xtest=xtest, correlation=correlation)
  } else if (package == "CGGP") {
    out <- run_CGGP(Ntotal=n, f=f, d=d, xtest=xtest, selection.method=selection.method, correlation=correlation)
  } else if (package == "laGP") {
    out <- run_lagp(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "aGP") {
    out <- run_lagp_bobby(Ntotal=n, f=f, d=d, xtest=xtest, use_agp=TRUE)
  } else if (package == "MRFA") {
    out <- run_MRFA(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "svm") {
    out <- run_svm(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "mlegp") {
    out <- run_mlegp(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "GPfit") {
    out <- run_GPfit(Ntotal=n, f=f, d=d, xtest=xtest)
  } else if (package == "BASS") {
    out <- run_BASS(Ntotal=n, f=f, d=d, xtest=xtest)
  } else {
    stop(paste("Package", package, "not recognized"))
  }
  # browser()
  outstats <- CGGP::valstats(predmean=out[[1]], predvar=out[[2]],Yval=ytest) #, bydim=FALSE)
  if (out$n > n) {warning(paste("n too big for", package, n, f, d))}
  outstats
}


expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

# eg1 <- expand.grid(selection.method = c("UCB", "Greedy"),
#                    correlation = c("CauchySQ", "Matern32", "PowerExp"), stringsAsFactors=F)
# eg2a <- expand.grid.df(eg1, data.frame(HandlingSuppData=c("Ignore", "Correct"), stringsAsFactors=F), data.frame(package="CGGPsupp", stringsAsFactors=F))
# eg2b <- expand.grid.df(eg1, data.frame(HandlingSuppData="NA", stringsAsFactors=F), data.frame(package=c("CGGP"), stringsAsFactors=F))
# eg2c <- expand.grid(selection.method="NA", correlation = c("CauchySQ", "Matern32", "PowerExp"), HandlingSuppData="NA", package=c("CGGPoneshot", "CGGPsupponly"), stringsAsFactors=F)
# eg2d <- data.frame(selection.method="NA", correlation="NA", HandlingSuppData="NA",
#                    package=c("MRFA", "svm", "aGP", "laGP", "mlegp", "GPfit"), stringsAsFactors=F)
# eg3 <- rbind(eg2a, eg2b, eg2c, eg2d)
eg4 <- data.frame(selection.method="NA", correlation="NA", HandlingSuppData="NA",
                  package=c("BASS"), stringsAsFactors=F)

require("comparer")

excomp <- ffexp$new(
  eval_func = run_one,
  varlist = c("decentLHS", "run_BASS"),
  fd=data.frame(f=c("beambending","OTL_Circuit","piston","borehole","wingweight"),
                d=c(3,6,7,8,10),
                row.names = c("beam","OTL","piston","borehole","wingweight"), stringsAsFactors = F),
  # package=c("CGGPsupp", "CGPPoneshot", "CGGPsupponly", "CGGP",
  #           "MRFA", "svm", "aGP", "laGP", "mlegp", "GPfit"),
  psch=eg4,
  npd=c(10, 30, 100, 300, 1000, 3000, 10000),
  parallel=T,
  parallel_cores = 20,
  replicate=1:10, #:5,
  folder_path= "/home/collin/scratch/SGGP/scratch/ExternalComparison/ExComp4bass"
  # folder_path="./scratch/ExternalComparison/ExComp4bass/"
)
# Remove ones that can't do full size
package.name <- excomp$arglist$psch$package[excomp$rungrid$psch]
npd.excomp <- excomp$arglist$npd[excomp$rungrid$npd]
n.excomp <- npd.excomp * excomp$arglist$fd$d[excomp$rungrid$fd]
# excomp$completed_runs[package.name == "CGGPsupponly" & n.excomp > 1000] <- TRUE
# excomp$completed_runs[package.name == "laGP" & n.excomp > 1000] <- TRUE
# excomp$completed_runs[package.name == "mlegp" & n.excomp > 400] <- TRUE
# excomp$completed_runs[package.name == "GPfit" & n.excomp > 100] <- TRUE
# # agp is giving errors
# excomp$completed_runs[package.name == "aGP"] <- TRUE
table(paste(package.name, excomp$completed_runs))
# plyr::ddply(data.frame(package.name, compruns=excomp$completed_runs), "package.name", function(d) {data.frame(done=sum(d$compruns), notdone=sum(!d$compruns))})
rbind(plyr::ddply(data.frame(package.name, compruns=excomp$completed_runs), "package.name", function(d) {data.frame(done=sum(d$compruns), notdone=sum(!d$compruns))}), data.frame(package.name="all", done=sum(excomp$completed_runs), notdone=sum(!excomp$completed_runs)))


# excomp$run_one(1752)
# excomp$run_one(7299)
# excomp$run_all(save_output = T, parallel = F, parallel_temp_save = T, run_order = "random")

# excomp$rungrid
# try because it gave delete error before, but shouldn't need it now
table(excomp$completed_runs)
try(excomp$recover_parallel_temp_save(delete_after = FALSE))
table(excomp$completed_runs)
# plyr::ddply(data.frame(package.name, compruns=excomp$completed_runs), "package.name", function(d) {data.frame(done=sum(d$compruns), notdone=sum(!d$compruns))})
rbind(plyr::ddply(data.frame(package.name, compruns=excomp$completed_runs), "package.name", function(d) {data.frame(done=sum(d$compruns), notdone=sum(!d$compruns))}), data.frame(package.name="all", done=sum(excomp$completed_runs), notdone=sum(!excomp$completed_runs)))
excomp$save_self()
excomp$run_all(parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE,
               write_start_files=T, write_error_files=T)
# # Getting errors, run by package.name to see which is causing it
# excomp$parallel_cores <- 10
# excomp$run_all(to_run = which(!excomp$completed_runs & (package.name == "mlegp")),
#                parallel_temp_save = TRUE, delete_parallel_temp_save_after=FALSE,
#                write_start_files=T, write_error_files=T)
# excomp$run_all()

cat("Completed all runs in ExternalComparer4.R\n")

excomp$save_self()

cat("Saved self\n")

if (F) {
  excomp$plot_run_times()
  plyr::dlply(excomp$outcleandf, "d")
  require('ggplot2')
  ecdf <- excomp$outcleandf[excomp$completed_runs,]
  ecdf$n <- ecdf$npd * ecdf$d
  ggplot(data=ecdf, mapping=aes(n, RMSE, color=as.factor(n))) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  ggplot(data=ecdf, mapping=aes(n, score, color=as.factor(n))) + geom_point() + facet_grid(f ~ package, scales="free_y")
  ggplot(data=ecdf, mapping=aes(n, CRPscore)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  ggplot(data=ecdf, mapping=aes(n, runtime)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10() + scale_x_log10()
  # saveRDS(excomp, "./scratch/ExternalComparison/ExComp1_completed.rds")
}
