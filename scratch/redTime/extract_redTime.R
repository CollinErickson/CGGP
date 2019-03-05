
extract_redTime_output <- function(fileID, outpath=NULL) {
  if (is.null(outpath)) {
    outpath <- paste0("/home/collin/scratch/redTime_v0.1/output_files/out_", fileID, ".out")
  }
  # outpath <- "../../../Downloads/out_LHS1LCompv1_1L_1rs_15.out"
  
  # First find where data part starts
  rl <- readLines(outpath)
  mainline <- which(sapply(rl, function(x) {substring(x, 1, 16) == "### main: output"}))
  if (length(mainline) != 1) {stop(paste(c("mainline length is not 1", outpath, mainline)))}
  
  # Need to read in relevant part of file as csv, it is space/tab separated
  df <- read.csv(outpath, header = FALSE, sep = "", nrows = 100, skip = mainline)
  
  # We have the 8th column, which is P_dd = Time-RG density-density power spectrum
  df[,8]
}

extract_redTime_from_completed_LHS <- function(n, prefix='') {
  # Will be DF to write out
  Ynew <- NULL
  
  for (i in 1:n) {
    # Get row from single LHS run
    # outputpath <- paste0(outpathbase, fileID, ".out")
    newrow <- extract_redTime_output(fileID= paste0(prefix, i)) #outpath=outpath)
    
    # Add it to df
    if (is.null(Ynew)) {Ynew <- newrow}
    else {Ynew <- rbind(Ynew, newrow)}
  }
  
  rownames(Ynew) <- NULL
  
  outpath <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/important_files/", prefix, "all_output.csv")
  write.csv(x=Ynew, file = outpath)
}

extract_redTime_from_SGGP_run <- function(prefix, Nmax) {
  # This function will get all redTime output for a SGGP run.
  # If the SGGP object isn't saved, or if it only keeps certain dimensions,
  #  then this will recover everything from the .out files.
  # Set Nmax if currently running so it doesn't try to get files
  #  that are still being run.
  # 
  # prefix="redTimeTestSup2o50"
  
  # Get all output files
  lfout <- list.files(paste0("~/scratch/redTime_v0.1/SGGPruns/", prefix,"/output_files/"), pattern=".out")
  # Param files
  # lfpar <- list.files(paste0("~/scratch/redTime_v0.1/SGGPruns/", prefix,"/param_files/"), pattern=".dat")
  
  # if (length(lfout) != length(lfpar)) {
  #   stop(paste("lfout and lfpar have different lengths", length(lfout), length(lfpar)))
  # }
  
  # dfpar <- cbind(lfpar, do.call(rbind, lapply(strsplit(unlist(strsplit(lfpar, ".dat")), "_"), function (ss) data.frame(N1=as.integer(ss[length(ss)-1]), N2=as.integer(ss[length(ss)])))))
  dfout <- cbind(lfout, do.call(rbind, lapply(strsplit(unlist(strsplit(lfout, ".out")), "_"), function (ss) data.frame(N1=as.integer(ss[length(ss)-1]), N2=as.integer(ss[length(ss)])))), stringsAsFactors=F)
  # dfall <- merge(dfpar, dfout, by=c("N1", "N2"))
  dfall <- dfout
  rm(dfout, lfout)
  
  
  # if ((nrow(dfpar) != nrow(dfall)) || (nrow(dfout) != nrow(dfall))) {
  #   stop("Problem merging, losing rows")
  # }
  dfall2 = plyr::ddply(dfall, "N1", function(x) {x$N2max = max(x$N2); x})
  rm(dfall)
  dfall2$iii = dfall2$N1-dfall2$N2max + dfall2$N2
  
  
  # source("extract_redTime.R")
  
  # Want it in order by iii
  dfall3 <- dfall2[order(dfall2$iii),]
  rm(dfall2)
  # Want to only keep those under or at Nmax
  dfall3 <- dfall3[dfall3$iii<=Nmax,]
  if (max(dfall3$iii) != nrow(dfall3)) {stop("Not all rows here???")}
  
  Ynew <- matrix(NaN, nrow(dfall3), 100)
  
  for (i in 1:nrow(dfall3)) {
    # Get row from single LHS run
    # outputpath <- paste0(outpathbase, fileID, ".out")
    newoutput <- extract_redTime_output(
      outpath=paste0("~/scratch/redTime_v0.1/SGGPruns/",
                     prefix,"/output_files/", dfall3$lfout[i])) #outpath=outpath)
    
    # # Add it to df
    # if (is.null(Ynew)) {Ynew <- newrow}
    # else {Ynew <- rbind(Ynew, newrow)}
    Ynew[dfall3$iii[i],] <- newoutput
  }
  
  rownames(Ynew) <- NULL
  
  outpath <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/important_files/", prefix, "_all_SGGP_output.csv")
  if (file.exists(outpath)) {stop(paste("outpath already exists:", outpath))}
  write.csv(x=Ynew, file = outpath)
}
