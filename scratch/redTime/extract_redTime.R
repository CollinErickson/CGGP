
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

extract_redTime_from_SGGP_run <- function(n, prefix) {
  # 
  
  # Get all output files
  lfout <- list.files(paste0("~/scratch/redTime_v0.1/SGGPruns/", prefix,"/output_files/"), pattern=".out")
  # Param files
  lfpar <- list.files(paste0("~/scratch/redTime_v0.1/SGGPruns/", prefix,"/param_files/"), pattern=".dat")

  if (length(lfout) != length(lfpar)) {
    stop(paste("lfout and lfpar have different lengths", length(lfout), length(lfpar)))
  }
  
  dfpar <- cbind(lfpar, do.call(rbind, lapply(strsplit(unlist(strsplit(lfpar, ".dat")), "_"), function (ss) data.frame(N1=as.integer(ss[length(ss)-1]), N2=as.integer(ss[length(ss)])))))
  dfout <- cbind(lfout, do.call(rbind, lapply(strsplit(unlist(strsplit(lfout, ".out")), "_"), function (ss) data.frame(N1=as.integer(ss[length(ss)-1]), N2=as.integer(ss[length(ss)])))))
  dfall <- merge(dfpar, dfout, by=c("N1", "N2"))

  
  
  if ((nrow(dfpar) != nrow(dfall)) || (nrow(dfout) != nrow(dfall))) {
    stop("Problem merging, losing rows")
  }
  dfall2 = plyr::ddply(dfall, "N1", function(x) {x$N2max = max(x$N2); x})
  dfall2$iii = dfall2$N1-dfall2$N2max + dfall2$N2
  
  Ynew <- matrix(NaN, nrow(dfall2), 100)
  
  for (i in 1:nrow(dfall2)) {
    # Get row from single LHS run
    # outputpath <- paste0(outpathbase, fileID, ".out")
    newoutput <- extract_redTime_output(outpath = paste0("~/scratch/redTime_v0.1/SGGPruns/", prefix,"/output_files/", dfall2[i,lfpar],".out")) #outpath=outpath)
    
    # # Add it to df
    # if (is.null(Ynew)) {Ynew <- newrow}
    # else {Ynew <- rbind(Ynew, newrow)}
    Ynew[dfall2$iii[i],] <- newoutput
  }
  
  rownames(Ynew) <- NULL
  
  outpath <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/important_files/", prefix, "all_output.csv")
  if (file.exists(outpath)) {stop(paste("outpath already exists:", outpath))}
  write.csv(x=Ynew, file = outpath)
}
