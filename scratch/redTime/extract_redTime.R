
extract_redTime_output <- function(fileID, outpath=NULL) {
  if (is.null(outpath)) {
    outpath <- paste0("/home/collin/scratch/redTime_v0.1/output_files/out_", fileID, ".out")
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
  for (i in 1:n) {
    # outputpath <- paste0(outpathbase, fileID, ".out")
    extract_redTime_output(fileID= paste0(prefix, i)) #outpath=outpath)
  }
}
