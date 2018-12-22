# Continue redTime experiment
# This runs after batch of points are evaluated.
# It is in charge of fitting, saving object, checking if done, appending and submitting next batch if not done.
cat("Starting SGGP_redTime_continue.R\n")

# This file will reload number_cores, d, SGGP_RDS_path, ???
# WILL NEED TO LET THIS FILE NAME CHANGE LATER
source("/home/collin/scratch/SGGP/scratch/redTime/SGGP_redTime_parameters.R")
cat("Reloaded parameters successfully\n")

# Function to check if all runs are completed
# This is terrible method now, since unrelated jobs of mine could still be running
all_done <- function() {
  # If this is only job running, then qstat will have three lines, something like
  # job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
  # -----------------------------------------------------------------------------------------------------------------
  #   168501 0.55500 testsub.sh collin       r     12/21/2018 09:00:20 all.q@crunch.local                 1
  
  timestamp()
  qstatout <- system('qstat', intern=TRUE)
  print(qstatout)
  length(qstatout) <= 3
}

# Check if all evaluations are done, if not, sleep until it is
while (!all_done()) {
  # system('sleep 100')
  Sys.sleep(10)
  # cat("Still files running\n")
  # timestamp()
  # cat(system('qstat', intern=TRUE))
}

if (!all_done()) {
  stop("Not possible, must be all_done")
}

cat("All other jobs done, moving on\n")

# Load SGGP object, the one that had new design points appended but not measured yet
SG <- readRDS(SGGP_after_append_RDS_path)
cat("Read back in RDS successfully\n")

# Read in all new output values
print(getwd())
source(paste0(sourcefilepath, "extract_redTime.R"))
Ynew <- NULL
for (i in 1:nrow(SG$design_unevaluated)) {
  newrow <- extract_redTime_output(outpath = paste0(outpathbase, "_", SG$ss, "_", i, ".out"))
  if (is.null(Ynew)) {Ynew <- newrow}
  else {Ynew <- rbind(Ynew, newrow)}
}
rownames(Ynew) <- NULL
cat("Extracted all values\n", Ynew, "\n")
cat("Class of Ynew is ", class(Ynew), " dim of Ynew is ", dim(Ynew), "\n")

# Update params with new data
library("SGGP")
SG <- SGGPfit(SG, Ynew=Ynew)
cat("SGGPfit successful\n")

# Save SGGP object
saveRDS(object = SG, file = paste0(SGGP_after_fit_RDS_path))
cat("Saved SG successfully\n")

# Check if done
if (SG$ss >= Nfinal) {
  cat("Simulation completed!")
  stop("All done")
}

# Append new values
SG <- SGGPappend(SG, batchsize)
cat("SGGPappend successful\n")

# Save version with appended design
saveRDS(object = SG, file = paste0(SGGP_after_append_RDS_path))
cat("saveRDS successful\n")

# write params, write .sh, qsub, and prepare next R script
source(paste0(sourcefilepath, "SGGP_redTime_qsub_unevaluated.R"))
cat("SGGP_redTime_qsub_unevaluated.R successful, this ends SGGP_redTime_continue.R\n")

# # Start master now
# source()
# cat("Master has completed\n")
# cat(timestamp(), '\n')
