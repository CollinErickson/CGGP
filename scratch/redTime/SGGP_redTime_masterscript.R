# This is the master file that will always be running.
# It forces itself on crunch.local so it is able to qsub other jobs.
# If I were to stop it and restart it with holds, it might lose its place
#  to other jobs and be slow getting the next batch set up.
# The main problem is that qsub can only be run from crunch.local as far as I can tell.

cat("Starting SGGP_redTime_masterscript.R\n")

source("/home/collin/scratch/SGGP/scratch/redTime/SGGP_redTime_parameters.R")
cat("Reloaded parameters successfully in masterscript, beginning loop\n")


while (TRUE) {
  # Sleep 60 seconds between iterations to make sure they have time to show up in queue, could be less
  Sys.sleep(60)
  # Run the script
  source(paste0(sourcefilepath, "SGGP_redTime_continue.R"))
}