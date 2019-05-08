# This will continue a redTime CGGP experiment that completed
#  but you want to run further.
# This is similar to _restart and _continue
# Make sure you have increased Nfinal in _parameters


# This file will reload number_cores, d, SGGP_RDS_path, parameters, etc
source("SGGP_redTime_parameters.R")
cat("Loading SGGP_redTime_parameters.R successfully\n")

cat("Assuming redTime param files were all completed and added to model already\n")
cat("\n")

devtools::build()
cat("devtools::build() ran successfully\n")

install.packages(repos=NULL, pkgs="/home/collin/scratch/CGGP_1.0.tar.gz")
cat("CGGP installed correctly\n")

library('CGGP')
cat("CGGP loaded successfully\n")

# Load SGGP object, the one that had completed, no points unfitted
SG <- readRDS(SGGP_after_fit_RDS_path)
cat("Read back in RDS successfully, has ", nrow(SG$design), "points\n")
if (nrow(SG$design) != nrow(SG$Y)) {
  stop("Not all points have been evaluated")
} else {
  cat("All points have been evaluated, looks good\n")
}


# Check if done
if (SG$ss >= Nfinal) {
  cat("Simulation completed!")
  stop("All done")
}

# Set seed here to check refitting on laptop
set.seed(SG$ss)

# Append new values
# batchsize <- if (SG$ss<1000) {batchsize1} else if (SG$ss<10000) {batchsize2} else {batchsize3}
batchsize <- if (SG$ss<1000) {batchsize1} else if (SG$ss<10000) {batchsize2} else if (SG$ss<20000) {batchsize3} else {batchsize4}
SG <- CGGPappend(SG, batchsize, selectionmethod=selectionmethod)
cat("CGGPappend successful, now has ", nrow(SG$design), "points\n")

# Save version with appended design
saveRDS(object = SG, file = paste0(SGGP_after_append_RDS_path))
cat("saveRDS successful\n")

# write params, write .sh, qsub, and prepare next R script
source(paste0(sourcefilepath, "SGGP_redTime_qsub_unevaluated.R"))
cat("SGGP_redTime_qsub_unevaluated.R successful, this ends SGGP_redTime_restart_completed.R\n")



# Start master job on crunch.local. This will pick up next iteration
system(paste("chmod 777 ", "SGGP_redTime_mastersub.sh"))
system('qsub -N SGGPmaster -l h=crunch SGGP_redTime_mastersub.sh')
cat("mastersub.sh submitted successfully\n")
cat(timestamp(), '\n')
