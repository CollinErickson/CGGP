# This file is run first to create the SGGP object
# that will be used for the rest of the simulation.

# Run this file on Crunch command line. It begins there,
# not with a submission. Run:
# Rscript SGGP_redTime_begin.R
# It will begin it and submit the necessary jobs after.

# Begin SGGP redTime experiment
cat("Starting SGGP_redTime_begin.R\n")

# This file will reload number_cores, d, SGGP_RDS_path, ???
source("SGGP_redTime_parameters.R")
cat("Loading SGGP_redTime_parameters.R successfully\n")

devtools::build()
cat("devtools::build() ran successfully\n")

install.packages(repos=NULL, pkgs="/home/collin/scratch/SGGP_1.0.tar.gz")
cat("SGGP installed correctly\n")

library('SGGP')
cat("SGGP loaded successfully\n")

# Create SGGP object
SG <- SGGPcreate(d=d, batchsize=N0, corr=corr, grid_sizes=grid_sizes)
print(SG)
cat("SGGPcreate successful\n")

# Save SG, nothing evaluated yet
saveRDS(object = SG, file = paste0(SGGP_after_append_RDS_path))
cat("saveRDS successful\n")

# Not able to fit until SG data is added, should be fixed
# # If Xsup is already given, fit to it.
# SG <- SGGPfit(SG, Ynew=Ynew, use_PCA=use_PCA,
#               Xs=Xs, Ys=Ys,
#               separateoutputparameterdimensions=separateoutputparameterdimensions)

# write params, write .sh, qsub, and prepare next R script
source("SGGP_redTime_qsub_unevaluated.R")
cat("qsub_unevaluated.R succesful\n")

# Start master job on crunch.local
system(paste("chmod 777 ", "SGGP_redTime_mastersub.sh"))
system('qsub -N SGGPmaster -l h=crunch SGGP_redTime_mastersub.sh')
cat("mastersub.sh submitted successfully\n")
cat(timestamp(), '\n')
