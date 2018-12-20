# Begin SGGP redTime experiment

# This file will reload number_cores, d, SGGP_RDS_path, ???
source("SGGP_redTime_parameters.R")

SG <- SGGPcreate(d=d, batchsize=N0, corr=corr)

# Save SG, nothing evaluated yet
saveRDS(object = SG, file = paste0(SGGP_after_append_RDS_path))

# write params, write .sh, qsub, and prepare next R script
source("SGGP_redTime_qsub_unevaluated.R")
