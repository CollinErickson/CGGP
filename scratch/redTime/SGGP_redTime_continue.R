# Continue redTime experiment
# This runs after batch of points are evaluated.
# It is in charge of fitting, saving object, checking if done, appending and submitting next batch if not done.

# This file will reload number_cores, d, SGGP_RDS_path, ???
source("SGGP_redTime_parameters.R")

# Function to check if all runs are completed
# This is terrible method now, since unrelated jobs of mine could still be running
all_done <- function() {
  timestamp()
  qstatout <- system('qstat', intern=TRUE)
  print(qstatout)
  length(qstatout) < 3
}

# Check if all evaluations are done, if not, sleep until it is
while (!all_done()) {
  system('sleep 100')
}

if (!all_done()) {
  stop("Not possible, must be all_done")
}

# Load SGGP object, the one that had new design points appended but not measured yet
readRDS(SGGP_after_append_RDS_path)

# Read in all new output values
print(getwd())
source("extract_redTime.R")
Ynew <- NULL
for (i in 1:nrow(SG$design_unevaluated)) {
  newrow <- extract_redTime_output(paste0(fileIDbase, "_", SG$ss, "_", i))
  if (is.null(Ynew)) {Ynew <- newrow}
  else {Ynew <- rbind(Ynew, newrow)}
}

# Update params with new data
SG <- SGGPfit(SG, Ynew=Ynew)

# Save SGGP object
saveRDS(object = SG, file = paste0(SGGP_after_fit_RDS_path))

# Check if done
if (SGGP$ss >= Nfinal) {
  stop("All done")
}

# Append new values
SG <- SGGPappend(SG, batchsize)

# Save version with appended design
saveRDS(object = SG, file = paste0(SGGP_after_append_RDS_path))


# write params, write .sh, qsub, and prepare next R script
source("SGGP_redTime_qsub_unevaluated.R")
