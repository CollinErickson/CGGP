# Parameters for the SGGP/redTime run

groupID <- "redTimeTest1"
groupID_short <- "T1"

# Number of cores to use at a time. Keep <= 40 so others can use server.
number_cores <- 40

# Input dimensions
d <- 9

# Initial sample size
N0 <- 40

# Correlation function
corr <- "CauchySQT"

# Number of points to add in each batch
batchsize <- 40

# Number of points after which to stop (will go up to batchsize-1 over)
Nfinal <- 8000

sourcefilepath <- "/home/collin/scratch/SGGP/scratch/redTime/"

dir.create(paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID))
dir.create(paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/param_files"))

dir.create(paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/sh_files"))
dir.create(paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/output_files"))



parampathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/param_files/params_", groupID)
shpathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/sh_files/sh_", groupID)

outpathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/output_files/out_", groupID)

SGGP_after_fit_RDS_path    <- paste0(outpathbase, "_SGGP_after_fit.rds")
SGGP_after_append_RDS_path <- paste0(outpathbase, "_SGGP_after_append.rds")
