# Parameters for the SGGP/redTime run

groupID <- "redTimeTest1"

number_cores <- 10
d <- 8
N0 <- 10
corr <- "CauchySQT"

batchsize <- 10
Nfinal <- 50

sourcefilepath <- "/home/collin/scratch/SGGP/scratch/redTime/"

parampathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/param_files/params_", groupID)
shpathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/sh_files/sh_", groupID)

outpathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/output_files/out_", groupID)

SGGP_after_fit_RDS_path    <- paste0(outpathbase, "_SGGP_after_fit.out")
SGGP_after_append_RDS_path <- paste0(outpathbase, "_SGGP_after_append.out")