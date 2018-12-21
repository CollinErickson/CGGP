# Parameters for the SGGP/redTime run

groupID <- "redTimeTest1"

number_cores <- 10
d <- 8
N0 <- 400
corr <- "CauchySQT"

batchsize <- 400
Nfinal <- 4000

fileIDbase
SGGP_after_fit_RDS_path <- paste0("/home/collin/scratch/redTime_v0.1/output_files/out_", groupID,
                                  "_SGGP_after_fit.out")
SGGP_after_append_RDS_path <- paste0("/home/collin/scratch/redTime_v0.1/output_files/out_", groupID,
                                     "_SGGP_after_append.out")