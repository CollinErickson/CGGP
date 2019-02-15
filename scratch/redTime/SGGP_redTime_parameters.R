# Parameters for the SGGP/redTime run

# T1 was first test, used PCA and shared parameters. Y not on log scale.
# T2 was also PCA and shared parameters b/c I forgot to change it. Y on log scale.
# T3 was separate opd and no PCA.
# T4 was same as T3 but with adding sigma2hat in append, so it should focus more on bad spots
# Sup1 is with supplementary data
groupID <- "redTimeTestO50_1k"
groupID_short <- "O50_1k"

# Number of cores to use at a time. Keep <= 40 so others can use server.
number_cores <- 50

# Input dimensions
d <- 9

# Initial sample size
N0 <- 200
# Number of points to add in each batch
batchsize <- 200

# Correlation function
corr <- "CauchySQ"

# Number of points after which to stop (will go up to batchsize-1 over)
Nfinal <- 1100

# Grid size to use. This option wasn't included in Test1
# grid_size <- c(1, 2, 2, 2, 4, 4, 4, 4, 4, 6, 32)
grid_size <- c(1, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 8, 12)

# use_PCA, 100 outputs, PCA can reduce to 37 I think.
use_PCA <- FALSE

# Much slower to fit parameters separately
separateoutputparameterdimensions <- TRUE

# outdims: output dimensions to use. Default should be 1:100.
outdims <- 50

# To use supplementary data
if (TRUE) {
  Xsup <- NULL
  Ysup <- NULL
} else {
  Xall <- unname(as.matrix(read.csv("/home/collin/scratch/SGGP/scratch/redTime/data/LHS1L_n8039_s1226_Xmatrix.csv")[,-1]))
  Yall <- log(unname(as.matrix(read.csv("/home/collin/scratch/SGGP/scratch/redTime/data/LHS1L_n8039_s1226_all_output.csv")[,-1])))
  set.seed(100) # Set seed for reproducibility
  SupRows <- sample(1:nrow(Xall), 100, replace=FALSE)
  set.seed(Sys.time()) # Clear seed in case running a repeat
  Xsup <- Xall[SupRows,]
  Ysup <- Yall[SupRows,]
  cat("Str of Xsup is\n", str(Xsup), "\n")
  cat("Str of Ysup is\n", str(Ysup), "\n")
  rm(Xall, Yall)
}

# Should log of redTime output be used? Test1 didn't, Test2 and on will
use_log_redTime <- TRUE

# When should the object be saved
# save_after <- c(200, 400, 1000, 2000, 4000, 8000)
save_after <- c(100, 300, 500, 700, 900, 1100)

sourcefilepath <- "/home/collin/scratch/SGGP/scratch/redTime/"

dir.create(paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID))
dir.create(paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/param_files"))

dir.create(paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/sh_files"))
dir.create(paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/output_files"))



parampathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/param_files/params_") #, groupID)
shpathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/sh_files/sh_") #, groupID)

outpathbase <- paste0("/home/collin/scratch/redTime_v0.1/SGGPruns/", groupID, "/output_files/out_")# , groupID)

SGGP_after_fit_RDS_path    <- paste0(outpathbase, groupID_short, "_SGGP_after_fit.rds")
SGGP_save_after_fit_RDS_path    <- paste0(outpathbase, groupID_short, "_SGGP") # For saving intermediate objects
SGGP_after_append_RDS_path <- paste0(outpathbase, groupID_short, "_SGGP_after_append.rds")
