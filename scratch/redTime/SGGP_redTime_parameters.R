# Parameters for the SGGP/redTime run

# T1 was first test, used PCA and shared parameters. Y not on log scale.
# T2 was also PCA and shared parameters b/c I forgot to change it. Y on log scale.
# T3 was separate opd and no PCA.
# T4 was same as T3 but with adding sigma2hat in append, so it should focus more on bad spots
# Sup1 is with supplementary data
# Sup2 is with supp data, but on expandedranges2
# Big1 is all output dimensions, no pca, shared params. Meant to be for paper. Not good b/c of UCB error
# Big2 is same as Big1 except using Greedy in append
# Big3 is back same as Big1, i.e. uses UCB, we think we fixed the error.
# ER3a is using ExpandedRanges3, Greedy, Power Exp, 
groupID <- "redTimeTestER3a"
groupID_short <- "ER3a"

# Number of cores to use at a time. Keep <= 40 so others can use server.
number_cores <- 125
hold_in_groups <- TRUE

# Input dimensions
d <- 9

# Initial sample size
N0 <- 200
# Number of points to add in each batch
# batchsize <- 100 #100
batchsize1 <- 200
batchsize2 <- 500
batchsize3 <- 2000
batchsize4 <- 10000

# Correlation function
corr <- "PowerExp" # "CauchySQ"

# Number of points after which to stop (will go up to batchsize-1 over)
Nfinal <- 25000

# Grid size to use. This option wasn't included in Test1
# grid_size <- c(1, 2, 2, 2, 4, 4, 4, 4, 4, 6, 32)
# grid_sizes <- c(1, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 8, 12)
# grid_sizes <- c(1,2,4,4,8,12,32) # This is default, was "fast" in our internal comparisons. Used for Bigx
grid_sizes <- c(1,2,4,4,8,12,20,28,32) # new default, goes up to 111. Used for ER3.

# PCA no longer an option
# # use_PCA, 100 outputs, PCA can reduce to 37 I think.
# use_PCA <- FALSE

# Much slower to fit parameters separately
separateoutputparameterdimensions <- FALSE

# outdims: output dimensions to use. Default should be 1:100.
outdims <- 1:80 # 1:100 is all, 1:76 is up to .2, 1:80 is up to .25

# append selectionmethod. Used to just use UCB, but UCB/TS are bad with MV out. Use Greedy instead.
selectionmethod <- "Greedy"

# To use supplementary data
if (TRUE) {
  # stop("use supp")
  Xsup <- NULL
  Ysup <- NULL
} else {
  stop("Not using supp any more")
  # Xall <- unname(as.matrix(read.csv("/home/collin/scratch/SGGP/scratch/redTime/data/LHS1L_n8039_s1226_Xmatrix.csv")[,-1]))
  # Yall <- log(unname(as.matrix(read.csv("/home/collin/scratch/SGGP/scratch/redTime/data/LHS1L_n8039_s1226_all_output.csv")[,-1])))
  # Xall <-     unname(as.matrix(read.csv("/home/collin/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n90_s0315_all_input.csv")[,-1]))
  # Yall <- log(unname(as.matrix(read.csv("/home/collin/scratch/redTime_v0.1/SGGPruns/important_files/ExpandedRanges2_LHS1L_n90_s0315_all_output.csv")[,-1])))
  # set.seed(100) # Set seed for reproducibility
  # SupRows <- sample(1:nrow(Xall), 100, replace=FALSE)
  SupRows <- 1:nrow(Xall)
  # set.seed(Sys.time()) # Clear seed in case running a repeat
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
save_after <- c(100, 200, 300, 500, 700, 900, 1100,1500,1900,2400,3000,4000,5000,6000,7000,8000,10000,
                12000,15000,18000,21000,23000,25000,30000,35000,40000)
# save_after <- c(50,150,250,350,450,550)

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
