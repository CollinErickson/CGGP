# Trying to figure out how I'll run a bunch of jobs
# 1. Get design points (e.g. from LH)
# 2. Write out param files
# 3. Write out sh/pbs files
# 4. qsub the files

create_LHS_and_submit <- function(n, prefix='', holdnum=NULL, seed=NULL) {
  if (!is.null(seed)) {set.seed(seed)}
  X <- lhs::maximinLHS(n=n, k=8)
  sapply(1:n,
         function(i) {
           write_params_file(x01 = X[i,], fileID = paste0(prefix,i))
           write_sh_file(fileID=paste0(prefix,i))
           if (is.null(holdnum)) {
             qsub_sh_file(fileID=paste0(prefix,i))
           } else {
             if (as.integer(holdnum) != holdnum) {stop("holdnum isn't integer")}
             qsub_sh_file(fileID=paste0(prefix,i), holdID=paste0(prefix, i-holdnum))
           }
         }
  )
}

write_sh_file <- function(fileID,
                          parampathbase="/home/collin/scratch/redTime_v0.1/sub_files/params_redTime_",
                          shpathbase="/home/collin/scratch/redTime_v0.1/sub_files/sub_",
                          outpathbase="/home/collin/scratch/redTime_v0.1/output_files/out_",
                          additional_command="",
                          email="a"
                          ) {
  param_path <- paste0(parampathbase, fileID, ".dat")
  shpath <- paste0(shpathbase, fileID, ".sh")
  if (file.exists(shpath)) {stop(paste("File already exists", shpath))}
  outputpath <- paste0(outpathbase, fileID, ".out")
  if (file.exists(outputpath)) {stop(paste("File already exists", outputpath))}
  #shpath <- ""
  cout <- function(...) {cat(..., '\n', file=shpath)}
  cout("
#!/bin/bash
#

########################
##### qsub options #####
########################
#
# Any line starting with #$ is interpreted as command line option for qsub

#### Working directory
#  These options specify the directory where the job is to be executed.
#  The default is your home directory (which you don't want).
#  If you use the -wd option, you need to make sure that directory exists.
#
#$ -cwd                           # Current working directory (where you run qsub)
#                                 # MUST BE UNDER $HOME/scratch !!!!!
#
# #$ -wd /home/me/scratch/CoolStuffIsHere # Run job /home/me/scratch/CoolStuffIsHere

##### Shell that is used
#$ -S /bin/bash

#### Number of threads
#  If you have a multi-threaded application, you need to specify here how many
#  cores your process uses.
#  Note: You explicitly have to tell you program how many threads to use
#$ -pe smp 1

#### Run time limit
#  Specify maximum CPU time after which job is to be killed (format HH:MM:SS).
#$ -l h_rt=4:00:00    # in this example, we set 10 minutes

#### Memory limit
#  specifies the maximum amount of memory this job can take
#  This is per thread, so the total amount is this number times the number
#  of threads. The default value is 2g.
#$ -l h_vmem=2g  # here we choose 4g, so that overall we reserve up to 8*4g=32g total

#### Email after done, -abe is abort, begin, end
##### $ -m ae no longer want emails, running too many
##### $ -M collinerickson@u.northwestern.edu

#########################
##### Your commands #####
#########################

cd /home/collin/scratch/redTime_v0.1
date
")
       cout(paste("date >> ", outputpath),append=T)
       cout(paste("./redTime.out ", param_path, " >> ", outputpath), append=T)
       cout(paste("date >> ", outputpath),append=T)
       cout("
date
       ", append=T)
       cout(additional_command, append=T)
}

qsub_sh_file <- function(fileID, holdID=NULL, shpathbase="/home/collin/scratch/redTime_v0.1/sub_files/sub_") {
  
  shpath <- paste0(shpathbase, fileID, ".sh")
  system(paste("chmod +x ", shpath))
  print(paste("About to qsub", shpath))
  qsub_string <- paste("qsub -N ", fileID) # Give it a name
  if (!is.null(holdID)) {
    qsub_string <- paste(qsub_string, "-hold_jid", holdID)
  }
  qsub_string <- paste(qsub_string, shpath) # shpath to run
  system(qsub_string)
  print(paste("qsubbed ", shpath))
}
