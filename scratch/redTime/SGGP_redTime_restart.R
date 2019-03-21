
restart_on <- 26155 # 26155
run_redTime <- 501:2000 # 501:2000
cat("Restarting redTime with SGGP from ", restart_on, "\n")



# This file will reload number_cores, d, SGGP_RDS_path, parameters, etc
source("SGGP_redTime_parameters.R")
cat("Loading SGGP_redTime_parameters.R successfully\n")

cat("Assuming redTime param files were created but not all run\n")
cat("qsubbing those set as run_redTime")
# source(paste0(sourcefilepath, 'write_param_files.R'))
source(paste0(sourcefilepath, 'write_sh_files.R'))

for (i in run_redTime) {
  # qsub .sh files
  if (hold_in_groups) {
    # This will hold until all of last number_cores groups finishes.
    # Makes it easier for others to get onto Crunch.
    qsub_sh_file(fileID=paste0(groupID_short, "_", restart_on, "_", i), 
                 holdID=paste0(groupID_short, "_", restart_on, "_",
                               ceiling(i/number_cores)*number_cores - number_cores + (1:number_cores),
                               collapse = ','),
                 shpathbase=shpathbase)
  } else {
    qsub_sh_file(fileID=paste0(groupID_short, "_", restart_on, "_", i), 
                 holdID=paste0(groupID_short, "_", restart_on, "_", i - number_cores),
                 shpathbase=shpathbase)
  }
}

cat("Finished qsubbing all the previous files.\n")


devtools::build()
cat("devtools::build() ran successfully\n")

install.packages(repos=NULL, pkgs="/home/collin/scratch/CGGP_1.0.tar.gz")
cat("SGGP installed correctly\n")

library('CGGP')
cat("CGGP loaded successfully\n")


# Start master job on crunch.local
system(paste("chmod 777 ", "SGGP_redTime_mastersub.sh"))
system('qsub -N SGGPmaster -l h=crunch SGGP_redTime_mastersub.sh')
cat("mastersub.sh submitted successfully\n")
cat(timestamp(), '\n')
