# This will submit all the sh files to run the
# points that the SGGP object chooses to run.

# Now need to evaluate SG$design_unevaluated
# cat(getwd(), '\n')

# Write param files and sh files, then qsub them with hold to avoid overload
source(paste0(sourcefilepath, 'write_param_files.R'))
source(paste0(sourcefilepath, 'write_sh_files.R'))
for (i in 1:nrow(SG$design_unevaluated)) {
  # Write param files
  write_params_file(x01=SG$design_unevaluated[i,],
                    parampathbase=parampathbase,
                    fileID=paste0("_", SG$ss, "_", i))
  
  # Write sh files
  write_sh_file(fileID=paste0("_", SG$ss, "_", i),
                shpathbase=shpathbase,
                parampathbase=parampathbase,
                outpathbase=outpathbase#,
                # No longer having last one gather results, have separate master
                # additional_command=if (i == nrow(SG$design_unevaluated)) {paste0("Rscript ", sourcefilepath, "SGGP_redTime_continue.R")} else {""},
                # node= if (i == nrow(SG$design_unevaluated)) {"crunch.local"} else {NULL}
                )
  
  # qsub .sh files
  qsub_sh_file(fileID=paste0("_", SG$ss, "_", i), 
               holdID=paste0("_", SG$ss, "_", i - number_cores),
               shpathbase=shpathbase)
}

# Below is with a master for each iteration.
# I'm setting up a permanent master in SGGP_redTime_masterscript/sub
if (FALSE) {
  # Now set up master file to run after all of these finish
  source(paste0(sourcefilelocation, "write_sh_master_file.R"))
  write_sh_master_file(fileID=paste0("_", SG$ss, "_SGGPmaster"))
  shpath <- paste0(shpathbase, "_", SG$ss, "_SGGPmaster")
  system(paste("chmod +x ", shpath))
  
  ## qsub continue redTime (this file) so it will run after rest are finished
  qsub_string <- "qsub "
  qsub_string <- paste0(qsub_string, " -N SGGPmaster", SG$ss) # Identify it with size of design
  qsub_string <- paste0(qsub_string, " -hold_jid ", "_", SG$ss, "_", 1) # Hold until last one submitted finishes, might still be others running though
  for (i in 2:nrow(SG$design_unevaluated)) {
    qsub_string <- paste0(qsub_string, ",_",SG$ss,"_",i)
  }
  qsub_string <- paste(qsub_string, " SGGPmaster.sh")
  # qsub_string <- paste(qsub_string, )
  # system('qsub -N SGGPmaster -hold_jid paste0(prefix, i - number_cores) thisfile.sh')
  cat("qsub_string is ", qsub_string, " \n")
  system(qsub_string)
  cat("qsub_string submitted\n")
}