

# Now need to evaluate SG$design_unevaluated

# Write param files and sh files, then qsub them with hold to avoid overload
source('write_param_files.R')
source('write_sh_files.R')
for (i in 1:nrow(SG$design_unevaluated)) {
  # Write param files
  write_params_file(x01=SG$design_unevaluated[i,],
                    fileID=paste0(fileIDbase, "_", SG$ss, "_", i))
  
  # Write sh files
  write_sh_file(fileID=paste0(fileIDbase, "_", SG$ss, "_", i))
  
  # qsub .sh files
  qsub_sh_file(fileID=paste0(fileIDbase, "_", SG$ss, "_", i), 
               holdID=paste0(fileIDbase, "_", SG$ss, "_", i - number_cores))
}

# qsub continue redTime (this file) so it will run after rest are finished
qsub_string <- "qsub "
qsub_string <- paste0(qsub_string, " -N SGGPmaster", SG$ss) # Identify it with size of design
qsub_string <- paste0(qsub_string, " -hold_jid ", prefix, i) # Hold until last one submitted finishes, might still be others running though
qsub_string <- paste(qsub_string, " SGGPmaster.sh")
# qsub_string <- paste(qsub_string, )
# system('qsub -N SGGPmaster -hold_jid paste0(prefix, i - number_cores) thisfile.sh')
system(qsub_string)