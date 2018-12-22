write_sh_master_file <- function(fileID,
                          # parampathbase="/home/collin/scratch/redTime_v0.1/sub_files/params_redTime_",
                          shpathbase="/home/collin/scratch/redTime_v0.1/sub_files/sub_",
                          outpathbase="/home/collin/scratch/redTime_v0.1/output_files/out_",
                          # additional_command="",
                          sourcefilepath,
                          email="a",
                          node=NULL
) {
  # param_path <- paste0(parampathbase, fileID, ".dat")
  shpath <- paste0(shpathbase, fileID, ".sh")
  if (file.exists(shpath)) {stop(paste("File already exists", shpath))}
  outputpath <- paste0(outpathbase, fileID, ".out")
  if (file.exists(outputpath)) {stop(paste("File already exists", outputpath))}
  # shpath <- paste0(shpathbase, "_SGGPmaster.sh")
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
       ")
  if (!is.null(node)) {
    cout("#$ -l h=crunch.local", append=TRUE)
  }
  cout("
       
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
       ", append=TRUE)
  cout(paste("date >> ", outputpath),append=T)
  cout(paste("pwd >> ", outputpath), append=T)
  # cout("source ~/.bashrc", append=T)
  # cout(paste("./redTime.out ", param_path, " >> ", outputpath), append=T)
  cout(paste0("Rscript ", sourcefilepath, "SGGP_redTime_continue.R"), append=T)
  cout(paste("date >> ", outputpath),append=T)
  cout("
       date
       ", append=T)
  cout(additional_command, append=T)
}
