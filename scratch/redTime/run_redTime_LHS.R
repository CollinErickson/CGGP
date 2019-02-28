# Run this file on Crunch to run a bunch of redTime runs

# Original unix command
# Rscript -e "source('write_param_files.R'); source('write_sh_files.R'); create_LHS_and_submit(n=8039, prefix='LHS1L_n8039_s1226_', holdnum=80, seed=1226);"

# Source needed files
source('write_param_files.R')
source('write_sh_files.R')

# create_LHS_and_submit(n=8039, prefix='LHS1L_n8039_s1226_', holdnum=80, seed=1226)
# Running with expanded parameter ranges
# create_LHS_and_submit(n=1000, prefix='ExpandedRanges_LHS1L_n1000_s0227_', holdnum=40, seed=0227)
# create_LHS_and_submit(n=1000, prefix='ExpandedRanges2_LHS1L_n1000_s0227_', holdnum=40, seed=0227)
create_LHS_and_submit(n=100, prefix='ExpandedRanges2_LHS1L_n100_s0228_', holdnum=50, seed=0228)


# After done, run
# extract_redTime_from_completed_LHS(n=, prefix=)
# extract_redTime_from_completed_LHS(n=1000, prefix='ExpandedRanges2_LHS1L_n1000_s0227_')