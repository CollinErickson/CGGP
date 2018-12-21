# Write out params_redTimexx.out file to run on cluster

convert_x_from_01_to_ranges <- function(x,
                                        low= c(.85 ,.7,.55,.12,.0215,0,-1.3,-1.5),
                                        high=c(1.05,.9,.85,.155,.0235,.01,-.7,1.15),
                                        low_redshift=0, high_redshift=5
                                        ) {
  if (any(x<0) || any(x>1)) {stop("x must be in range [0,1]^n")}
  if (length(x)==9) {low <- c(low, low_redshift); high <- c(high, high_redshift)}
  else if (length(x)==8) {}
  else {stop(paste("x in convert_x_from_01_to_ranges is bad length", length(x), "values are", x))}
  low + (high - low) * x
}

write_params_file <- function(..., x01, fileID, overwrite=F,
                              parampathbase="/home/collin/scratch/redTime_v0.1/sub_files/params_redTime_") {
  x <- convert_x_from_01_to_ranges(x01)
  n_s <- x[1]
  sigma_8<- x[2]
  h<- x[3]
  Omega_m<- x[4]
  Omega_b<- x[5]
  Omega_nu<- x[6]
  w0<- x[7]
  wa<- x[8]
  redshift <- x[9] # This is new, used to be set for whole system
  paste0(rep('0',1),as.character(3))
  outpath <- paste0(parampathbase, fileID, ".dat")
  if (outpath == "/home/collin/scratch/redTime_v0.1/params_redTime.dat") {
    stop("Pick a new outpath, don't overwrite params_redTime.dat")
  }
  if (file.exists(outpath)) {stop(paste("File already exists", outpath))}
  #outpath <- ""
  cout <- function(..., append=T) {cat(..., '\n', file=outpath, append=append)}
  cout("# ==============================================================================
# =================== Parameter file for Time-RG computation ===================
# ==============================================================================
#
# (Lines beginning with a '#' are ignored by the code.  Parameters must be
# listed in the following order:
#   n_s
#   sigma_8
#   h
#   Omega_m
#   Omega_b
#   Omega_nu
#   T_cmb[Kelvins]
#   w_0
#   w_a
#   switch_nonlinear
#   switch_1loop
#   switch_print_linear
#   switch_print_rsd
#   z_initial
#   num_z_outputs
#   z_out_1 z_out_2 ... z_out_{num_z_outputs}
#   file_transfer_function
#   num_massive_nu_approx
#   file_nu_transfer_root
#   num_interp_redshifts
#   z_1 z_2 ... z_num_interp_redshifts
# where parameters beginning with \"switch\" are 0 or 1, those beginning with
# \"num\" are integers, and those beginning with \"file\" are character strings.)
#
# -------------------------- cosmologial parameters ----------------------------
#", append=T)
  cout("# n_s: scalar spectral index", append=T)
  cout(n_s)
  cout("#")
  cout("# sigma_8: z=0 normalization of linear power spectrum")
  cout(sigma_8)
  cout("#")
  cout("# h: Hubble parameter today, H_0 / (100 km/sec/Mpc)")
  cout(h)
  cout("#")
  cout("# Omega_m: total matter fraction today (cdm, baryons, massive nu)")
  cout(Omega_m)
  cout("#")
  cout("# Omega_b: baryon fraction today")
  cout(Omega_b)
  cout("#")
  cout("# Omega_nu: massive neutrino fraction today (0 if massless)")
  cout(Omega_nu)
  cout("#")
  cout("# T_cmb_K: CMB temperature today, in units of Kelvins")
  cout(2.726) #T_cmb_K)
  cout("#")
  cout("# w0: dark energy equation of state today")
  cout(w0)
  cout("#")
  cout("# wa: derivative -dw/da in CPL parameterization")
  cout(wa)
  cout("#")
  cout("# ----------------------------- code switches ----------------------------------
#
# nonlinear computation: 0 for linear, 1 for nonlinear
1
#
# 1-loop approximation to nonlinear computation (0 for full NL, 1 for 1-loop)
# 0 # Default is this, the full calc, but takes 8 hours.
1
#
# print linear evolution information: D, f, P_lin for cb and nu
1
#
# print RSD information: P_{B,j} {j=2,4,6} and P_{T,j} {j=2,4,6,8}
1
#
# -------------------------------- outputs -------------------------------------
#
# initial redshift
200
#
# number of redshifts at which to print results
# need to match number of redshifts given below 7
1
#
# redshifts of outputs (arranged from greatest to least)
# trying single redshift below: 4.95 4 3.04 2.02 1.006 0.511 0
# 2.02
# Now setting redshift as input")
  cout(redshift)
  cout("#
# ---------------------------- transfer inputs ---------------------------------
#
# Transfer function at z=0, in CAMB standard format (7 columns:
#   k[h/Mpc], delta_cdm/k^2, delta_b/k^2, delta_gamma/k^2, delta_massless/k^2,
#   delta_massive/k^2, delta_tot/k^2
# where massive and massless refer to neutrinos.
camb_transfer_z0.dat
#
# Massive neutrino approx: 0 for CAMB interpolation, which is the only method
# implemented in the code so far.
0
#
# If using CAMB interpolation, provide the transfer filename root, the number
# of redshifts, and the list of redshifts from greatest to least.  A transfer
# file for a given redshift must be named {transfer root}{redshift}.dat.  For
# example, the transfer file corresponding to a transfer root of
# camb_transfer_z and a redshift of 3 is camb_transfer_z3.dat.
# transfer root:
camb_transfer_z
#
# Number of interpolation redshifts:
12
#
# List of redshifts, from greatest to least, separated by whitespace.  These
# should span the range z_initial <= z <= 0.
200 100 50 20 10 5 4 3 2 1 .5 0
#
# ------------------------------------------------------------------------------")

  # HAVING ISSUE with spaces at ends of lines, aborts if spaces on end of numeric lines.
  # Use sed to remove trailing spaces
  system(paste("sed -i 's/\\s*$//' ", outpath))
}
