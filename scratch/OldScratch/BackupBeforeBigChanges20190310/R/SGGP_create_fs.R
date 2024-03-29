#' Create sparse grid GP
#'
#' @param d Input dimension
# @param xmin Min x values, vector. Must be rep(0,d).
# @param xmax Max x values, vector. Must be rep(1,d).
#' @param batchsize Number added to design each batch
# @param nugget Nugget term added to diagonal of correlation matrix,
#' for now only on predictions
#' @param corr Name of correlation function to use. Must be one of "CauchySQT",
#' "CauchySQ", "Cauchy", "Gaussian", "PowerExp", "Matern32", "Matern52".
#' @param grid_sizes Size of grid refinements.
#' @param Xs Supplemental X data
#' @param Ys Supplemental Y data
#' @param supp_args Arguments used to fit if Xs and Ys are given
#' @param ... Force named arguments
#' @param HandlingSuppData How should supplementary data be handled?
#' * Correct: full likelihood with grid and supplemental data
#' * Only: only use supplemental data
#' * Ignore: ignore supplemental data
#' * Mixture: sum of grid LLH and supplemental LLH, not statistically valid
#' * MarginalValidation: a validation shortcut
#' * FullValidation: a validation shortcut
#' 
#' @importFrom stats rbeta
#'
#' @return SGGP
#' @export
#' @family SGGP core functions
#'
#' @examples
#' d <- 8
#' SG = SGGPcreate(d,200)
SGGPcreate <- function(d, batchsize, corr="CauchySQ",
                       ...,
                       grid_sizes=c(1,2,4,4,8,12,32),
                       Xs=NULL, Ys=NULL,
                       HandlingSuppData=if (is.null(Xs)) {"Ignore"} else {"Correct"},
                       supp_args=list()
) {
  if (d <= 1) {stop("d must be at least 2")}
  if (length(list(...))>0) {stop("Unnamed arguments given to SGGPcreate")}
  # This is list representing our GP object
  SGGP = list()
  class(SGGP) <- c("SGGP", "list") # Give it class SGGP
  
  SGGP$d <- d
  SGGP$HandlingSuppData <- HandlingSuppData
  SGGP <- SGGP_internal_set_corr(SGGP, corr)
  
  # Partial matching is very bad! Keep these as length 0 instead of NULL,
  #  otherwise SGGP$Y can return SGGP$Ys
  SGGP$Y <- numeric(0)
  SGGP$y <- numeric(0)
  
  # If supplemental data is given, fit it here
  if (!is.null(Xs) && !is.null(Ys)) {
    if (!is.null(supp_args) && length(supp_args) > 0 && is.null(names(supp_args))) {stop("Give names for supp_args")}
    supp_args$SGGP <- SGGP
    supp_args$Xs <- Xs
    supp_args$Ys <- Ys
    SGGP <- do.call(SGGP_internal_fitwithonlysupp, supp_args)
  }
  
  # Levels are blocks. Level is like eta from paper.
  SGGP$ML = min(choose(SGGP$d + 6, SGGP$d), 10000) #max levels
  
  # Track evaluated blocks, aka used levels
  SGGP$uo = matrix(0, nrow = SGGP$ML, ncol = SGGP$d) # blocks that have been selected
  SGGP$uoCOUNT = 0 # number of selected blocks
  
  # Track the blocks that are allowed to be evaluated
  SGGP$po = matrix(0, nrow = 4 * SGGP$ML, ncol = SGGP$d) #proposed levels tracker
  # Only option at first is initial block (1,1,...,1)
  SGGP$po[1, ] <- rep(1, SGGP$d)
  SGGP$poCOUNT <- 1
  
  # Ancestors are blocks one level down in any dimension.
  SGGP$maxgridsize = 400
  SGGP$pila = matrix(0, nrow = SGGP$ML, ncol =SGGP$maxgridsize ) #proposed immediate level ancestors
  SGGP$pala = matrix(0, nrow = SGGP$ML, ncol =SGGP$maxgridsize ) #proposedal all level ancestors
  SGGP$uala = matrix(0, nrow = SGGP$ML, ncol =SGGP$maxgridsize ) #used all level ancestors
  SGGP$pilaCOUNT = rep(0, SGGP$ML) #count of number of pila
  SGGP$palaCOUNT = rep(0, SGGP$ML) #count of number of pala
  SGGP$ualaCOUNT = rep(0, SGGP$ML) #count of number of uala
  
  # Initial block (1,1,...,1) has no ancestors
  SGGP$pilaCOUNT[1] <- 0
  SGGP$pila[1, 1] <- 0
  
  # SGGP$sizes = c(1,2,4,4,8,12,32) # Num of points added to 1D design as you go further in any dimension
  SGGP$sizes <- grid_sizes
  SGGP$maxlevel = length(SGGP$sizes)
  # Proposed grid size
  SGGP$pogsize = rep(0, 4 * SGGP$ML)
  SGGP$pogsize[1:SGGP$poCOUNT] = apply(matrix(SGGP$sizes[SGGP$po[1:SGGP$poCOUNT, ]], SGGP$poCOUNT, SGGP$d), 1, prod)
  # Selected sample size
  SGGP$ss = 0
  
  
  SGGP$w = rep(0, SGGP$ML) #keep track of + and - for prediction
  SGGP$uoCOUNT = 0 ###1 # Number of used levels
  # While number selected + min sample size <= batch size, i.e., still have enough spots for a block
  while (batchsize > (SGGP$ss + min(SGGP$pogsize[1:SGGP$poCOUNT]) - 0.5)) {
    SGGP$uoCOUNT = SGGP$uoCOUNT + 1 #increment used count
    
    if (SGGP$uoCOUNT < 1.5) { # Nothing picked yet, so take base block (1,1,...,1)
      pstar <- 1
    } else if (SGGP$uoCOUNT < (SGGP$d + 1.5)) {
      # Next d iterations pick the (2,1,1,1,1),(1,2,1,1,1) blocks b/c we need info on each dimension before going adaptive
      pstar = 1
    } else{ # The next d iterations randomly pick from the boxes with minimal number of points
      if (SGGP$uoCOUNT < (2 * SGGP$d + 1.5)) {
        pstar = sample(which(SGGP$pogsize[1:SGGP$poCOUNT] <= 0.5 + min(SGGP$pogsize[1:SGGP$poCOUNT])), 1)
      } else{ # After that randomly select from blocks that still fit
        pstar = sample(which(SGGP$pogsize[1:SGGP$poCOUNT] < min(batchsize - SGGP$ss + 0.5,SGGP$maxgridsize)), 1)
      }
    }
    
    l0 =  SGGP$po[pstar, ] # Selected block e.g. (2,1,1,2)
    # Need to make sure there is still an open row in uo to set with new values
    if (SGGP$uoCOUNT > nrow(SGGP$uo)) {
      SGGP <- SGGP_internal_addrows(SGGP)
    }
    SGGP$uo[SGGP$uoCOUNT, ] = l0 # Store new block
    SGGP$ss =  SGGP$ss + SGGP$pogsize[pstar] # Update selected sample size
    
    # Ancestors of block just selected
    # Need to give possibility for initial block, has no ancestors, and 1:0 is bad
    new_an = SGGP$pila[pstar, if (SGGP$pilaCOUNT[pstar]>.5) {1:SGGP$pilaCOUNT[pstar]} else {numeric(0)}]
    total_an = new_an
    
    # Loop over ancestors of block just selected
    if (length(total_an) > .5) { # Initial block has no total_an , need this to avoid 1:0
      for (anlcv in 1:length(total_an)) {
        # If more than one ancestor, update with unique ones.
        if (total_an[anlcv] > 1.5) {
          total_an = unique(c(total_an, SGGP$uala[total_an[anlcv], 1:SGGP$ualaCOUNT[total_an[anlcv]]]))
        }
      }
      # Update storage of ancestors
      SGGP$ualaCOUNT[SGGP$uoCOUNT]  = length(total_an)
      SGGP$uala[SGGP$uoCOUNT, 1:length(total_an)] = total_an
    }
    
    # Loop over ancestors, update coefficient
    if (length(total_an) > .5) { # Initial block has no total_an , need this to avoid 1:0
      for (anlcv in 1:length(total_an)) {
        lo = SGGP$uo[total_an[anlcv], ]
        if (max(abs(lo - l0)) < 1.5) {
          SGGP$w[total_an[anlcv]] = SGGP$w[total_an[anlcv]] + (-1) ^ abs(round(sum(l0-lo)))
          
        }
      }
    }
    SGGP$w[SGGP$uoCOUNT] = SGGP$w[SGGP$uoCOUNT] + 1
    
    SGGP$po[pstar,] <- 0 # Clear just used row
    if (SGGP$poCOUNT > 1.5) { # Move up other options if there are others left
      po_rows_to_move <- (1:SGGP$poCOUNT)[-pstar] # Moving all po rows except selected
      SGGP$po[1:(SGGP$poCOUNT - 1), ] = SGGP$po[po_rows_to_move, ]
      SGGP$pila[1:(SGGP$poCOUNT - 1), ] = SGGP$pila[po_rows_to_move, ]
      SGGP$pilaCOUNT[1:(SGGP$poCOUNT - 1)] = SGGP$pilaCOUNT[po_rows_to_move]
      SGGP$pogsize[1:(SGGP$poCOUNT - 1)] = SGGP$pogsize[po_rows_to_move]
    }
    
    # One less option now
    SGGP$poCOUNT = SGGP$poCOUNT - 1
    
    # Loop over dimensions to add new possible blocks
    for (dimlcv in 1:SGGP$d) {
      # The block e.g. (1,2,1,1,3) just selected
      lp = l0
      
      lp[dimlcv] = lp[dimlcv] + 1 # Increase single dimension by 1, will see if it is possible
      
      # Check if within some bounds??
      if (max(lp) <= SGGP$maxlevel && SGGP$poCOUNT < 4*SGGP$ML) {
        # Dimensions which are past first design level
        kvals = which(lp > 1.5)
        
        canuse = 1 # Can this block be used? Will be set to 0 below if not.
        ap = rep(0, SGGP$d) # ????
        nap = 0 # ?????
        
        # Loop over dims at 2+ and do what?
        for (activedimlcv in 1:length(kvals)) {
          lpp = lp # The block selected with 1 dim incremented
          lpp[kvals[activedimlcv]] = lpp[kvals[activedimlcv]] - 1 # ????
          
          ismem = rep(1, SGGP$uoCOUNT) # Boolean???
          # Loop over dimensions
          for (dimdimlcv in 1:SGGP$d) { # Set to 0 or 1 if all points already selected have same value???????
            ismem  = ismem * (SGGP$uo[1:SGGP$uoCOUNT, dimdimlcv] == lpp[dimdimlcv])
          }
          # If any are still 1,
          if (max(ismem) > 0.5) {
            ap[activedimlcv] = which(ismem > 0.5)
            nap = nap + 1 # Count number that are >=1
          } else{ # All are 0, so can't use
            canuse = 0
          }
        }
        # If it can be used, add to possible blocks
        if (canuse > 0.5) {
          SGGP$poCOUNT = SGGP$poCOUNT + 1
          SGGP$po[SGGP$poCOUNT, ] = lp
          SGGP$pogsize[SGGP$poCOUNT] = prod(SGGP$sizes[lp])
          SGGP$pila[SGGP$poCOUNT, 1:nap] = ap[1:nap]
          SGGP$pilaCOUNT[SGGP$poCOUNT] = nap
          
        }
      }
    }
  }
  
  # Create points for design
  #  These are distances from the center 0.5.
  xb = rep(
    c(
      1 / 2, # 0, 1
      3 / 8, # 1/8, 7/8
      1 / 4, # 1/4, 3/4
      1 / 8,
      15 / 32,
      7 / 16,
      3 / 16,
      5 / 16,
      7 / 32,
      11 / 32,
      3 / 32,
      13 / 32,
      9 / 32,
      5 / 32,
      1 / 32,
      1 / 16,
      seq(31,1,-2)/64
    ),
    "each" = 2
  )
  SGGP$xb = 0.5 + c(0, xb * rep(c(-1, 1), length(xb) / 2))
  SGGP$xindex = 1:length(xb) # Why not length(SGGP$xb), which is one longer than xb?
  # After this xb is
  #  [1] 0.50000 0.12500 0.87500 0.25000 0.75000 0.37500 0.62500 0.28125 0.71875 0.31250 0.68750 0.00000 1.00000 0.18750 0.81250
  # [16] 0.06250 0.93750 0.43750 0.56250 0.40625 0.59375 0.09375 0.90625 0.21875 0.78125 0.34375 0.65625 0.46875 0.53125 0.15625
  # [31] 0.84375 0.03125 0.96875
  SGGP$sizest = cumsum(SGGP$sizes) # Total number of points in 1D design as you go along axis
  
  
  # This is all to create design from uo.
  # If only supp data is given, don't run it.
  if (SGGP$uoCOUNT > 0) {
    
    # Get design from uo and other data
    SGGP <- SGGP_internal_getdesignfromSGGP(SGGP)
    
    SGGP$design_unevaluated <- SGGP$design
  }
  return(SGGP)
}