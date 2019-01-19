#' Create sparse grid GP
#'
#' @param d Input dimension
# @param xmin Min x values, vector. Must be rep(0,d).
# @param xmax Max x values, vector. Must be rep(1,d).
#' @param batchsize Number added to design each batch
# @param nugget Nugget term added to diagonal of correlation matrix,
#' for now only on predictions
#' @param corr Name of correlation function to use. Must be one of "CauchySQT", "CauchySQ".
#' @importFrom stats rbeta
#'
#' @return SGGP
#' @export
#'
#' @examples
#' d <- 8
#' SG = SGGPcreate(d,200)
SGGPcreate <- function(d, batchsize, corr="CauchySQT") {
  if (d <= 1) {stop("d must be at least 2")}
  # This is list representing our GP object
  SGGP = list()
  class(SGGP) <- c("SGGP", "list") # Give it class SGGP
  
  SGGP$d <- d
  # SGGP$CorrMat <- SGGP_internal_CorrMatCauchySQT
  if (tolower(corr) %in% c("cauchysqt")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatCauchySQT
  } else if (tolower(corr) %in% c("cauchysq")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatCauchySQ
  } else {
    stop("corr given to SGGPcreate should be one of CauchySQT or CauchySQ")
  }
  SGGP$numpara <- SGGP$CorrMat(return_numpara=TRUE)
  # print("Fix thetaMAP and thetaPostSamples in create")
  SGGP$thetaMAP <- rep(0,d*SGGP$numpara)
  SGGP$numPostSamples <- 100
  SGGP$thetaPostSamples  <- matrix(2*rbeta(d*SGGP$numpara*SGGP$numPostSamples , 0.5, 0.5)-1,ncol=SGGP$numPostSamples )
  SGGP$y <- rep(0,0)
  
  # Levels are blocks. Level is like eta from paper.
  SGGP$ML = min(choose(SGGP$d + 6, SGGP$d), 10000) #max levels
  
  # What is levelpoint? The current point? This is not used again in this file!
  # Not used ever, so I am removing it.
  # SGGP$levelpoint = rep(0, SGGP$ML)
  
  # Track evaluated blocks, aka used levels
  SGGP$uo = matrix(0, nrow = SGGP$ML, ncol = SGGP$d) #used levels tracker
  # First observation is always (1,1,...,1)
  SGGP$uo[1, ] = rep(1, SGGP$d) #first observation in middle of space
  SGGP$uoCOUNT = 1 #number of used levels
  
  # Track the blocks that are allowed to be evaluated
  SGGP$po = matrix(0, nrow = 4 * SGGP$ML, ncol = SGGP$d) #proposed levels tracker
  # Now possible blocks are (2,1,1,1,1), (1,2,1,1,1), (1,1,2,1,1), etc
  SGGP$po[1:SGGP$d, ] = matrix(1, nrow = SGGP$d, ncol = SGGP$d) + diag(SGGP$d) #one at a time
  SGGP$poCOUNT = SGGP$d #number of proposed levels
  
  
  # What are ancestors? Why are we doing this?
  # Are ancestors support blocks? Is this for calculating coefficient?
  # How any of its ancestors be proposed? They should all be used already?
  # Might be transposed??? Is this right?
  SGGP$maxgridsize = 400
  SGGP$pila = matrix(0, nrow = SGGP$ML, ncol =SGGP$maxgridsize ) #proposed immediate level ancestors
  SGGP$pala = matrix(0, nrow = SGGP$ML, ncol =SGGP$maxgridsize ) #proposedal all level ancestors
  SGGP$uala = matrix(0, nrow = SGGP$ML, ncol =SGGP$maxgridsize ) #used all level ancestors
  SGGP$pilaCOUNT = rep(0, SGGP$ML) #count of number of pila
  SGGP$palaCOUNT = rep(0, SGGP$ML) #count of number of pala
  SGGP$ualaCOUNT = rep(0, SGGP$ML) #count of number of uala
  
  SGGP$pilaCOUNT[1:SGGP$d] = 1
  SGGP$pila[1:SGGP$d, 1] = 1
  
  
  SGGP$bss = batchsize#1+4*SGGP$d  #must be at least 3*d
  SGGP$sizes = c(1,2,4,4,8,12,32) # Num of points added to 1D design as you go further in any dimension
  SGGP$maxlevel = length(SGGP$sizes)
  # Proposed grid size? More points further along the blocks?
  SGGP$pogsize = rep(0, 4 * SGGP$ML)
  SGGP$pogsize[1:SGGP$poCOUNT] = apply(matrix(SGGP$sizes[SGGP$po[1:SGGP$poCOUNT, ]], SGGP$poCOUNT, SGGP$d), 1, prod)
  # Selected sample size?
  SGGP$ss = 1
  
  
  SGGP$w = rep(0, SGGP$ML) #keep track of + and - for prediction
  SGGP$w[1] = 1 #keep track of + and - for prediction
  SGGP$uoCOUNT = 1 # Number of used levels
  # While number selected + min sample size <= batch size, i.e., still have enough spots for a block
  while (SGGP$bss > (SGGP$ss + min(SGGP$pogsize[1:SGGP$poCOUNT]) - 0.5)) {
    SGGP$uoCOUNT = SGGP$uoCOUNT + 1 #increment used count
    
    # First d iterations pick the (2,1,1,1,1),(1,2,1,1,1) blocks b/c we need info on each dimension before going adaptive
    if (SGGP$uoCOUNT < (SGGP$d + 1.5)) {
      pstar = 1 #pick a proposed to add
    } else{ # The next d iterations randomly pick from the boxes with minimal number of points, not sure this makes sense
      if (SGGP$uoCOUNT < (2 * SGGP$d + 1.5)) {
        pstar = sample(which(SGGP$pogsize[1:SGGP$poCOUNT] <= 0.5 + min(SGGP$pogsize[1:SGGP$poCOUNT])), 1)
      } else{ # After that randomly select from blocks that still fit
        pstar = sample(which(SGGP$pogsize[1:SGGP$poCOUNT] < min(SGGP$bss - SGGP$ss + 0.5,SGGP$maxgridsize)), 1)
      }
    }
    
    l0 =  SGGP$po[pstar, ] # Selected block e.g. (2,1,1,2)
    # Need to make sure there is still an open row in uo to set with new values
    if (SGGP$uoCOUNT > nrow(SGGP$uo)) {
      numrowstoadd <- 20
      SGGP$uo <- rbind(SGGP$uo, numrowstoadd)
      SGGP$ML <- nrow(SGGP$uo)
      
      # Need to get everything else upsized too
      SGGP$po = rbind(SGGP$po, matrix(0, nrow = 4 * numrowstoadd, ncol = ncol(SGGP$po))) #proposed levels tracker
      SGGP$pila = rbind(SGGP$pila, matrix(0, nrow = numrowstoadd, ncol=ncol(SGGP$pila))) #proposed immediate level ancestors
      SGGP$pala = rbind(SGGP$pala, matrix(0, nrow = numrowstoadd, ncol=ncol(SGGP$pala))) #proposedal all level ancestors
      SGGP$uala = rbind(SGGP$uala, matrix(0, nrow = numrowstoadd, ncol=ncol(SGGP$uala))) #used all level ancestors
      SGGP$pilaCOUNT = c(SGGP$pilaCOUNT, rep(0, numrowstoadd)) #count of number of pila
      SGGP$palaCOUNT = c(SGGP$palaCOUNT, rep(0, numrowstoadd)) #count of number of pala
      SGGP$ualaCOUNT = c(SGGP$ualaCOUNT, rep(0, numrowstoadd)) #count of number of uala
      SGGP$pogsize = c(SGGP$pogsize, rep(0, 4 * numrowstoadd))
      SGGP$w = c(SGGP$w, rep(0, numrowstoadd))
    }
    SGGP$uo[SGGP$uoCOUNT, ] = l0 # Store new block
    SGGP$ss =  SGGP$ss + SGGP$pogsize[pstar] # Update selected sample size
    
    # New ancestors?
    new_an = SGGP$pila[pstar, 1:SGGP$pilaCOUNT[pstar]]
    total_an = new_an
    
    # Loop over ancestors???
    for (anlcv in 1:length(total_an)) {
      # If more than one ancestor, update with unique ones.
      if (total_an[anlcv] > 1.5) {
        total_an = unique(c(total_an, SGGP$uala[total_an[anlcv], 1:SGGP$ualaCOUNT[total_an[anlcv]]]))
      }
    }
    # Update storage of ancestors
    SGGP$ualaCOUNT[SGGP$uoCOUNT]  = length(total_an)
    SGGP$uala[SGGP$uoCOUNT, 1:length(total_an)] = total_an
    
    # Loop over ancestors, update coefficient
    for (anlcv in 1:length(total_an)) {
      lo = SGGP$uo[total_an[anlcv], ]
      if (max(abs(lo - l0)) < 1.5) {
        SGGP$w[total_an[anlcv]] = SGGP$w[total_an[anlcv]] + (-1) ^ abs(round(sum(l0-lo)))
        
      }
    }
    SGGP$w[SGGP$uoCOUNT] = SGGP$w[SGGP$uoCOUNT] + 1
    
    # Update block tracking
    if (pstar < 1.5) { # If you picked the first block, update like this 
      SGGP$po[1:(SGGP$poCOUNT - 1), ] = SGGP$po[2:SGGP$poCOUNT, ]
      SGGP$pila[1:(SGGP$poCOUNT - 1), ] = SGGP$pila[2:SGGP$poCOUNT, ]
      SGGP$pilaCOUNT[1:(SGGP$poCOUNT - 1)] = SGGP$pilaCOUNT[2:SGGP$poCOUNT]
      SGGP$pogsize[1:(SGGP$poCOUNT - 1)] = SGGP$pogsize[2:SGGP$poCOUNT]
    }
    if (pstar > (SGGP$poCOUNT - 0.5)) { # If you picked the last block, do this
      SGGP$po[1:(SGGP$poCOUNT - 1), ] = SGGP$po[1:(pstar - 1), ]
      SGGP$pila[1:(SGGP$poCOUNT - 1), ] = SGGP$pila[1:(pstar - 1), ]
      SGGP$pilaCOUNT[1:(SGGP$poCOUNT - 1)] = SGGP$pilaCOUNT[1:(pstar - 1)]
      SGGP$pogsize[1:(SGGP$poCOUNT - 1)] = SGGP$pogsize[1:(pstar - 1)]
    }
    if (pstar < (SGGP$poCOUNT - 0.5) && pstar > 1.5) { # If in between, do this
      SGGP$po[1:(SGGP$poCOUNT - 1), ] = SGGP$po[c(1:(pstar - 1), (pstar + 1):SGGP$poCOUNT), ]
      SGGP$pila[1:(SGGP$poCOUNT - 1), ] = SGGP$pila[c(1:(pstar - 1), (pstar +
                                                                        1):SGGP$poCOUNT), ]
      SGGP$pilaCOUNT[1:(SGGP$poCOUNT - 1)] = SGGP$pilaCOUNT[c(1:(pstar - 1), (pstar +
                                                                                1):SGGP$poCOUNT)]
      SGGP$pogsize[1:(SGGP$poCOUNT - 1)] = SGGP$pogsize[c(1:(pstar - 1), (pstar +
                                                                            1):SGGP$poCOUNT)]
    }
    # One less option now???
    SGGP$poCOUNT = SGGP$poCOUNT - 1
    
    # Loop over dimensions WHY???
    for (dimlcv in 1:SGGP$d) {
      # The block e.g. (1,2,1,1,3) just selected
      lp = l0
      
      lp[dimlcv] = lp[dimlcv] + 1 # Increase THIS dim by 1. This is a new possibility, are we adding it?
      
      # Check if within some bounds??
      if (max(lp) <= SGGP$maxlevel && SGGP$poCOUNT < 4*SGGP$ML) {
        # Dimensions which are past first design level
        kvals = which(lp > 1.5)
        
        canuse = 1 # ????
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
  SGGP$xindex = 1:length(xb)
  # After this xb is
  #  [1] 0.50000 0.12500 0.87500 0.25000 0.75000 0.37500 0.62500 0.28125 0.71875 0.31250 0.68750 0.00000 1.00000 0.18750 0.81250
  # [16] 0.06250 0.93750 0.43750 0.56250 0.40625 0.59375 0.09375 0.90625 0.21875 0.78125 0.34375 0.65625 0.46875 0.53125 0.15625
  # [31] 0.84375 0.03125 0.96875
  SGGP$sizest = cumsum(SGGP$sizes) # Total number of points in 1D design as you go along axis
  
  
  SGGP$gridsizes = matrix(SGGP$sizes[SGGP$uo[1:SGGP$uoCOUNT, ]], SGGP$uoCOUNT, SGGP$d)
  SGGP$gridsizest = matrix(SGGP$sizest[SGGP$uo[1:SGGP$uoCOUNT, ]], SGGP$uoCOUNT, SGGP$d)
  SGGP$gridsize = apply(SGGP$gridsizes, 1, prod)
  SGGP$gridsizet = apply(SGGP$gridsizest, 1, prod)
  
  SGGP$di = matrix(0, nrow = SGGP$uoCOUNT, ncol = max(SGGP$gridsize))
  SGGP$dit = matrix(0, nrow = SGGP$uoCOUNT, ncol = sum((SGGP$gridsize)))
  
  SGGP$design = matrix(0, nrow = sum(SGGP$gridsize), ncol = SGGP$d)
  SGGP$designindex = matrix(0, nrow = sum(SGGP$gridsize), ncol = SGGP$d) # Use this to track which indices have been used
  tv = 0
  for (blocklcv in 1:SGGP$uoCOUNT) {
    SGGP$di[blocklcv, 1:SGGP$gridsize[blocklcv]] = (tv + 1):(tv + SGGP$gridsize[blocklcv])
    for (dimlcv in 1:SGGP$d) {
      levelnow = SGGP$uo[blocklcv, dimlcv]
      if (levelnow < 1.5) {
        SGGP$design[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(SGGP$xb[1], SGGP$gridsize[blocklcv])
        SGGP$designindex[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(SGGP$xindex[1], SGGP$gridsize[blocklcv])
      } else{
        x0 = SGGP$xb[(SGGP$sizest[levelnow - 1] + 1):SGGP$sizest[levelnow]]
        xi0 = SGGP$xindex[(SGGP$sizest[levelnow - 1] + 1):SGGP$sizest[levelnow]]
        if (dimlcv < 1.5) {
          SGGP$design[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(x0, "each" = SGGP$gridsize[blocklcv] /
                                                                               SGGP$gridsizes[blocklcv, dimlcv])
          SGGP$designindex[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(xi0, "each" = SGGP$gridsize[blocklcv] /
                                                                                    SGGP$gridsizes[blocklcv, dimlcv])
        }
        if (dimlcv > (SGGP$d - 0.5)) {
          SGGP$design[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(x0, SGGP$gridsize[blocklcv] /
                                                                               SGGP$gridsizes[blocklcv, dimlcv])
          SGGP$designindex[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(xi0, SGGP$gridsize[blocklcv] /
                                                                                    SGGP$gridsizes[blocklcv, dimlcv])
        }
        if (dimlcv < (SGGP$d - 0.5)  && dimlcv > 1.5) {
          SGGP$design[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(rep(x0, "each" =
                                                                                   prod(SGGP$gridsizes[blocklcv, (dimlcv + 1):SGGP$d])), prod(SGGP$gridsizes[blocklcv, 1:(dimlcv -1)]))
          SGGP$designindex[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(rep(xi0, "each" =
                                                                                        prod(SGGP$gridsizes[blocklcv, (dimlcv + 1):SGGP$d])), prod(SGGP$gridsizes[blocklcv, 1:(dimlcv -1)]))
        }
      }
    }
    
    tvv = 0
    if (blocklcv > 1.5) {
      for (ances in SGGP$uala[blocklcv, 1:SGGP$ualaCOUNT[blocklcv]]) {
        SGGP$dit[blocklcv, (tvv + 1):(tvv + SGGP$gridsize[ances])] = SGGP$di[ances, 1:SGGP$gridsize[ances]]
        tvv = tvv + SGGP$gridsize[ances]
      }
      SGGP$dit[blocklcv, (tvv + 1):(tvv + SGGP$gridsize[blocklcv])] = SGGP$di[blocklcv, 1:SGGP$gridsize[blocklcv]]
      Xset = SGGP$design[SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]], ]
      reorder = do.call(order, lapply(1:NCOL(Xset), function(kvt) Xset[, kvt]))
      SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]] = SGGP$dit[blocklcv, reorder]
    } else{
      SGGP$dit[blocklcv, 1:SGGP$gridsize[blocklcv]] = SGGP$di[blocklcv, 1:SGGP$gridsize[blocklcv]]
    }
    
    SGGP$design_unevaluated <- SGGP$design
    
    tv = tv + SGGP$gridsize[blocklcv]
  }
  return(SGGP)
}