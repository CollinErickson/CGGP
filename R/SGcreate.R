#' Create sparse grid GP
#'
#' @param d Input dimension
#' @param xmin Min x values, vector. Must be rep(0,d).
#' @param xmax Max x values, vector. Must be rep(1,d).
#' @param batchsize Number added to design each batch
#' @param nugget Nugget term added to diagonal of correlation matrix,
#' for now only on predictions
#' @param corr Correlation function to use.
#'
#' @return SGGP
#' @export
#'
#' @examples
#' d <- 8
#' SG = SGcreate(d,201)
SGcreate <- function(d, batchsize, corr="Matern32", nugget=0, xmin=rep(0,d), xmax=rep(1,d)) {
  # This is list representing our GP object
  SG = list("xmin" = xmin, "xmax" = xmax)
  class(SG) <- c("SGGP", "list")
  if (tolower(corr) %in% c('matern32', 'mat32', 'm32', 'matern3', 'mat3', 'm3')) {
    SG$CorrMat <- CorrMatMatern32
    SG$dCorrMat <- dCorrMatMatern32
    SG$diag_corrMat <- diag_corrMatMatern32
    SG$ddiag_corrMat <- ddiag_corrMatMatern32
  } else if (tolower(corr) %in% c('gaussian', 'gauss')) {
    SG$CorrMat <- gausscorr
    SG$dCorrMat <- dgausscorr
    SG$diag_corrMat <- diag_gausscorr
  } else {
    stop(paste0("corr given was ', corr,', must be one of matern32."))
  }
  SG$nugget <- nugget
  
  if (any(xmin!=0) || any(xmax!=1)) {stop("For now must use xmin=0 and xmax=1")}
  SG$d = length(xmin) # input dimension
  
  # Levels are blocks. Level is like eta from paper.
  SG$ML = min(choose(SG$d + 6, SG$d), 100000) #max levels
  
  # What is levelpoint? The current point? This is not used again in this file!
  SG$levelpoint = rep(0, SG$ML)
  
  # Track evaluated blocks, aka used levels
  SG$uo = matrix(0, nrow = SG$ML, ncol = SG$d) #used levels tracker
  # First observation is always (1,1,...,1)
  SG$uo[1, ] = rep(1, SG$d) #first observation in middle of space
  SG$uoCOUNT = 1 #number of used levels
  
  # Track the blocks that are allowed to be evaluated
  SG$po = matrix(0, nrow = 4 * SG$ML, ncol = SG$d) #proposed levels tracker
  # Now possible blocks are (2,1,1,1,1), (1,2,1,1,1), (1,1,2,1,1), etc
  SG$po[1:SG$d, ] = matrix(1, nrow = SG$d, ncol = SG$d) + diag(SG$d) #one at a time
  SG$poCOUNT = SG$d #number of proposed levels
  
  
  # What are ancestors? Why are we doing this?
  # Are ancestors support blocks? Is this for calculating coefficient?
  # How any of its ancestors be proposed? They should all be used already?
  # Might be transposed??? Is this right?
  SG$pila = matrix(0, nrow = SG$ML, ncol = 1000) #proposed immediate level ancestors
  SG$pala = matrix(0, nrow = SG$ML, ncol = 1000) #proposedal all level ancestors
  SG$uala = matrix(0, nrow = SG$ML, ncol = 1000) #used all level ancestors
  SG$pilaCOUNT = rep(0, SG$ML) #count of number of pila
  SG$palaCOUNT = rep(0, SG$ML) #count of number of pala
  SG$ualaCOUNT = rep(0, SG$ML) #count of number of uala
  
  SG$pilaCOUNT[1:SG$d] = 1
  SG$pila[1:SG$d, 1] = 1
  
  SG$bss = batchsize#1+4*SG$d  #must be at least 3*d
  SG$sizes = c(1, 2, 2, 3, 3, 4, 4, 6, 8) # Num of points added to 1D design as you go further in any dimension
  # Proposed grid size? More points further along the blocks?
  SG$pogsize = rep(0, 4 * SG$ML)
  SG$pogsize[1:SG$poCOUNT] = apply(matrix(SG$sizes[SG$po[1:SG$poCOUNT, ]], SG$poCOUNT, SG$d), 1, prod)
  # Selected sample size?
  SG$ss = 1
  
  
  SG$w = rep(0, SG$ML) #keep track of + and - for prediction
  SG$w[1] = 1 #keep track of + and - for prediction
  SG$uoCOUNT = 1 # Number of used levels
  # While number selected + min sample size <= batch size, i.e., still have enough spots for a block
  while (SG$bss > (SG$ss + min(SG$pogsize[1:SG$poCOUNT]) - 0.5)) {
    SG$uoCOUNT = SG$uoCOUNT + 1 #increment used count
    
    # First d iterations pick the (2,1,1,1,1),(1,2,1,1,1) blocks b/c we need info on each dimension before going adaptive
    if (SG$uoCOUNT < (SG$d + 1.5)) {
      pstar = 1 #pick a proposed to add
    } else{ # The next d iterations randomly pick from the boxes with minimal number of points, not sure this makes sense
      if (SG$uoCOUNT < (2 * SG$d + 1.5)) {
        pstar = sample(which(SG$pogsize[1:SG$poCOUNT] <= 0.5 + min(SG$pogsize[1:SG$poCOUNT])), 1)
      } else{ # After that randomly select from blocks that still fit
        pstar = sample(which(SG$pogsize[1:SG$poCOUNT] < (SG$bss - SG$ss + 0.5)), 1)
      }
    }
    
    l0 =  SG$po[pstar, ] # Block name e.g. (2,1,1,2)
    SG$uo[SG$uoCOUNT, ] = l0 # Store new block
    SG$ss =  SG$ss + SG$pogsize[pstar] # Update selected sample size
    
    # New ancestors?
    new_an = SG$pila[pstar, 1:SG$pilaCOUNT[pstar]]
    total_an = new_an
    
    # Loop over ancestors???
    for (lcv2 in 1:length(total_an)) {
      # If more than one ancestor, update with unique ones.
      if (total_an[lcv2] > 1.5) {
        total_an = unique(c(total_an, SG$uala[total_an[lcv2], 1:SG$ualaCOUNT[total_an[lcv2]]]))
      }
    }
    # Update storage of ancestors
    SG$ualaCOUNT[SG$uoCOUNT]  = length(total_an)
    SG$uala[SG$uoCOUNT, 1:length(total_an)] = total_an
    
    # Loop over ancestors, update coefficient
    for (lcv2 in 1:length(total_an)) {
      lo = SG$uo[total_an[lcv2], ]
      if (max(abs(lo - l0)) < 1.5) {
        SG$w[total_an[lcv2]] = SG$w[total_an[lcv2]] + (-1) ^ abs(round(sum(l0 -
                                                                             lo)))
        
      }
    }
    SG$w[SG$uoCOUNT] = SG$w[SG$uoCOUNT] + 1
    
    # Update block tracking
    if (pstar < 1.5) { # If you picked the first block, update like this 
      SG$po[1:(SG$poCOUNT - 1), ] = SG$po[2:SG$poCOUNT, ]
      SG$pila[1:(SG$poCOUNT - 1), ] = SG$pila[2:SG$poCOUNT, ]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[2:SG$poCOUNT]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[2:SG$poCOUNT]
    }
    if (pstar > (SG$poCOUNT - 0.5)) { # If you picked the last block, do this
      SG$po[1:(SG$poCOUNT - 1), ] = SG$po[1:(pstar - 1), ]
      SG$pila[1:(SG$poCOUNT - 1), ] = SG$pila[1:(pstar - 1), ]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[1:(pstar - 1)]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[1:(pstar - 1)]
    }
    if (pstar < (SG$poCOUNT - 0.5) && pstar > 1.5) { # If in between, do this
      SG$po[1:(SG$poCOUNT - 1), ] = SG$po[c(1:(pstar - 1), (pstar + 1):SG$poCOUNT), ]
      SG$pila[1:(SG$poCOUNT - 1), ] = SG$pila[c(1:(pstar - 1), (pstar +
                                                                  1):SG$poCOUNT), ]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[c(1:(pstar - 1), (pstar +
                                                                          1):SG$poCOUNT)]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[c(1:(pstar - 1), (pstar +
                                                                      1):SG$poCOUNT)]
    }
    # One less option now???
    SG$poCOUNT = SG$poCOUNT - 1
    
    # Loop over dimensions WHY???
    for (lcv2 in 1:SG$d) {
      # The block e.g. (1,2,1,1,3) just selected
      lp = l0
      
      lp[lcv2] = lp[lcv2] + 1 # Increase THIS dim by 1. This is a new possibility, are we adding it?
      
      # Check if within some bounds??
      if (max(lp) < 7.5 && SG$poCOUNT < 4 * SG$ML) {
        # Dimensions which are past first design level
        kvals = which(lp > 1.5)
        
        canuse = 1 # ????
        ap = rep(0, SG$d) # ????
        nap = 0 # ?????
        
        # Loop over dims at 2+ and do what?
        for (lcv3 in 1:length(kvals)) {
          lpp = lp # The block selected with 1 dim incremented
          lpp[kvals[lcv3]] = lpp[kvals[lcv3]] - 1 # ????
          
          ismem = rep(1, SG$uoCOUNT) # Boolean???
          # Loop over dimensions
          for (lcv4 in 1:SG$d) { # Set to 0 or 1 if all points already selected have same value???????
            ismem  = ismem * (SG$uo[1:SG$uoCOUNT, lcv4] == lpp[lcv4])
          }
          # If any are still 1,
          if (max(ismem) > 0.5) {
            ap[lcv3] = which(ismem > 0.5)
            nap = nap + 1 # Count number that are >=1
          } else{ # All are 0, so can't use
            canuse = 0
          }
        }
        # If it can be used, add to possible blocks
        if (canuse > 0.5) {
          SG$poCOUNT = SG$poCOUNT + 1
          SG$po[SG$poCOUNT, ] = lp
          SG$pogsize[SG$poCOUNT] = prod(SG$sizes[lp])
          SG$pila[SG$poCOUNT, 1:nap] = ap[1:nap]
          SG$pilaCOUNT[SG$poCOUNT] = nap
          
        }
      }
    }
  }
  
  # Create points for design
  xb = rep(
    c(
      3 / 8,
      1 / 4,
      1 / 8,
      7 / 32,
      3 / 16,
      1 / 2,
      5 / 16,
      7 / 16,
      1 / 16,
      3 / 32,
      13 / 32,
      9 / 32,
      5 / 32,
      1 / 32,
      11 / 32,
      15 / 32
    ),
    "each" = 2
  )
  SG$xb = 0.5 + c(0, xb * rep(c(-1, 1), length(xb) / 2))
  SG$xindex = 1:length(xb)
  # After this xb is
  #  [1] 0.50000 0.12500 0.87500 0.25000 0.75000 0.37500 0.62500 0.28125 0.71875 0.31250 0.68750 0.00000 1.00000 0.18750 0.81250
  # [16] 0.06250 0.93750 0.43750 0.56250 0.40625 0.59375 0.09375 0.90625 0.21875 0.78125 0.34375 0.65625 0.46875 0.53125 0.15625
  # [31] 0.84375 0.03125 0.96875
  SG$sizest = cumsum(SG$sizes) # Total number of points in 1D design as you go along axis
  
  
  SG$gridsizes = matrix(SG$sizes[SG$uo[1:SG$uoCOUNT, ]], SG$uoCOUNT, SG$d)
  SG$gridsizest = matrix(SG$sizest[SG$uo[1:SG$uoCOUNT, ]], SG$uoCOUNT, SG$d)
  SG$gridsize = apply(SG$gridsizes, 1, prod)
  SG$gridsizet = apply(SG$gridsizest, 1, prod)
  
  SG$di = matrix(0, nrow = SG$uoCOUNT, ncol = max(SG$gridsize))
  SG$dit = matrix(0, nrow = SG$uoCOUNT, ncol = sum((SG$gridsize)))
  
  SG$design = matrix(0, nrow = sum(SG$gridsize), ncol = SG$d)
  SG$designindex = matrix(0, nrow = sum(SG$gridsize), ncol = SG$d) # Use this to track which indices have been used
  tv = 0
  for (lcv1 in 1:SG$uoCOUNT) {
    SG$di[lcv1, 1:SG$gridsize[lcv1]] = (tv + 1):(tv + SG$gridsize[lcv1])
    for (lcv2 in 1:SG$d) {
      levelnow = SG$uo[lcv1, lcv2]
      if (levelnow < 1.5) {
        SG$design[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(SG$xb[1], SG$gridsize[lcv1])
        SG$designindex[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(SG$xindex[1], SG$gridsize[lcv1])
      } else{
        x0 = SG$xb[(SG$sizest[levelnow - 1] + 1):SG$sizest[levelnow]]
        xi0 = SG$xindex[(SG$sizest[levelnow - 1] + 1):SG$sizest[levelnow]]
        if (lcv2 < 1.5) {
          SG$design[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(x0, "each" = SG$gridsize[lcv1] /
                                                                     SG$gridsizes[lcv1, lcv2])
          SG$designindex[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(xi0, "each" = SG$gridsize[lcv1] /
                                                                          SG$gridsizes[lcv1, lcv2])
        }
        if (lcv2 > (SG$d - 0.5)) {
          SG$design[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(x0, SG$gridsize[lcv1] /
                                                                     SG$gridsizes[lcv1, lcv2])
          SG$designindex[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(xi0, SG$gridsize[lcv1] /
                                                                          SG$gridsizes[lcv1, lcv2])
        }
        if (lcv2 < (SG$d - 0.5)  && lcv2 > 1.5) {
          SG$design[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(rep(x0, "each" =
                                                                         prod(SG$gridsizes[lcv1, (lcv2 + 1):SG$d])), prod(SG$gridsizes[lcv1, 1:(lcv2 -
                                                                                                                                                  1)]))
          SG$designindex[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(rep(xi0, "each" =
                                                                              prod(SG$gridsizes[lcv1, (lcv2 + 1):SG$d])), prod(SG$gridsizes[lcv1, 1:(lcv2 -
                                                                                                                                                       1)]))
        }
      }
    }
    
    tvv = 0
    if (lcv1 > 1.5) {
      for (ances in SG$uala[lcv1, 1:SG$ualaCOUNT[lcv1]]) {
        SG$dit[lcv1, (tvv + 1):(tvv + SG$gridsize[ances])] = SG$di[ances, 1:SG$gridsize[ances]]
        tvv = tvv + SG$gridsize[ances]
      }
      SG$dit[lcv1, (tvv + 1):(tvv + SG$gridsize[lcv1])] = SG$di[lcv1, 1:SG$gridsize[lcv1]]
      Xset = SG$design[SG$dit[lcv1, 1:SG$gridsizet[lcv1]], ]
      reorder = do.call(order, lapply(1:NCOL(Xset), function(kvt)
        Xset[, kvt]))
      SG$dit[lcv1, 1:SG$gridsizet[lcv1]] = SG$dit[lcv1, reorder]
    } else{
      SG$dit[lcv1, 1:SG$gridsize[lcv1]] = SG$di[lcv1, 1:SG$gridsize[lcv1]]
    }
    
    tv = tv + SG$gridsize[lcv1]
  }
  
  # Save predictive weights Rinv*y
  SG$pw <- NULL
  
  return(SG)
}