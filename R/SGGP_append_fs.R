#' Calculate MSE over single dimension for Matern 3/2
#' 
#' Calcaluted using grid of integration points.
#' Can be calculated exactly, but not much reason in 1D.
#'
#' @param xl Vector of points in 1D
#' @param theta Log of correlation parameters.
#' @param theta Correlation parameters
#' @param nugget Nugget to add to diagonal of correlation matrix.
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#' @param diag_corrMat Function that gives diagonal of correlation matrix for vector of 1D points.
#' @param ... Don't use, just forces theta to be named
#'
#' @return MSE value
#' @export
#'
#' @examples
#' MSE_calc(xl=c(0,.5,.9), theta=1, nugget=.001,
#'          CorrMat=CorrMatMatern32,
#'          diag_corrMat=diag_corrMatMatern32)
SGGP_internal_MSEcalc <- function(xl, theta, CorrMat) {
  S = CorrMat(xl, xl, theta)
  xp = seq(0,1,l=101)
  Cp = CorrMat(xp,xl,theta)
  n = length(xl)
  cholS = chol(S)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  
  MSE_val = mean(1 - rowSums(t(CiCp)*Cp))
  
  MSE_val
}


#' Calculate MSE over blocks
#' 
#' Delta of adding block is product over i=1..d of IMSE(i,j-1) - IMSE(i,j)
#'
#' @param valsinds Block levels to calculate MSEs for
#' @param MSE_v Matrix of MSE values
#'
#' @return All MSE values
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' theta <- c(.1,.1,.1)
#' MSE_v <- outer(1:SGGP$d, 1:8, 
#'      Vectorize(function(lcv1, dimlcv) {
#'         MSE_calc(SGGP$xb[1:SGGP$sizest[dimlcv]], theta=theta[lcv1], nugget=0,
#'          CorrMat=CorrMatMatern32,
#'          diag_corrMat=diag_corrMatMatern32)
#'  }))
#' MSE_de(SGGP$po[1:SGGP$poCOUNT, ], MSE_v)
SGGP_internal_MSEde <- function(valsinds, MSE_v) {
  if(is.matrix(valsinds)){
    MSE_de = rep(0, dim(valsinds)[1])
    
    for (lcv1 in 1:dim(valsinds)[1]) {
      MSE_de[lcv1] = 0
      
      for (dimlcv in 1:dim(valsinds)[2]) {
        if (valsinds[lcv1, dimlcv] > 1.5) {
          MSE_de[lcv1] = MSE_de[lcv1] + log(-MSE_v[dimlcv, valsinds[lcv1, dimlcv]] + MSE_v[dimlcv, valsinds[lcv1, dimlcv] - 1])
          
        } else {
          # This is when no ancestor block, 1 comes from when there is no data. 
          # 1 is correlation times integrated value over range.
          # This depends on correlation function.
          MSE_de[lcv1] = MSE_de[lcv1] + log(-MSE_v[dimlcv, valsinds[lcv1, dimlcv]] + 1)
          
        }
      }
    }
  } else {
    MSE_de = 0
    
    for (dimlcv in 1:length(valsinds)) {
      if (valsinds[dimlcv] > 1.5) {
        MSE_de = MSE_de + log(-MSE_v[dimlcv, valsinds[dimlcv]] + MSE_v[dimlcv, valsinds[dimlcv] -1])
        
      } else {
        MSE_de = MSE_de + log(-MSE_v[dimlcv, valsinds[dimlcv]] + 1)
        
      }
    }}
  MSE_de = exp(MSE_de)
  
  MSE_de # CBE added this so it will return normally.
}



#' Add points to SGGP
#' 
#' Add `batchsize` points to `SG` using `theta`.
#'
#' @param SG Sparse grid object
#' @param batchsize Number of points to add
#' @param theta Correlation parameters
#' @param theta Log of theta, give one of theta and theta
#' @param ... Don't use, just forces theta to be named
#'
#' @return SG with new points added.
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' SG <- SGappend(theta=c(.1,.1,.1), SG=SG, batchsize=20)
SGGPappend <- function(SGGP,batchsize){
  n_before <- nrow(SGGP$design)
  
  # Set up blank matrix to store MSE values
  MSE_v = matrix(0, SGGP$d, SGGP$maxgridsize) # 8 because he only defined the 1D designs up to 8.
  # Why do we consider dimensions independent of each other?
  # Loop over dimensions and design refinements
  for (dimlcv in 1:SGGP$d) {
    for (levellcv in 1:SGGP$maxlevel) {
      # Calculate some sort of MSE from above, not sure what it's doing
      MSE_v[dimlcv, levellcv] = max(0, abs(MSE_calc(SGGP$xb[1:SGGP$sizest[levellcv]],SGGP$theta[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],SGGP$CorrMat)))
      if (levellcv > 1.5) { # If past first level, it is as good as one below it. Why isn't this a result of calculation?
        MSE_v[dimlcv, levellcv] = min(MSE_v[dimlcv, levellcv], MSE_v[dimlcv, levellcv - 1])
      }
    }
  }
  
  # What is this? Integrate MSE
  I_mes = rep(0, SGGP$ML)
  
  # For all possible blocks, calculate MSE_v? Is that all that MSE_de does?
  I_mes[1:SGGP$poCOUNT] = MSE_de(SGGP$po[1:SGGP$poCOUNT, ], MSE_v)
  
  # Increase count of points evaluated. Do we check this if not reached exactly???
  SGGP$bss = SGGP$bss + batchsize
  
  # Keep adding points until reaching bss
  while (SGGP$bss > (SGGP$ss + min(SGGP$pogsize[1:SGGP$poCOUNT]) - 0.5)) {
    SGGP$uoCOUNT = SGGP$uoCOUNT + 1 #increment used count
    # Find the best one that still fits
    M_comp = max(I_mes[which(SGGP$pogsize[1:SGGP$poCOUNT] < (SGGP$bss - SGGP$ss + 0.5))])
    # Find which ones are close to M_comp and
    possibleO =which((I_mes[1:SGGP$poCOUNT] >= 0.5*M_comp)&(SGGP$pogsize[1:SGGP$poCOUNT] < (SGGP$bss - SGGP$ss + 0.5)))
    # If more than one is possible, randomly pick among them.
    if(length(possibleO)>1.5){
      pstar = sample(possibleO,1)
    } else{
      pstar = possibleO
    }
    
    l0 =  SGGP$po[pstar,] # Selected block
    SGGP$uo[SGGP$uoCOUNT,] = l0 # Save selected block
    SGGP$ss =  SGGP$ss + SGGP$pogsize[pstar] # Update selected size
    
    # New ancestors???
    new_an = SGGP$pila[pstar, 1:SGGP$pilaCOUNT[pstar]]
    total_an = new_an
    for (anlcv in 1:length(total_an)) { # Loop over ancestors
      if (total_an[anlcv] > 1.5) { # If there's more than 1, do ???
        total_an = unique(c(total_an, SGGP$uala[total_an[anlcv], 1:SGGP$ualaCOUNT[total_an[anlcv]]]))
      }
    }
    SGGP$ualaCOUNT[SGGP$uoCOUNT]  = length(total_an)
    SGGP$uala[SGGP$uoCOUNT, 1:length(total_an)] = total_an
    
    # Loop over all ancestors, why???
    for (anlcv in 1:length(total_an)) {
      lo = SGGP$uo[total_an[anlcv],]
      if (max(abs(lo - l0)) < 1.5) {
        SGGP$w[total_an[anlcv]] = SGGP$w[total_an[anlcv]] + (-1) ^ abs(round(sum(l0-lo)))
        
      }
    }
    SGGP$w[SGGP$uoCOUNT] = SGGP$w[SGGP$uoCOUNT] + 1
    
    
    # If on the first block
    if (pstar < 1.5) {
      SGGP$po[1:(SGGP$poCOUNT - 1),] = SGGP$po[2:SGGP$poCOUNT,]
      SGGP$pila[1:(SGGP$poCOUNT - 1),] = SGGP$pila[2:SGGP$poCOUNT,]
      SGGP$pilaCOUNT[1:(SGGP$poCOUNT - 1)] = SGGP$pilaCOUNT[2:SGGP$poCOUNT]
      SGGP$pogsize[1:(SGGP$poCOUNT - 1)] = SGGP$pogsize[2:SGGP$poCOUNT]
      I_mes[1:(SGGP$poCOUNT - 1)] = I_mes[2:SGGP$poCOUNT]
    }
    if (pstar > (SGGP$poCOUNT - 0.5)) {
      SGGP$po[1:(SGGP$poCOUNT - 1),] = SGGP$po[1:(pstar - 1),]
      SGGP$pila[1:(SGGP$poCOUNT - 1),] = SGGP$pila[1:(pstar - 1),]
      SGGP$pilaCOUNT[1:(SGGP$poCOUNT - 1)] = SGGP$pilaCOUNT[1:(pstar - 1)]
      SGGP$pogsize[1:(SGGP$poCOUNT - 1)] = SGGP$pogsize[1:(pstar - 1)]
      I_mes[1:(SGGP$poCOUNT - 1)] = I_mes[1:(pstar - 1)]
    }
    if (pstar < (SGGP$poCOUNT - 0.5) && pstar > 1.5) {
      SGGP$po[1:(SGGP$poCOUNT - 1),] = SGGP$po[c(1:(pstar - 1), (pstar + 1):SGGP$poCOUNT),]
      SGGP$pila[1:(SGGP$poCOUNT - 1),] = SGGP$pila[c(1:(pstar - 1), (pstar +1):SGGP$poCOUNT),]
      SGGP$pilaCOUNT[1:(SGGP$poCOUNT - 1)] = SGGP$pilaCOUNT[c(1:(pstar - 1), (pstar + 1):SGGP$poCOUNT)]
      SGGP$pogsize[1:(SGGP$poCOUNT - 1)] = SGGP$pogsize[c(1:(pstar - 1), (pstar + 1):SGGP$poCOUNT)]
      I_mes[1:(SGGP$poCOUNT - 1)] = I_mes[c(1:(pstar - 1), (pstar + 1):SGGP$poCOUNT)]
    }
    SGGP$poCOUNT = SGGP$poCOUNT - 1
    
    for (dimlcv in 1:SGGP$d) {
      lp = l0
      
      lp[dimlcv] = lp[dimlcv] + 1
      
      if (max(lp) < SGGP$maxlevel && SGGP$poCOUNT < 4 * SGGP$ML) {
        kvals = which(lp > 1.5)
        
        canuse = 1
        ap = rep(0, SGGP$d)
        nap = 0
        for (activedimlcv in 1:length(kvals)) {
          lpp = lp
          lpp[kvals[activedimlcv]] = lpp[kvals[activedimlcv]] - 1
          
          ismem = rep(1, SGGP$uoCOUNT)
          for (dimdimlcv in 1:SGGP$d) {
            ismem  = ismem * (SGGP$uo[1:SGGP$uoCOUNT, dimdimlcv] == lpp[dimdimlcv])
          }
          
          if (max(ismem) > 0.5) {
            ap[activedimlcv] = which(ismem > 0.5)
            nap = nap + 1
          } else{
            canuse = 0
          }
        }
        if (canuse > 0.5) {
          SGGP$poCOUNT = SGGP$poCOUNT + 1
          SGGP$po[SGGP$poCOUNT,] = lp
          SGGP$pogsize[SGGP$poCOUNT] = prod(SGGP$sizes[lp])
          SGGP$pila[SGGP$poCOUNT, 1:nap] = ap[1:nap]
          SGGP$pilaCOUNT[SGGP$poCOUNT] = nap
          I_mes[SGGP$poCOUNT] =  MSE_de(as.vector(SGGP$po[SGGP$poCOUNT, ]), MSE_v)
        }
      }
    }
  }
  
  SGGP$gridsizes = matrix(SGGP$sizes[SGGP$uo[1:SGGP$uoCOUNT, ]], SGGP$uoCOUNT, SGGP$d)
  SGGP$gridsizest = matrix(SGGP$sizest[SGGP$uo[1:SGGP$uoCOUNT, ]], SGGP$uoCOUNT, SGGP$d)
  SGGP$gridsize = apply(SGGP$gridsizes, 1, prod)
  SGGP$gridsizet = apply(SGGP$gridsizest, 1, prod)
  
  SGGP$di = matrix(0, nrow = SGGP$uoCOUNT, ncol = max(SGGP$gridsize))
  SGGP$dit = matrix(0, nrow = SGGP$uoCOUNT, ncol = sum((SGGP$gridsize)))
  
  # THIS OVERWRITES AND RECALCULATES design EVERY TIME, WHY NOT JUST DO FOR NEW ROWS?
  SGGP$design = matrix(0, nrow = sum(SGGP$gridsize), ncol = SGGP$d)
  SGGP$designindex = matrix(0, nrow = sum(SGGP$gridsize), ncol = SGGP$d)
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
                                                                         prod(SGGP$gridsizes[blocklcv, (dimlcv + 1):SGGP$d])), prod(SGGP$gridsizes[blocklcv, 1:(dimlcv -
                                                                                                                                                  1)]))
          SGGP$designindex[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(rep(xi0, "each" =
                                                                              prod(SGGP$gridsizes[blocklcv, (dimlcv + 1):SGGP$d])), prod(SGGP$gridsizes[blocklcv, 1:(dimlcv -
                                                                                                                                                       1)]))
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
      reorder = do.call(order, lapply(1:NCOL(Xset), function(kvt)
        Xset[, kvt]))
      SGGP$dit[blocklcv, 1:SGGP$gridsizet[blocklcv]] = SGGP$dit[blocklcv, reorder]
    } else{
      SGGP$dit[blocklcv, 1:SGGP$gridsize[blocklcv]] = SGGP$di[blocklcv, 1:SGGP$gridsize[blocklcv]]
    }
    
    tv = tv + SGGP$gridsize[blocklcv]
  }
  
  # Save Xnew to make it easy to know which ones to add
  SGGP$Xnew <- SGGP$design[(n_before+1):nrow(SGGP$design),]
  
  return(SG)
}
