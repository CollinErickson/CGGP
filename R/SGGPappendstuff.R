# No clue what this one is, looks like it's a specific function, but a real mess.
# But it's actually used in SGappend. This is insane.
# Maybe it's a convoluted way to calculated MSE over entire region when using given correlation function.
# Only 1 dimension.

#' Calculate MSE over single dimension for Matern 3/2
#' 
#' Calcaluted using grid of integration points.
#' Can be calculated exactly, but not much reason in 1D.
#'
#' @param xl Vector of points in 1D
#' @param logtheta Log of correlation parameters.
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
MSE_calc <- function(xl, ..., logtheta, theta, nugget, CorrMat, diag_corrMat) {
  if (missing(theta)) {theta <- exp(logtheta)}
  S = CorrMat(xl, xl, theta=theta)
  xp = seq(0,1,l=101)
  Sn = CorrMat(xp, xl, theta=theta)
  diag(S) <- diag(S) + nugget
  
  n = length(xl)
  Ci = solve(S, t(Sn))
  
  #MSE = diag_corrMat(.5, theta=theta, nugget=nugget) - sum(diag(out1 %*% Ci))
  MSE = mean(diag_corrMat(xp, theta=theta, nugget=nugget) - rowSums(Sn * t(Ci)))
  
  MSE
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
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' theta <- c(.1,.1,.1)
#' MSE_v <- outer(1:SG$d, 1:8, 
#'      Vectorize(function(lcv1, lcv2) {
#'         MSE_calc(SG$xb[1:SG$sizest[lcv2]], theta=theta[lcv1], nugget=0,
#'          CorrMat=CorrMatMatern32,
#'          diag_corrMat=diag_corrMatMatern32)
#'  }))
#' MSE_de(SG$po[1:SG$poCOUNT, ], MSE_v)
MSE_de <- function(valsinds, MSE_v) {
  if(is.matrix(valsinds)){
    MSE_de = rep(0, dim(valsinds)[1])
    
    for (lcv1 in 1:dim(valsinds)[1]) {
      MSE_de[lcv1] = 0
      
      for (lcv2 in 1:dim(valsinds)[2]) {
        if (valsinds[lcv1, lcv2] > 1.5) {
          MSE_de[lcv1] = MSE_de[lcv1] + log(-MSE_v[lcv2, valsinds[lcv1, lcv2]] + MSE_v[lcv2, valsinds[lcv1, lcv2] - 1])
          
        } else {
          # This is when no ancestor block, 1 comes from when there is no data. 
          # 1 is correlation times integrated value over range.
          # This depends on correlation function.
          MSE_de[lcv1] = MSE_de[lcv1] + log(-MSE_v[lcv2, valsinds[lcv1, lcv2]] + 1)
          
        }
      }
    }
  } else {
    MSE_de = 0
    
    for (lcv2 in 1:length(valsinds)) {
      if (valsinds[lcv2] > 1.5) {
        MSE_de = MSE_de + log(-MSE_v[lcv2, valsinds[lcv2]] + MSE_v[lcv2, valsinds[lcv2] -1])
        
      } else {
        MSE_de = MSE_de + log(-MSE_v[lcv2, valsinds[lcv2]] + 1)
        
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
#' @param ... Don't use, just forces theta to be named
#'
#' @return SG with new points added.
#' @export
#'
#' @examples
#' SG <- SGcreate(c(0,0,0), c(1,1,1), batchsize=100)
#' SG <- SGappend(theta=c(.1,.1,.1), SG=SG, batchsize=20)
SGappend <- function(SG,batchsize,..., theta){
  
  # Set up blank matrix to store MSE values
  MSE_v = matrix(0, SG$d, 8) # 8 because he only defined the 1D designs up to 8.
  # Why do we consider dimensions independent of each other?
  # Loop over dimensions and design refinements
  for (lcv1 in 1:SG$d) {
    for (lcv2 in 1:8) {
      # Calculate some sort of MSE from above, not sure what it's doing
      MSE_v[lcv1, lcv2] = max(0, abs(MSE_calc(SG$xb[1:SG$sizest[lcv2]], theta=theta[lcv1], nugget=SG$nugget,
                                              CorrMat=SG$CorrMat,
                                              diag_corrMat=SG$diag_corrMat)))
      if (lcv2 > 1.5) { # If past first level, it is as good as one below it. Why isn't this a result of calculation?
        MSE_v[lcv1, lcv2] = min(MSE_v[lcv1, lcv2], MSE_v[lcv1, lcv2 - 1])
      }
    }
  }
  
  # What is this? Integrate MSE
  I_mes = rep(0, SG$ML)
  
  # For all possible blocks, calculate MSE_v? Is that all that MSE_de does?
  I_mes[1:SG$poCOUNT] = MSE_de(SG$po[1:SG$poCOUNT, ], MSE_v)
  
  # Increase count of points evaluated. Do we check this if not reached exactly???
  SG$bss = SG$bss + batchsize
  
  # Keep adding points until reaching bss
  while (SG$bss > (SG$ss + min(SG$pogsize[1:SG$poCOUNT]) - 0.5)) {
    SG$uoCOUNT = SG$uoCOUNT + 1 #increment used count
    # Find the best one that still fits
    M_comp = max(I_mes[which(SG$pogsize[1:SG$poCOUNT] < (SG$bss - SG$ss + 0.5))])
    # Find which ones are close to M_comp and
    possibleO =which((I_mes[1:SG$poCOUNT] >= 0.5*M_comp)&(SG$pogsize[1:SG$poCOUNT] < (SG$bss - SG$ss + 0.5)))
    # If more than one is possible, randomly pick among them.
    if(length(possibleO)>1.5){
      pstar = sample(possibleO,1)
    } else{
      pstar = possibleO
    }
    
    l0 =  SG$po[pstar,] # Selected block
    SG$uo[SG$uoCOUNT,] = l0 # Save selected block
    SG$ss =  SG$ss + SG$pogsize[pstar] # Update selected size
    
    # New ancestors???
    new_an = SG$pila[pstar, 1:SG$pilaCOUNT[pstar]]
    total_an = new_an
    for (lcv2 in 1:length(total_an)) { # Loop over ancestors
      if (total_an[lcv2] > 1.5) { # If there's more than 1, do ???
        total_an = unique(c(total_an, SG$uala[total_an[lcv2], 1:SG$ualaCOUNT[total_an[lcv2]]]))
      }
    }
    SG$ualaCOUNT[SG$uoCOUNT]  = length(total_an)
    SG$uala[SG$uoCOUNT, 1:length(total_an)] = total_an
    
    # Loop over all ancestors, why???
    for (lcv2 in 1:length(total_an)) {
      lo = SG$uo[total_an[lcv2],]
      if (max(abs(lo - l0)) < 1.5) {
        SG$w[total_an[lcv2]] = SG$w[total_an[lcv2]] + (-1) ^ abs(round(sum(l0 -
                                                                             lo)))
        
      }
    }
    SG$w[SG$uoCOUNT] = SG$w[SG$uoCOUNT] + 1
    
    
    # If on the first block
    if (pstar < 1.5) {
      SG$po[1:(SG$poCOUNT - 1),] = SG$po[2:SG$poCOUNT,]
      SG$pila[1:(SG$poCOUNT - 1),] = SG$pila[2:SG$poCOUNT,]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[2:SG$poCOUNT]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[2:SG$poCOUNT]
      I_mes[1:(SG$poCOUNT - 1)] = I_mes[2:SG$poCOUNT]
    }
    if (pstar > (SG$poCOUNT - 0.5)) {
      SG$po[1:(SG$poCOUNT - 1),] = SG$po[1:(pstar - 1),]
      SG$pila[1:(SG$poCOUNT - 1),] = SG$pila[1:(pstar - 1),]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[1:(pstar - 1)]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[1:(pstar - 1)]
      I_mes[1:(SG$poCOUNT - 1)] = I_mes[1:(pstar - 1)]
    }
    if (pstar < (SG$poCOUNT - 0.5) && pstar > 1.5) {
      SG$po[1:(SG$poCOUNT - 1),] = SG$po[c(1:(pstar - 1), (pstar + 1):SG$poCOUNT),]
      SG$pila[1:(SG$poCOUNT - 1),] = SG$pila[c(1:(pstar - 1), (pstar +1):SG$poCOUNT),]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[c(1:(pstar - 1), (pstar + 1):SG$poCOUNT)]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[c(1:(pstar - 1), (pstar + 1):SG$poCOUNT)]
      I_mes[1:(SG$poCOUNT - 1)] = I_mes[c(1:(pstar - 1), (pstar + 1):SG$poCOUNT)]
    }
    SG$poCOUNT = SG$poCOUNT - 1
    
    for (lcv2 in 1:SG$d) {
      lp = l0
      
      lp[lcv2] = lp[lcv2] + 1
      
      if (max(lp) < 7.5 && SG$poCOUNT < 4 * SG$ML) {
        kvals = which(lp > 1.5)
        
        canuse = 1
        ap = rep(0, SG$d)
        nap = 0
        for (lcv3 in 1:length(kvals)) {
          lpp = lp
          lpp[kvals[lcv3]] = lpp[kvals[lcv3]] - 1
          
          ismem = rep(1, SG$uoCOUNT)
          for (lcv4 in 1:SG$d) {
            ismem  = ismem * (SG$uo[1:SG$uoCOUNT, lcv4] == lpp[lcv4])
          }
          
          if (max(ismem) > 0.5) {
            ap[lcv3] = which(ismem > 0.5)
            nap = nap + 1
          } else{
            canuse = 0
          }
        }
        if (canuse > 0.5) {
          SG$poCOUNT = SG$poCOUNT + 1
          SG$po[SG$poCOUNT,] = lp
          SG$pogsize[SG$poCOUNT] = prod(SG$sizes[lp])
          SG$pila[SG$poCOUNT, 1:nap] = ap[1:nap]
          SG$pilaCOUNT[SG$poCOUNT] = nap
          I_mes[SG$poCOUNT] =  MSE_de(as.vector(SG$po[SG$poCOUNT, ]), MSE_v)
        }
      }
    }
  }
  
  SG$gridsizes = matrix(SG$sizes[SG$uo[1:SG$uoCOUNT, ]], SG$uoCOUNT, SG$d)
  SG$gridsizest = matrix(SG$sizest[SG$uo[1:SG$uoCOUNT, ]], SG$uoCOUNT, SG$d)
  SG$gridsize = apply(SG$gridsizes, 1, prod)
  SG$gridsizet = apply(SG$gridsizest, 1, prod)
  
  SG$di = matrix(0, nrow = SG$uoCOUNT, ncol = max(SG$gridsize))
  SG$dit = matrix(0, nrow = SG$uoCOUNT, ncol = sum((SG$gridsize)))
  
  SG$design = matrix(0, nrow = sum(SG$gridsize), ncol = SG$d)
  SG$designindex = matrix(0, nrow = sum(SG$gridsize), ncol = SG$d)
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
  
  return(SG)
}
