#' Calculate MSE over single dimension
#' 
#' Calcaluted using grid of integration points.
#' Can be calculated exactly, but not much reason in 1D.
#'
#' @param xl Vector of points in 1D
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#'
#' @return MSE value
#' @export
#'
#' @examples
#' SGGP_internal_calcMSE(xl=c(0,.5,.9), theta=c(1,2,3),
#'          CorrMat=SGGP_internal_CorrMatCauchySQT)
SGGP_internal_calcMSE <- function(xl, theta, CorrMat) {
  S = CorrMat(xl, xl, theta)
  xp = seq(0,1,l=101)
  Cp = CorrMat(xp,xl,theta)
  n = length(xl)
  cholS = chol(S)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  
  MSE_MAPal = mean(1 - rowSums(t(CiCp)*Cp))
  
  MSE_MAPal
}


#' Calculate MSE over blocks
#' 
#' Delta of adding block is product over i=1..d of IMSE(i,j-1) - IMSE(i,j)
#'
#' @param valsinds Block levels to calculate MSEs for
#' @param MSE_MAP Matrix of MSE values
#'
#' @return All MSE values
#' @export
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG, Y=y)
#' MSE_MAP <- outer(1:SG$d, 1:8, 
#'      Vectorize(function(dimlcv, lcv1) {
#'         SGGP_internal_calcMSE(SG$xb[1:SG$sizest[dimlcv]],
#'         theta=SG$thetaMAP[(dimlcv-1)*SG$numpara+1:SG$numpara],
#'         CorrMat=SG$CorrMat)
#'  }))
#' SGGP_internal_calcMSEde(SG$po[1:SG$poCOUNT, ], MSE_MAP)
SGGP_internal_calcMSEde <- function(valsinds, MSE_MAP) {
  if(is.matrix(valsinds)){
    MSE_de = rep(0, dim(valsinds)[1])
    
    for (levellcv2 in 1:dim(valsinds)[1]) {
      MSE_de[levellcv2] = 0
      for (levellcv in 1:dim(valsinds)[2]) {
        if (valsinds[levellcv2, levellcv] > 1.5) {
          MSE_de[levellcv2] = MSE_de[levellcv2] + max(log(-MSE_MAP[levellcv, valsinds[levellcv2, levellcv]] + 
                                                        MSE_MAP[levellcv, valsinds[levellcv2, levellcv] - 1]),-10)
          
        } else {
          # This is when no ancestor block, 1 comes from when there is no data. 
          # 1 is correlation times integrated value over range.
          # This depends on correlation function.
          MSE_de[levellcv2] = MSE_de[levellcv2] + max(log(-MSE_MAP[levellcv, valsinds[levellcv2, levellcv]] + 1),-10)
          
        }
      }
    }
  } else {
    MSE_de = 0
    
    for (levellcv in 1:length(valsinds)) {
      if (valsinds[levellcv] > 1.5) {
        MSE_de = MSE_de + max(log(-MSE_MAP[levellcv, valsinds[levellcv]] + MSE_MAP[levellcv, valsinds[levellcv] -1]),-10)
        
      } else {
        MSE_de = MSE_de + max(log(-MSE_MAP[levellcv, valsinds[levellcv]] + 1),-10)
        
      }
    }
  }
  
  MSE_de = exp(MSE_de)
  
  return(MSE_de)
}



#' Add points to SGGP
#' 
#' Add `batchsize` points to `SG` using `theta`.
#'
#' @param SGGP Sparse grid object
#' @param batchsize Number of points to add
#' @param selectionmethod How points will be selected: one of `UCB`, `TS`, or `Greedy`
#' @param RIMSEperpoint Should RIMSE per point be used?
#' @param multioutputdim_weights Weights for each output dimension.
#' @importFrom stats quantile sd var
#'
#' @return SG with new points added.
#' @export
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- SGGPfit(SG, Y=y)
#' SG <- SGGPappend(SGGP=SG, batchsize=20)
#' # UCB,TS,Greedy
SGGPappend <- function(SGGP,batchsize, selectionmethod = "UCB", RIMSEperpoint=TRUE,
                       multioutputdim_weights=1){
  if (!(selectionmethod %in% c("UCB", "TS", "Greedy"))) {
    stop("selectionmethod in SGGPappend must be one of UCB, TS, or Greedy")
  }
  if (is.numeric(multioutputdim_weights)) {
    if (length(multioutputdim_weights) != 1 && length(multioutputdim_weights) != ncol(SGGP$y)) {
      stop("multioutputdim_weights if numeric must have length 1 or number of outputs")
    }
  } else if (multioutputdim_weights == "sd") {
    multioutputdim_weights <- apply(SGGP$Y, 2, sd)
  } else if (multioutputdim_weights == "var") {
    multioutputdim_weights <- apply(SGGP$Y, 2, var)
  } else {
    stop("multioutputdim_weights not acceptable")
  }
  
  n_before <- nrow(SGGP$design)
  
  max_polevels = apply(SGGP$po[1:SGGP$poCOUNT,], 2, max)
  
  separateoutputparameterdimensions <- is.matrix(SGGP$thetaMAP)
  # nnn is numberofoutputparameterdimensions
  nnn <- if (separateoutputparameterdimensions) {
    ncol(SGGP$y)
  } else {
    1
  }
  
  if(selectionmethod=="Greedy"){
    # Set up blank matrix to store MSE values
    # MSE_MAP = matrix(0, SGGP$d, SGGP$maxlevel)
    # Now use an array for nnn
    MSE_MAP = array(0, dim=c(SGGP$d, SGGP$maxlevel,nnn))
    
    # Loop over dimensions and design refinements
    for (opdlcv in 1:nnn) {
      thetaMAP.thisloop <- if (nnn==1) SGGP$thetaMAP else SGGP$thetaMAP[, opdlcv]
      for (dimlcv in 1:SGGP$d) {
        for (levellcv in 1:max_polevels[dimlcv]) {
          # Calculate some sort of MSE from above, not sure what it's doing
          MSE_MAP[dimlcv, levellcv, opdlcv] = max(0, 
                                                  abs(SGGP_internal_calcMSE(SGGP$xb[1:SGGP$sizest[levellcv]],
                                                                            thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
                                                                            SGGP$CorrMat)))
          if (levellcv > 1.5) { # If past first level, it is as good as one below it. Why isn't this a result of calculation?
            MSE_MAP[dimlcv, levellcv, opdlcv] = min(MSE_MAP[dimlcv, levellcv, opdlcv], MSE_MAP[dimlcv, levellcv - 1, opdlcv])
          }
        }
      }
    }
    
    # What is this? Integrate MSE
    IMES_MAP = rep(0, SGGP$ML)
    
    # For all possible blocks, calculate MSE_MAP? Is that all that MSE_de does?
    # IMES_MAP[1:SGGP$poCOUNT] = SGGP_internal_calcMSEde(SGGP$po[1:SGGP$poCOUNT, ], MSE_MAP)
    # Need to apply it over nnn
    IMES_MAP_beforemean = apply(MSE_MAP, 3, function(x) SGGP_internal_calcMSEde(SGGP$po[1:SGGP$poCOUNT, ], x))
    # IMES_MAP_apply has a column for each opdlcv, so take rowMeans
    IMES_MAP[1:SGGP$poCOUNT] = rowMeans(IMES_MAP_beforemean)
    
    # Clean up to avoid silly errors
    rm(opdlcv, thetaMAP.thisloop)
  } else { # selectionmethod is UCB or TS
    # MSE_PostSamples = array(0, c(SGGP$d, SGGP$maxlevel,SGGP$numPostSamples))
    # Array needs another dimension for nnn
    MSE_PostSamples = array(0, c(SGGP$d, SGGP$maxlevel,SGGP$numPostSamples, nnn))
    #  MSE_UCB = matrix(0, SGGP$d, SGGP$maxlevel)
    # Dimensions can be considered independently
    # Loop over dimensions and design refinements
    for (opdlcv in 1:nnn) { # Loop over output parameter dimensions
      thetaPostSamples.thisloop <- if (nnn==1) SGGP$thetaPostSamples else SGGP$thetaPostSamples[, , opdlcv]
      for (dimlcv in 1:SGGP$d) {
        for (levellcv in 1:max_polevels[dimlcv]) {
          for(samplelcv in 1:SGGP$numPostSamples){
            # Calculate some sort of MSE from above, not sure what it's doing
            MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv] = 
              max(0, 
                  abs(
                    SGGP_internal_calcMSE(
                      SGGP$xb[1:SGGP$sizest[levellcv]],
                      thetaPostSamples.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara,
                                                samplelcv],
                      SGGP$CorrMat)
                  )
              )
            if (levellcv > 1.5) { # If past first level, it is as good as one below it. Why isn't this a result of calculation?
              MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv] = min(MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv],
                                                                        MSE_PostSamples[dimlcv, levellcv - 1,samplelcv, opdlcv])
            }
          }
          #  done below MSE_UCB[dimlcv, levellcv] = quantile(MSE_PostSamples[dimlcv, levellcv,],0.99)
        }
      }
    }
    IMES_PostSamples = matrix(0, SGGP$ML,SGGP$numPostSamples)
    for(samplelcv in 1:SGGP$numPostSamples){
      # IMES_PostSamples[1:SGGP$poCOUNT,samplelcv] = SGGP_internal_calcMSEde(SGGP$po[1:SGGP$poCOUNT,], MSE_PostSamples[,,samplelcv])
      # It's an array, need to average over opdlcv. Over 3rd dim since samplelcv removes 3rd dim of array.
      if (nnn == 1) { # Will be a matrix
        IMES_PostSamples[1:SGGP$poCOUNT,samplelcv] = SGGP_internal_calcMSEde(SGGP$po[1:SGGP$poCOUNT,], MSE_PostSamples[,,samplelcv,])
      } else { # Is a 3d array, need to use an apply and then apply again with mean
        IMES_PostSamples_beforemean <- apply(MSE_PostSamples[,,samplelcv,], 3, function(x){SGGP_internal_calcMSEde(SGGP$po[1:SGGP$poCOUNT,], x)})
        IMES_PostSamples[1:SGGP$poCOUNT,samplelcv] <- apply(IMES_PostSamples_beforemean, 1,
                                                            function(x) {
                                                              mean(multioutputdim_weights*x)
                                                            })
      }
    }
    IMES_UCB = matrix(0, SGGP$ML)
    IMES_UCB[1:SGGP$poCOUNT] = apply(IMES_PostSamples[1:SGGP$poCOUNT,],1,quantile, probs=0.9) 
    # Clean up to avoid silly errors
    rm(opdlcv, thetaPostSamples.thisloop)
  }
  
  
  # Increase count of points evaluated. Do we check this if not reached exactly???
  SGGP$bss = SGGP$bss + batchsize
  
  # Keep adding points until reaching bss
  while (SGGP$bss > (SGGP$ss + min(SGGP$pogsize[1:SGGP$poCOUNT]) - 0.5)) {
    
    if(selectionmethod=="Greedy"){
      IMES = IMES_MAP
    } else if(selectionmethod=="UCB"){
      IMES = IMES_UCB
    } else if(selectionmethod=="TS"){
      IMES = IMES_PostSamples[,sample(1:SGGP$numPostSamples,1)]
    } else {
      stop("Selection method not acceptable")
    }
    SGGP$uoCOUNT = SGGP$uoCOUNT + 1 #increment used count
    
    
    # Old way, no RIMSEperpoint option
    # # Find the best one that still fits
    # M_comp = max(IMES[which(SGGP$pogsize[1:SGGP$poCOUNT] < (SGGP$bss - SGGP$ss + 0.5))])
    # # Find which ones are close to M_comp and
    # possibleO =which((IMES[1:SGGP$poCOUNT] >= 0.99*M_comp)&(SGGP$pogsize[1:SGGP$poCOUNT] < (SGGP$bss - SGGP$ss + 0.5)))
    
    # New way, now you can pick best IMES per point in the block, more efficient
    stillpossible <- which(SGGP$pogsize[1:SGGP$poCOUNT] < (SGGP$bss - SGGP$ss + 0.5))
    
    # Either pick block with max IMES or with max IMES per point in the block.
    if (RIMSEperpoint) {
      metric <- IMES[1:SGGP$poCOUNT] / SGGP$pogsize[1:SGGP$poCOUNT]
    } else {
      metric <- IMES[1:SGGP$poCOUNT]
    }
    # Find the best one that still fits
    M_comp = max(metric[stillpossible])
    # Find which ones are close to M_comp and
    # possibleO =which((IMES[stillpossible] >= 0.99*M_comp)&(SGGP$pogsize[1:SGGP$poCOUNT] < (SGGP$bss - SGGP$ss + 0.5)))
    possibleO = stillpossible[metric[stillpossible] >= 0.99*M_comp]
    
    
    # If more than one is possible and near the best, randomly pick among them.
    if(length(possibleO)>1.5){
      pstar = sample(possibleO,1)
    } else{
      pstar = possibleO
    }
    
    l0 =  SGGP$po[pstar,] # Selected block
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
    # print(list(dim(SGGP$uo), SGGP$uoCOUNT, SGGP$uo[SGGP$uoCOUNT,], l0))
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
        SGGP$w[total_an[anlcv]] = SGGP$w[total_an[anlcv]] + (-1)^abs(round(sum(l0-lo)))
      }
    }
    SGGP$w[SGGP$uoCOUNT] = SGGP$w[SGGP$uoCOUNT] + 1
    
    
    # Update data. Remove selected item, move rest up.
    # First get correct indices to change
    new_indices <- 1:(SGGP$poCOUNT - 1)
    if (pstar < 1.5) {
      old_indices <- 2:SGGP$poCOUNT
    } else if (pstar > (SGGP$poCOUNT - 0.5)) {
      old_indices <- 1:(pstar - 1)
    } else if (pstar < (SGGP$poCOUNT - 0.5) && pstar > 1.5) {
      old_indices <- c(1:(pstar - 1), (pstar + 1):SGGP$poCOUNT)
    } else {stop("Not possible #729588")}
    # Then change the data
    SGGP$po[new_indices,] = SGGP$po[old_indices,]
    SGGP$pila[new_indices,] = SGGP$pila[old_indices,]
    SGGP$pilaCOUNT[new_indices] = SGGP$pilaCOUNT[old_indices]
    SGGP$pogsize[new_indices] = SGGP$pogsize[old_indices]
    if(selectionmethod=="Greedy"){
      IMES_MAP[new_indices] = IMES_MAP[old_indices]
    }
    if(selectionmethod=="UCB"){
      IMES_UCB[new_indices] = IMES_UCB[old_indices]
    }
    if(selectionmethod=="TS"){
      IMES_PostSamples[new_indices,] = IMES_PostSamples[old_indices,]
    }
    # And reduce number of available blocks by one.
    SGGP$poCOUNT = SGGP$poCOUNT - 1
    
    # Loop over possible descendents of selected block, add them if possible    
    for (dimlcv in 1:SGGP$d) {
      lp = l0
      
      lp[dimlcv] = lp[dimlcv] + 1
      
      if (max(lp) < SGGP$maxlevel && SGGP$poCOUNT < 4 * SGGP$ML) {
        kvals = which(lp > 1.5) # Dimensions above base level
        
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
        if (canuse > 0.5) { # If it can be used, add to possible blocks
          SGGP$poCOUNT = SGGP$poCOUNT + 1
          SGGP$po[SGGP$poCOUNT,] = lp
          SGGP$pogsize[SGGP$poCOUNT] = prod(SGGP$sizes[lp])
          SGGP$pila[SGGP$poCOUNT, 1:nap] = ap[1:nap]
          SGGP$pilaCOUNT[SGGP$poCOUNT] = nap
          
          max_polevels_old = max_polevels
          max_polevels = apply(SGGP$po[1:SGGP$poCOUNT,], 2, max)
          
          if(selectionmethod=="Greedy"){
            for (opdlcv in 1:nnn) { # Loop over output parameter dimensions
              thetaMAP.thisloop <- if (nnn==1) SGGP$thetaMAP else SGGP$thetaMAP[, opdlcv]
              for (dimlcv in 1:SGGP$d) {
                if((max_polevels_old[dimlcv]+0.5)<max_polevels[dimlcv]){
                  levellcv = max_polevels[dimlcv]
                  MSE_MAP[dimlcv, levellcv, opdlcv] = max(0, abs(SGGP_internal_calcMSE(SGGP$xb[1:SGGP$sizest[levellcv]],
                                                                                       thetaMAP.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara],
                                                                                       SGGP$CorrMat)))
                  if (levellcv > 1.5) { # If past first level, it is as good as one below it. Why isn't this a result of calculation?
                    MSE_MAP[dimlcv, levellcv, opdlcv] = min(MSE_MAP[dimlcv, levellcv, opdlcv], MSE_MAP[dimlcv, levellcv - 1, opdlcv])
                  }
                }
              }
            }
            # Clean up
            rm(thetaMAP.thisloop, opdlcv)
          }else{
            for (opdlcv in 1:nnn) {
              thetaPostSamples.thisloop <- if (nnn==1) SGGP$thetaPostSamples else SGGP$thetaPostSamples[, , opdlcv]
              for (dimlcv in 1:SGGP$d) {
                if((max_polevels_old[dimlcv]+0.5)<max_polevels[dimlcv]){
                  levellcv = max_polevels[dimlcv]
                  for(samplelcv in 1:SGGP$numPostSamples){
                    # Calculate some sort of MSE from above, not sure what it's doing
                    MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv] = max(0,
                                                                              abs(SGGP_internal_calcMSE(
                                                                                SGGP$xb[1:SGGP$sizest[levellcv]],
                                                                                thetaPostSamples.thisloop[(dimlcv-1)*SGGP$numpara+1:SGGP$numpara,
                                                                                                          samplelcv],
                                                                                SGGP$CorrMat)))
                    if (levellcv > 1.5) { # If past first level, it is as good as one below it. Why isn't this a result of calculation?
                      MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv] = min(MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv],
                                                                                MSE_PostSamples[dimlcv, levellcv - 1,samplelcv, opdlcv])
                    }
                  }
                }
              }
            }
            # Clean up
            rm(thetaPostSamples.thisloop, opdlcv)
          }
          
          if(selectionmethod=="Greedy"){
            # IMES_MAP[SGGP$poCOUNT] = SGGP_internal_calcMSEde(as.vector(SGGP$po[SGGP$poCOUNT, ]), MSE_MAP)
            # Need to apply first
            IMES_MAP_beforemeannewpoint <- apply(MSE_MAP, 3, function(x) {SGGP_internal_calcMSEde(as.vector(SGGP$po[SGGP$poCOUNT, ]), x)})
            IMES_MAP[SGGP$poCOUNT] <- mean(IMES_MAP_beforemeannewpoint)
          } else if (selectionmethod=="UCB" || selectionmethod=="TS"){
            for(samplelcv in 1:SGGP$numPostSamples){
              # IMES_PostSamples[SGGP$poCOUNT,samplelcv] = SGGP_internal_calcMSEde(as.vector(SGGP$po[SGGP$poCOUNT, ]),
              #                                                                    MSE_PostSamples[,,samplelcv])
              if (nnn == 1) { # is a matrix
                IMES_PostSamples[SGGP$poCOUNT,samplelcv] = SGGP_internal_calcMSEde(as.vector(SGGP$po[SGGP$poCOUNT, ]),
                                                                                   MSE_PostSamples[,,samplelcv,])
              } else { # is an array, need to apply
                IMES_PostSamples_beforemeannewpoint = apply(MSE_PostSamples[,,samplelcv,],
                                                            3, # 3rd dim since samplelcv removes 3rd
                                                            function(x) {SGGP_internal_calcMSEde(as.vector(SGGP$po[SGGP$poCOUNT, ]), x)}
                )
                IMES_PostSamples[SGGP$poCOUNT,samplelcv] <- mean(multioutputdim_weights * IMES_PostSamples_beforemeannewpoint)
              }
            }
            IMES_UCB[SGGP$poCOUNT] = quantile(IMES_PostSamples[SGGP$poCOUNT,],probs=0.9)
          } else {stop("Not possible #9235058")}
          # Removing below and adding it as part of UCB above
          # if(selectionmethod=="TS"){
          #   for(samplelcv in 1:SGGP$numPostSamples){
          #     IMES_PostSamples[SGGP$poCOUNT,samplelcv] = SGGP_internal_calcMSEde(as.vector(SGGP$po[SGGP$poCOUNT, ]),
          #                                                                        MSE_PostSamples[,,samplelcv])
          #   }
          # }
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
          SGGP$design[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] =
            rep(
              rep(x0,
                  each=prod(SGGP$gridsizes[blocklcv, (dimlcv + 1):SGGP$d])
              ),
              prod(SGGP$gridsizes[blocklcv, 1:(dimlcv - 1)])
            )
          SGGP$designindex[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] =
            rep(
              rep(xi0,
                  each=prod(SGGP$gridsizes[blocklcv, (dimlcv + 1):SGGP$d])
              ),
              prod(SGGP$gridsizes[blocklcv, 1:(dimlcv - 1)])
            )
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
  
  # Check if none were added, return warning/error
  if (n_before == nrow(SGGP$design)) {
    warning("No points could be added. You may need a larger batch size.")
  } else {
    # Save design_unevaluated to make it easy to know which ones to add
    SGGP$design_unevaluated <- SGGP$design[(n_before+1):nrow(SGGP$design),]
  }
  return(SGGP)
}
