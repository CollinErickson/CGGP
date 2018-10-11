SGcreate <- function(xmin, xmax,batchsize) {
  SG = list("xmin" = xmin, "xmax" = xmax)
  
  SG$d = length(xmin)
  SG$ML = min(choose(SG$d + 6, SG$d), 100000) #max levels
  
  SG$levelpoint = rep(0, SG$ML)
  
  SG$uo = matrix(0, nrow = SG$ML, ncol = SG$d) #used levels tracker
  SG$uo[1, ] = rep(1, SG$d) #first observation in middle of space
  SG$uoCOUNT = 1 #number of used levels
  
  SG$po = matrix(0, nrow = 4 * SG$ML, ncol = SG$d) #proposed levels tracker
  SG$po[1:SG$d, ] = matrix(1, nrow = SG$d, ncol = SG$d) + diag(SG$d) #one at a time
  SG$poCOUNT = SG$d #number of proposed levels
  
  
  
  SG$pila = matrix(0, nrow = SG$ML, ncol = 1000) #proposed immediate level ancestors
  SG$pala = matrix(0, nrow = SG$ML, ncol = 1000) #proposedal all level ancestors
  SG$uala = matrix(0, nrow = SG$ML, ncol = 1000) #used all level ancestors
  SG$pilaCOUNT = rep(0, SG$ML) #count of number of pila
  SG$palaCOUNT = rep(0, SG$ML) #count of number of pala
  SG$ualaCOUNT = rep(0, SG$ML) #count of number of uala
  
  SG$pilaCOUNT[1:SG$d] = 1
  SG$pila[1:SG$d, 1] = 1
  
  SG$bss = batchsize#1+4*SG$d  #must be at least 3*d
  SG$sizes = c(1, 2, 2, 2, 4, 4, 4, 6, 8)
  SG$pogsize = rep(0, 4 * SG$ML)
  SG$pogsize[1:SG$poCOUNT] = apply(matrix(SG$sizes[SG$po[1:SG$poCOUNT, ]], SG$poCOUNT, SG$d), 1, prod)
  SG$ss = 1
  
  
  SG$w = rep(0, SG$ML) #keep track of + and - for prediction
  SG$w[1] = 1 #keep track of + and - for prediction
  SG$uoCOUNT = 1
  while (SG$bss > (SG$ss + min(SG$pogsize[1:SG$poCOUNT]) - 0.5)) {
    SG$uoCOUNT = SG$uoCOUNT + 1 #increment used count
    if (SG$uoCOUNT < (SG$d + 1.5)) {
      pstar = 1 #pick a proposed to add
    } else{
      if (SG$uoCOUNT < (2 * SG$d + 1.5)) {
        pstar = sample(which(SG$pogsize[1:SG$poCOUNT] <= 0.5 + min(SG$pogsize[1:SG$poCOUNT])), 1)
      } else{
        pstar = sample(which(SG$pogsize[1:SG$poCOUNT] < (SG$bss - SG$ss + 0.5)), 1)
      }
    }
    
    l0 =  SG$po[pstar, ]
    SG$uo[SG$uoCOUNT, ] = l0
    SG$ss =  SG$ss + SG$pogsize[pstar]
    
    new_an = SG$pila[pstar, 1:SG$pilaCOUNT[pstar]]
    total_an = new_an
    # for(lcv5 in 1:10){
    for (lcv2 in 1:length(total_an)) {
      if (total_an[lcv2] > 1.5) {
        total_an = unique(c(total_an, SG$uala[total_an[lcv2], 1:SG$ualaCOUNT[total_an[lcv2]]]))
      }
      #}
    }
    SG$ualaCOUNT[SG$uoCOUNT]  = length(total_an)
    SG$uala[SG$uoCOUNT, 1:length(total_an)] = total_an
    
    for (lcv2 in 1:length(total_an)) {
      lo = SG$uo[total_an[lcv2], ]
      if (max(abs(lo - l0)) < 1.5) {
        SG$w[total_an[lcv2]] = SG$w[total_an[lcv2]] + (-1) ^ abs(round(sum(l0 -
                                                                             lo)))
        
      }
    }
    SG$w[SG$uoCOUNT] = SG$w[SG$uoCOUNT] + 1
    
    
    if (pstar < 1.5) {
      SG$po[1:(SG$poCOUNT - 1), ] = SG$po[2:SG$poCOUNT, ]
      SG$pila[1:(SG$poCOUNT - 1), ] = SG$pila[2:SG$poCOUNT, ]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[2:SG$poCOUNT]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[2:SG$poCOUNT]
    }
    if (pstar > (SG$poCOUNT - 0.5)) {
      SG$po[1:(SG$poCOUNT - 1), ] = SG$po[1:(pstar - 1), ]
      SG$pila[1:(SG$poCOUNT - 1), ] = SG$pila[1:(pstar - 1), ]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[1:(pstar - 1)]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[1:(pstar - 1)]
    }
    if (pstar < (SG$poCOUNT - 0.5) && pstar > 1.5) {
      SG$po[1:(SG$poCOUNT - 1), ] = SG$po[c(1:(pstar - 1), (pstar + 1):SG$poCOUNT), ]
      SG$pila[1:(SG$poCOUNT - 1), ] = SG$pila[c(1:(pstar - 1), (pstar +
                                                                  1):SG$poCOUNT), ]
      SG$pilaCOUNT[1:(SG$poCOUNT - 1)] = SG$pilaCOUNT[c(1:(pstar - 1), (pstar +
                                                                          1):SG$poCOUNT)]
      SG$pogsize[1:(SG$poCOUNT - 1)] = SG$pogsize[c(1:(pstar - 1), (pstar +
                                                                      1):SG$poCOUNT)]
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
          SG$po[SG$poCOUNT, ] = lp
          SG$pogsize[SG$poCOUNT] = prod(SG$sizes[lp])
          SG$pila[SG$poCOUNT, 1:nap] = ap[1:nap]
          SG$pilaCOUNT[SG$poCOUNT] = nap
          
        }
      }
    }
  }
  
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
  SG$sizest = cumsum(SG$sizes)
  
  
  SG$gridsizes = matrix(SG$sizes[SG$uo[1:SG$uoCOUNT, ]], SG$uoCOUNT, SG$d)
  SG$gridsizest = matrix(SG$sizest[SG$uo[1:SG$uoCOUNT, ]], SG$uoCOUNT, SG$d)
  SG$gridsize = apply(SG$gridsizes, 1, prod)
  SG$gridsizet = apply(SG$gridsizest, 1, prod)
  
  SG$di = matrix(0, nrow = SG$uoCOUNT, ncol = max(SG$gridsize))
  SG$dit = matrix(0, nrow = SG$uoCOUNT, ncol = sum((SG$gridsize)))
  
  SG$design = matrix(0, nrow = sum(SG$gridsize), ncol = SG$d)
  tv = 0
  for (lcv1 in 1:SG$uoCOUNT) {
    SG$di[lcv1, 1:SG$gridsize[lcv1]] = (tv + 1):(tv + SG$gridsize[lcv1])
    for (lcv2 in 1:SG$d) {
      levelnow = SG$uo[lcv1, lcv2]
      if (levelnow < 1.5) {
        SG$design[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(SG$xb[1], SG$gridsize[lcv1])
      } else{
        x0 = SG$xb[(SG$sizest[levelnow - 1] + 1):SG$sizest[levelnow]]
        if (lcv2 < 1.5) {
          SG$design[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(x0, "each" = SG$gridsize[lcv1] /
                                                                     SG$gridsizes[lcv1, lcv2])
        }
        if (lcv2 > (SG$d - 0.5)) {
          SG$design[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(x0, SG$gridsize[lcv1] /
                                                                     SG$gridsizes[lcv1, lcv2])
        }
        if (lcv2 < (SG$d - 0.5)  && lcv2 > 1.5) {
          SG$design[(tv + 1):(tv + SG$gridsize[lcv1]), lcv2] = rep(rep(x0, "each" =
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