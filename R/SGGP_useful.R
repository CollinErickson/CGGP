#' Set correlation function of SGGP object
#'
#' @param SGGP SGGP object
#' @param corr Correlation function
#'
#' @return SGGP object
#' @export
#'
#' @examples
#' obj <- SGGPcreate(3, 20, corr="matern52")
#' SGGP_internal_set_corr(obj, "gaussian")
SGGP_internal_set_corr <- function(SGGP, corr) {
  if (is.function(corr)) {
    SGGP$CorrMat <- corr
    SGGP$CorrName <- "UserDefined"
  } else if (tolower(corr) %in% c("cauchysqt")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatCauchySQT
    SGGP$CorrName <- "CauchySQT"
  } else if (tolower(corr) %in% c("cauchysq")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatCauchySQ
    SGGP$CorrName <- "CauchySQ"
  } else if (tolower(corr) %in% c("cauchy")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatCauchy
    SGGP$CorrName <- "Cauchy"
  } else if (tolower(corr) %in% c("gaussian", "gauss", "sqexp")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatGaussian
    SGGP$CorrName <- "Gaussian"
  } else if (tolower(corr) %in% c("powerexp", "pe", "powerexponential")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatPowerExp
    SGGP$CorrName <- "PowerExponential"
  } else if (tolower(corr) %in% c("matern32", "m32", "m3")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatMatern32
    SGGP$CorrName <- "Matern32"
  } else if (tolower(corr) %in% c("matern52", "m52", "m5")) {
    SGGP$CorrMat <- SGGP_internal_CorrMatMatern52
    SGGP$CorrName <- "Matern52"
  } else {
    stop(paste0("corr given to SGGPcreate should be one of CauchySQT, CauchySQ,", 
                " Cauchy, Gaussian, PowerExponential, Matern32, or Matern52.\n",
                "Given value was ", corr, "."))
  }
  
  # Fix related parameters stored with SGGP
  SGGP$numpara <- SGGP$CorrMat(return_numpara=TRUE)
  SGGP$thetaMAP <- rep(0,SGGP$d*SGGP$numpara)
  SGGP$numPostSamples <- 100
  SGGP$thetaPostSamples  <- matrix(2*rbeta(SGGP$d*SGGP$numpara*SGGP$numPostSamples,
                                           0.5, 0.5)-1,
                                   ncol=SGGP$numPostSamples )
  SGGP
}

SGGP_internal_addrows <- function(SGGP, numrowstoadd=20) {
  
  SGGP$uo <- rbind(SGGP$uo, matrix(0,numrowstoadd,ncol(SGGP$uo)))
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
  
  SGGP
}

SGGP_internal_getdesignfromSGGP <- function(SGGP) {
  
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
          SGGP$design[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(rep(x0, each =
                                                                                   prod(SGGP$gridsizes[blocklcv, (dimlcv + 1):SGGP$d])), prod(SGGP$gridsizes[blocklcv, 1:(dimlcv - 1)]))
          SGGP$designindex[(tv + 1):(tv + SGGP$gridsize[blocklcv]), dimlcv] = rep(rep(xi0, each =
                                                                                        prod(SGGP$gridsizes[blocklcv, (dimlcv + 1):SGGP$d])), prod(SGGP$gridsizes[blocklcv, 1:(dimlcv - 1)]))
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
    
    tv = tv + SGGP$gridsize[blocklcv]
  }
  SGGP
}