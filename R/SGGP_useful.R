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
                " Cauchy, Gaussian, PowerExponential, Matern32, or Matern52"))
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