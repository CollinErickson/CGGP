#' Calculate correlation matrix for two sets of points in one dimension
#' 
#' Note that this is not the correlation between two vectors.
#' It is for two
#'
#' @param x1 Vector of coordinates from same dimension
#' @param x2 Vector of coordinates from same dimension
#' @param ... Don't use, just forces theta to be named
#' @param LS Log of parameter that controls lengthscale
#' @param FD Logit of 0.5*parameter  that controls the fractal demension
#' @param HE Log of parameter that controls the hurst effect
#'
#' @return Matrix
#' @export
#'
#' @examples
#' CorrMatMatern32(c(0,.2,.4),c(.1,.3,.5), theta=.1)
CorrMatCauchy <- function(x1, x2,..., LS =0, FD = 0, HE = 0, 
                            return_C = TRUE, return_diagC = FALSE, 
                            return_dCdLS = FALSE, return_dCdFD = FALSE, return_dCdHE = FALSE) {
  d1 = length(x1)
  d2 = length(x2)
  # Would be easier to use outer
  diffmat =abs(outer(x1,x2,'-')); 
  
  expLS = exp(LS)
  expHE = exp(HE)
  h = diffmat/expLS
  alpha = 2*exp(FD)/(1+exp(FD))
  halpha = h^alpha
  pow = -expHE/alpha
  
  C = (1+halpha)^pow
  if(return_dCdLS || return_dCdFD || return_dCdHE){
  if(return_dCdLS){
    dCdLS = -expHE*(pow*(1+halpha)^(pow-1))*(halpha)
  } else{
    dCdLS = NA
  }
  if(return_dCdFD){
    dCdFD = C*(log(halpha+1)/alpha-halpha*log(h)/(halpha+1))*(expHE/(exp(FD)+1))
  } else{
    dCdFD = NA
  }
  if(return_dCdHE){
    dCdHE = C*pow*log(1+halpha)
  } else{
    dCdHE = NA
  }
    
    out <- list(C=C,
                dCdLS=dCdLS,
                dCdFD=dCdFD,
                dCdHE=dCdHE)
    return(out)
  }else{
    return(C)
  }
  
  
}
