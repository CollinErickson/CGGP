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
CorrMatCauchy <- function(x1, x2,theta, ..., return_dCdtheta = FALSE) {
  diffmat =abs(outer(x1,x2,'-')); 
  
  expLS = exp(theta[1])
  expHE = exp(theta[2])
  h = diffmat/expLS
  theta3=3
  alpha = 2*exp(theta3)/(1+exp(theta3))
  halpha = h^alpha
  pow = -expHE/alpha
  
  C = (1+halpha)^pow
  if(return_dCdtheta){
    dCdtheta = cbind(expHE*((1+halpha)^(pow-1))*(halpha),dCdHE = C*pow*log(1+halpha))
    dCdtheta[is.na(dCdtheta)] = 0
    out <- list(C=C,dCdtheta=dCdtheta)
    return(out)
  }else{
    return(C)
  }
  
  
}


#' #' Calculate correlation matrix for two sets of points in one dimension
#' #' 
#' #' Note that this is not the correlation between two vectors.
#' #' It is for two
#' #'
#' #' @param x1 Vector of coordinates from same dimension
#' #' @param x2 Vector of coordinates from same dimension
#' #' @param ... Don't use, just forces theta to be named
#' #' @param LS Log of parameter that controls lengthscale
#' #' @param FD Logit of 0.5*parameter  that controls the fractal demension
#' #' @param HE Log of parameter that controls the hurst effect
#' #'
#' #' @return Matrix
#' #' @export
#' #'
#' #' @examples
#' #' CorrMatMatern32(c(0,.2,.4),c(.1,.3,.5), theta=.1)
#' CorrMatCauchy <- function(x1, x2,theta, ..., return_dCdtheta = FALSE) {
#'   diffmat =abs(outer(x1,x2,'-')); 
#'   
#'   expLS = exp(theta[1])
#'   expHE = exp(theta[2])
#'   h = diffmat/expLS
#'   alpha = 2*exp(theta[3])/(1+exp(theta[3]))
#'   halpha = h^alpha
#'   pow = -expHE/alpha
#'   
#'   C = (1+halpha)^pow
#'   if(return_dCdtheta){
#'     dCdtheta = cbind(expHE*((1+halpha)^(pow-1))*(halpha),dCdHE = C*pow*log(1+halpha), C*(log(halpha+1)/alpha-halpha*log(h)/(halpha+1))*(expHE/(exp(theta[3])+1)))
#'     dCdtheta[is.na(dCdtheta)] = 0
#'     out <- list(C=C,dCdtheta=dCdtheta)
#'     return(out)
#'   }else{
#'     return(C)
#'   }
#'   
#'   
#' }