#' Calculate correlation matrix for two sets of points in one dimension
#' 
#' Note that this is not the correlation between two vectors.
#'
#' @param x1 Vector of coordinates from same dimension
#' @param x2 Vector of coordinates from same dimension
#' @param ... Don't use, just forces theta to be named
# ' @param LS Log of parameter that controls lengthscale
# ' @param FD Logit of 0.5*parameter  that controls the fractal demension
# ' @param HE Log of parameter that controls the hurst effect
#' @param theta Correlation parameter
#' @param return_dCdtheta Should dCdtheta be returned?
#' @param return_numpara Should it just return the number of parameters?
#'
#' @return Matrix
#' @export
#'
#' @examples
#' CorrMatMatern32(c(0,.2,.4),c(.1,.3,.5), theta=.1)
SGGP_internal_CorrMatCauchy <- function(x1, x2,theta, ..., return_dCdtheta = FALSE, return_numpara =FALSE) {
  if(return_numpara){
    return(3)
  }else{ 
  diffmat =abs(outer(x1,x2,'-')); 
  
  expLS = exp(3*(theta[1]))
  expHE = exp(3*(theta[2]))
  h = diffmat/expLS
  alpha = 2*exp(3*theta[3]+4)/(1+exp(3*theta[3]+4))
  halpha = h^alpha
  pow = -expHE/alpha
  
  C = (1+halpha)^pow
  if(return_dCdtheta){
    dCdtheta = cbind(3*expHE*((1+halpha)^(pow-1))*(halpha),dCdHE =3*C*pow*log(1+halpha), 3*C*(log(halpha+1)/alpha-halpha*log(h)/(halpha+1))*(expHE/(exp(4*theta[3]+4)+1)))
    dCdtheta[is.na(dCdtheta)] = 0
    out <- list(C=C,dCdtheta=dCdtheta)
    return(out)
  }else{
    return(C)
  }
  }
}

#' Calculate correlation matrix for two sets of points in one dimension
#' 
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams SGGP_internal_CorrMatCauchy
#'
#' @return Matrix
#' @export
#'
#' @examples
#' CorrMatMatern32(c(0,.2,.4),c(.1,.3,.5), theta=.1)
SGGP_internal_CorrMatCauchySQT <- function(x1, x2,theta, ..., return_dCdtheta = FALSE, return_numpara =FALSE) {
  if(return_numpara){
    return(3);
  } else{
    expTILT = exp((theta[3]))
    expLS = exp(3*(theta[1]))
    x1t = (x1+10^(-2))^expTILT
    x2t = (x2+10^(-2))^expTILT
    x1ts = x1t/expLS
    x2ts = x2t/expLS
    
    diffmat =abs(outer(x1ts,x2ts,'-')); 
    expHE = exp(3*(theta[2]))
    h = diffmat
    alpha = 2*exp(5)/(1+exp(5))
    halpha = h^alpha
    pow = -expHE/alpha
    
    C = (1+halpha)^pow
    
    if(return_dCdtheta){
       Q = ((1+halpha)^(pow-1))
       gt1 = x1t*log(x1+10^(-2))
       gt2 = x2t*log(x2+10^(-2))
       lh =outer(gt1,gt2,'-')
       hnabs = outer(x1ts,x2ts,'-')
       LO = alpha*expTILT*(pow/expLS)*(abs(h)^(alpha-1)*lh*sign(hnabs))
      dCdtheta = cbind(3*expHE*((1+halpha)^(pow-1))*(halpha),3*C*pow*log(1+halpha),LO*Q)
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}



#' Calculate correlation matrix for two sets of points in one dimension
#' 
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams SGGP_internal_CorrMatCauchy
#'
#' @return Matrix
#' @export
#'
#' @examples
#' CorrMatMatern32(c(0,.2,.4),c(.1,.3,.5), theta=.1)
SGGP_internal_CorrMatCauchySQ <- function(x1, x2,theta, ..., return_dCdtheta = FALSE, return_numpara =FALSE) {
  if(return_numpara){
    return(2);
  }else{ 
    diffmat =abs(outer(x1,x2,'-')); 
    
    expLS = exp(3*theta[1])
    expHE = exp(3*theta[2])
    h = diffmat/expLS
    alpha = 2*exp(0+6)/(1+exp(0+6))
    halpha = h^alpha
    pow = -expHE/alpha
    
    C = (1-10^(-10))*(1+halpha)^pow+10^(-10)*(diffmat<10^(-4))
    if(return_dCdtheta){
      dCdtheta = (1-10^(-10))*cbind(3*expHE*((1+halpha)^(pow-1))*(halpha),3*C*pow*log(1+halpha))
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}


