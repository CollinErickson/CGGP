# TODO add nugget as param, add to diag(S).
MSE_calc2 <- function(xl, theta, check_calcs=FALSE, methodA=FALSE) {
  xl = sort(xl)
  if (any(xl>1) || any(xl<0)) {stop("Can't be outside of 0,1")}
  S = CorrMat(xl, xl, theta)
  t = exp(theta)
  n = length(xl)
  Ci = solve(S)
  
  corrf <- function(a,b,t) {(1+abs(a-b)/t) * exp(-abs(a-b)/t)}
  # g1 <- function(cc,dd) {
  #   #exp(-(2*cc+dd)/t) * (-2*cc^2-2*cc*(dd+3*t)+t*(exp(2*cc/t)-1) * (3*dd+5*t)) / (4*t)
  #   exp(-(2*cc+dd)/t) * (t*(exp(2*cc/t)-1)*(3*dd+5*t)-2*cc*dd-3*cc*t*(exp(2*cc/t)+1)) / (4*t)
  #}
  if (methodA) {
    g2 <- function(cc,dd) {
      cc*exp(-dd/t)*(-2*cc^2+3*cc*dd+6*t*(dd+t)) / (6*t^2)
    }
    g3 <- function(cc,dd) {
      exp(-(2*cc+dd)/t) * (-2*cc^2-2*cc*(dd+3*t)+t*(exp(2*cc/t)-1) * (3*dd+5*t)) / (4*t)
    }
    
    t1 <- outer(1:length(xl), 1:length(xl), 
                Vectorize(function(i,j) {#browser()
                  a <- min(xl[i],xl[j])
                  b <- max(xl[i],xl[j])
                  #print(c(i,j,a,b))
                  #print(c(g3(a,b-a) , g2(b-a,b-a) , g3(1-b,b-a)))
                  #print(integrate(function(x) {corrf(x,a,t)*corrf(x,b,t)}, 0,1))
                  #print(c(integrate(function(x) {corrf(x,a,t)*corrf(x,b,t)}, 0,a)$value,
                  #        integrate(function(x) {corrf(x,a,t)*corrf(x,b,t)}, a,b)$value,
                  #        integrate(function(x) {corrf(x,a,t)*corrf(x,b,t)}, b,1)$value))
                  g3(a,b-a) + g2(b-a,b-a) + g3(1-b,b-a)
                })
    )
    t2 <- 1-sum(Ci * t1)
    return(t2)
  }
  
  A = matrix(rep(xl, each = n), nrow = n)
  a = pmin(A, t(A))
  b = pmax(A, t(A))
  bminusa <- b-a
  oneminusb <- 1-b
  g1 <- exp(-(2*a+bminusa)/t) * (-2*a^2-2*a*(bminusa+3*t)+t*(exp(2*a/t)-1) * (3*bminusa+5*t)) / (4*t)
  g2 <- bminusa*exp(-bminusa/t)*(-2*bminusa^2+3*bminusa*bminusa+6*t*(bminusa+t)) / (6*t^2)
  g3 <- exp(-(2*oneminusb+bminusa)/t) * (-2*oneminusb^2-2*oneminusb*(bminusa+3*t)+t*(exp(2*oneminusb/t)-1) * (3*bminusa+5*t)) / (4*t)
  if (check_calcs) {
    print(g1) #,g2,g3)
    print(outer(1:n,1:n, Vectorize(function(i,j) {
      a <- a[i,j]
      b <- b[i,j]
      print(c(i,j,a,b))
      c(integrate(function(x) {corrf(x,a,t)*corrf(x,b,t)}, 0,a)$value)})))
  }
  t1 <- g1 + g2 + g3
  t2 <- 1-sum(Ci * t1)
  return(t2)
  
  A = matrix(rep(xl, each = n), nrow = n)
  a = pmin(A, t(A))
  b = pmax(A, t(A))
  
  t2 = 1.0 / t
  t3 = a + b - 2.0
  t4 = t2 * t3
  t5 = exp(t4)
  t6 = a - b
  t7 = t2 * t6
  t8 = exp(t7)
  t9 = t ^ 2
  t10 = a * t2 * 2.0
  t11 = exp(t10)
  t12 = a * b * 2.0
  
  out1 = t5 * (-3.0 / 2.0) + a * t5 * (3.0 / 4.0) - a * t8 * (3.0 / 4.0) +
    b * t5 * (3.0 / 4.0) + b * t8 * (3.0 / 4.0) - t * t5 * (5.0 / 4.0) - t2 *
    t5 * (1.0 / 2.0) + t * t8 * (5.0 / 4.0) + a * t2 * t5 * (1.0 / 2.0) + b *
    t2 * t5 * (1.0 / 2.0) - t2 * exp(-t2 * (a + b)) * (t9 * 5.0 + t12 + a *
                                                         t * 3.0 + b * t * 3.0 - t9 * t11 * 5.0 + a * t * t11 * 3.0 - b * t * t11 *
                                                         3.0) * (1.0 / 4.0) - a * b * t2 * t5 * (1.0 / 2.0) - 1.0 / t ^ 2 * t6 *
    t8 * (t9 * 6.0 - t12 - a * t * 6.0 + b * t * 6.0 + a ^ 2 + b ^ 2) * (1.0 /
                                                                           6.0)
  
  out1
  # TODO Nugget should replace 1 below???
  MSE = 1 - sum(diag(out1 %*% Ci))
}