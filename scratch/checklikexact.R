CorrMat32Vecs <- function(u,v, theta) {
  d <- abs(u-v)
  prod((1+d/theta/sqrt(3)) * exp(-d/theta/sqrt(3)))
}
CorrMat32Full <- function(X, theta) {
  n <- nrow(X)
  outer(1:n, 1:n, Vectorize(function(i,j) {CorrMat32Vecs(X[i,],X[j,],theta)}))
}
# I feel like the product by dimension isn't equal to doing all dims at once.
debuglik <- T
SG = SGcreate(rep(0, d), rep(1, d),51) #create the design.  it has so many entries because i am sloppy
Y = testf(SG$design) #the design is $design, simple enough, right?
logtheta <- rep(0,8)
lik(logtheta = logtheta, SG = SG, y = Y)
#logthetaest = logthetaMLE(SG,Y)

lDet
log(det(CorrMat32Full(SG$design, theta=exp(logtheta))))

log(sigma_hat)
log(sum(y * solve(CorrMat32Full(SG$design, theta=exp(logtheta)), y))/length(y))
