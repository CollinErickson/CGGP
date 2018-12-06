# Check that the predicted values make sense.
# Predict on design points. It should interpolate exactly,
#  so plotting actual vs pred should lie on y=x line.

# Use borehole
borehole <- function(x) {
  rw <- x[, 1] * (0.15 - 0.05) + 0.05
  r <-  x[, 2] * (50000 - 100) + 100
  Tu <- x[, 3] * (115600 - 63070) + 63070
  Hu <- x[, 4] * (1110 - 990) + 990
  Tl <- x[, 5] * (116 - 63.1) + 63.1
  Hl <- x[, 6] * (820 - 700) + 700
  L <-  x[, 7] * (1680 - 1120) + 1120
  Kw <- x[, 8] * (12045 - 9855) + 9855
  
  m1 <- 2 * pi * Tu * (Hu - Hl)
  m2 <- log(r / rw)
  m3 <- 1 + 2 * L * Tu / (m2 * rw ^ 2 * Kw) + Tu / Tl
  return(m1 / m2 / m3)
}
d = 8
testf<-function (x) {  return(borehole(x))} 

Ns <- 10001
par(mfrow=c(2,2))
#Ns <- c(1001,2001,4001,8001)
for (N in Ns) {
  SG = SGcreate(rep(0, d), rep(1, d),N, nugget=1e-3) #create the design.  it has so many entries because i am sloppy
  Y = testf(SG$design) #the design is $design, simple enough, right?
  logtheta <- rep(0,8)
  logtheta <- logthetaMLE(SG=SG,y=Y)
  inds <- sample(1:nrow(SG$design), 200, replace = F)
  plot(Y[inds], SGGPpred(xp=SG$design[inds,], SG=SG, y=Y, logtheta=logtheta)$me,
       xlab="Actual", ylab='Predicted', main=paste("Predictions of 200 design points, N =",N)); abline(a=0,b=1,col=2)
  
}
