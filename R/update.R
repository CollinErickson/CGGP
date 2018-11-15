updateSG <- function(SG, y, ...) {
  if (is.matrix(y) && ncol(y)>1) {
    #multiY
    logthetaMLEMV(SG=SG, yMV=y, ...)
  } else {
    logthetaMLE(SG=SG, y=y, ...)
  }
}