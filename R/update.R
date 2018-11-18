updateSG <- function(SG, y, ...) {
  # Can be multioutput or uni
  # Can be updating params or not
  # Can be returning new Y or not
  if (is.matrix(y) && ncol(y)>1) {
    #multiY
    logthetaMLEMV(SG=SG, yMV=y, ...)
  } else {
    logthetaMLE(SG=SG, y=y, ...)
  }
}

updateSG2 <- function(SG, y, restarts=4, ..., ynew, method="MLE") {
  # Can be multioutput or uni
  # Can be updating params or not
  # Can be returning new Y or not
  
  # If only giving in ynew, create y 
  if (!is.missing(ynew)) {
    if (!missing(y)) {if (any(y != SG$y)) {stop("ynew was given but y != SG$y")}}
    if (is.matrix(SG$y)) {
      y <- rbind(SG$y, ynew)
    } else { # Just a vector
      y <- c(SG$y, ynew)
    }
  }
  
  if (is.matrix(y) && ncol(y)>1) {
    #multiY
    if (tolower(method)=="mle") {f <- logthetaMLEMV}
    else if (tolower(method)=="validation") {stop('logthetaVALIDMV not implemented')}
    else {stop("Method not recognized")}
  } else {
    if (tolower(method)=="mle") {f <- logthetaMLE}
    else if (tolower(method)=="validation") {stop('logthetaVALID not implemented')}
    else {stop("Method not recognized")}
  }
  
  # --------------------------------------
  # Now run optimization with restarts
  # --------------------------------------
  
  # Create initial parameters for optimization, lt0=logtheta0
  lt0 <- matrix(0, restarts, SG$d)
  rowsset <- 0
  # Set first row to current logtheta if possible
  if (!is.null(SG$logtheta)) {
    lt0[1,] <- SG$logtheta
    # rowsleft <- rowsleft - 1
    rowsset <- rowsset + 1
  }
  # Set next row to all 0's
  lt0[rowsset+1,] <- rep(0, SG$d)
  rowsset <- rowsset + 1
  
  # If still rows left, set using random uniform
  if (rowsset < restarts) {
    lt0[(rowsset+1):nrow(lt0),] <- runif((nrow(lt0)-(rowsset))*SG$d, -2, 2)
  }
  print("Going to run lt0"); print(lt0)

  # Run optimization. Use lapply so it can easily be parallelized later
  optout <- lapply(1:restarts,
         function(i) {f(SG, y, logtheta0=lt0[i,])})
  
  # Find best from optout
  # Set SG$logtheta, update SG$pw
  
  SG$y <- y
}
