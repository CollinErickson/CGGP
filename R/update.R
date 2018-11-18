#' Update SGGP object with new data
#'
#' @param SG SGGP object
#' @param y Measured output values
#' @param restarts Number of optimization restarts
#' @param ... Not used
#' @param ynew If only returning new output values from last append, use this
#' @param method Whether parameters should be updated using MLE or validation.
#'
#' @return Updated SG object
#' @export
#'
#' @examples
#' f <- function(x) {(x[1]-.5)^2 + exp(-x[2]) + sin(2*pi*x[3])}
#' SG <- SGcreate(d=3, batchsize=100)
#' SG <- updateSG(SG, apply(SG$design, 1, f))
#' #SG <- SGappend(SG=SG, batchsize=20)
updateSG <- function(SG, y, restarts=4, ..., ynew, method="MLE") {
  # Can be multioutput or uni
  # Can be updating params or not
  # Can be returning new Y or not
  
  # If only giving in ynew, create y 
  if (!missing(ynew)) {
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
    # lt0[(rowsset+1):nrow(lt0),] <- runif((nrow(lt0)-(rowsset))*SG$d, -2, 2)
    lt0[(rowsset+1):nrow(lt0),] <- lhs::maximinLHS(n=nrow(lt0)-(rowsset), k=SG$d)*4-2
  }
  print("Going to run lt0"); print(lt0)

  # Run optimization. Use lapply so it can easily be parallelized later
  optout <- lapply(1:restarts,
         function(i) {f(SG, y, logtheta0=lt0[i,], return_optim=T)})
  print("optout is"); print(optout)
  

  # Find best logtheta from all restarts
  bestlogthetaind <- which.min(sapply(optout, function(x) {x$value}))
  bestlogtheta <- optout[[bestlogthetaind]]$par
  cat('best ind is', bestlogthetaind, 'bestlogtheta is', bestlogtheta, '\n')
  
  # Set new y's
  SG$y <- y
  
  # Set best logtheta
  SG$logtheta <- bestlogtheta

  # Save pw with SG
  SG$pw <- calculate_pw(SG=SG, y=y-mean(y), logtheta=SG$logtheta)
  
  SG
}
