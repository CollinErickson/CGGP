#' Print SGGP object
#' 
#' Default print as a list is bad since there's a lot of elements.
#'
#' @param x SGGP object
#' @param ... Passed to print
#'
#' @return String to be printed
#' @export
#'
#' @examples
#' SG = SGcreate(3,201)
#' print(SG)
print.SGGP <- function(x, ...) {
  s <- paste0(
    "SGGP object\n",
    "  d = ", x$d, '\n',
    "  number of design points = ", nrow(x$design), '\n',
    "  logtheta = ", if (is.null(x$logtheta)) "(not yet calculated)" else x$logtheta, '\n',
    "To update parameters, run updateSG()\n",
    "To predict, run SGGPpred()\n",
    "To get new design points, run SGappend()\n"
  )
  cat(s)
}

#' S3 predict method for SGGP
#' 
#' Passes to SGGPpred
#' 
#' @param object SGGP object
#'
#' @rdname SGGPpred
#' @export
predict.SGGP <- function(object, xp, y, logtheta, theta, ...) {
  SGGPpred(SG=object, xp=xp, y=y, ..., logtheta=logtheta, theta=theta)
}