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
#' SG = SGGPcreate(3,201)
#' print(SG)
#' f <- function(x) {x[1]+exp(x[2]) + log(x[3]+4)}
#' y <- apply(SG$design, 1, f)
#' SG <- SGGPfit(SG, y)
#' print(SG)
print.SGGP <- function(x, ...) {
  s <- paste0(
    c(
      "SGGP object\n",
      "   d = ", x$d, '\n',
      "   number of design points = ", nrow(x$design), '\n',
      "   number of unevaluated design points = ", if (is.null(x$design_unevaluated)) 0 else nrow(x$design_unevaluated), '\n',
      #"  theta = ", if (all(x$thetaMAP==0)) "(not yet calculated)" else x$thetaMAP, '\n',
      "   Available functions:\n",
      "     - SGGPfit(SGGP, Y) to update parameters with new data\n",
      "     - sGGPpred(SGGP, xp) to predict at new points\n",
      "     - SGGPappend(SGGP, batchsize) to add new design points\n"
    )
  )
  cat(s, sep="")
}

#' S3 predict method for SGGP
#' 
#' Passes to SGGPpred
#' 
#' @param object SGGP object
#' @param ... Required for S3 consistency
#'
#' @rdname SGGPpred
#' @export
predict.SGGP <- function(object, xp, ...) {
  SGGPpred(SGGP=object, xp=xp)
}

#' S3 plot method for SGGP
#' 
#' There are a few different plot functions for SGGP objects:
#' `SGGPblockplot`, `SGGPcorrplot`, `SGGPprojectionplot`,
#' `SGGPvalplot`, `SGGPheat`, and `SGGPhist`.
#' Currently `SGGPblockplot` is the default plot object.
#'
#' @param x SGGP object
#' @param y Don't use
#' @param ... Passed to SGGPblockplot
#'
#' @return Either makes plot or returns plot object
#' @export
#'
#' @examples
#' SG = SGGPcreate(3,100)
#' plot(SG)
plot.SGGP <- function(x, y, ...) {
  SGGPblockplot(x, ...)
}