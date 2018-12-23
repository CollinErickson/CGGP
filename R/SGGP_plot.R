#' Heatmap of SG design depth
#' 
#' The values on the diagonal are largest design depth for that dimension.
#' The off-diagonal values are the largest design depth that both dimensions
#' have been measured at simultaneously.
#' A greater depth means that more points have been measured along that
#' dimension or two-dimensional subspace.
#'
#' @param SGGP SGGP object
#'
#' @return A heat map made from ggplot2
#' @export
#' @references https://stackoverflow.com/questions/14290364/heatmap-with-values-ggplot2
#'
#' @examples
#' # All dimensions should look similar
#' d <- 8
#' SG = SGGPcreate(d,201)
#' SGGPheat(SG)
#' 
#' # The first and fourth dimensions are most active and will have greater depth
#' SG <- SGGPcreate(d=5, batchsize=50)
#' f <- function(x) {cos(2*pi*x[1]*3) + exp(4*x[4])}
#' for (i in 1:1) {
#'   SG <- SGGPfit(SG, Y=apply(SG$design, 1, f))
#'   SG <- SGGPappend(SGGP=SG, batchsize=200)
#' }
#' SG <- SGGPfit(SG, Y=apply(SG$design, 1, f))
#' SGGPheat(SG)
SGGPheat <- function(SGGP) {
  # heatmatrix <- matrix(NaN, SG$d, SG$d)
  skinny <- NULL
  for (i in 1:SGGP$d) {
    # heatmatrix[i,i] <- max(SG$designindex[,i])
    skinny <- rbind(skinny, c(i, i, max(SGGP$uo[,i])))
  }
  for (i in 1:(SGGP$d-1)) {
    for (j in (i+1):SGGP$d) {
      # heatmatrix[i,j] <- heatmatrix[j,i] <- max(apply(SG$designindex[,c(i,j)], 1, min))
      skinny <- rbind(skinny,
                      c(i, j, max(apply(SGGP$uo[,c(i,j)], 1, min))),
                      c(j, i, max(apply(SGGP$uo[,c(i,j)], 1, min)))
      )
    }
  }
  
  skdf <- data.frame(skinny)
  names(skdf) <- c('Var1', 'Var2', 'value')
  ggplot2::ggplot(skdf, ggplot2::aes_string('Var1', 'Var2')) +
    ggplot2::geom_tile(ggplot2::aes_string(fill = 'value')) + 
    ggplot2::geom_text(ggplot2::aes_string(label = 'round(value, 1)')) +
    ggplot2::scale_fill_gradient(low = "white", high = "red") 
}


#' Histogram of measurements at each design depth of each input dimension
#' 
#' A greater design depth signifies a more important dimension.
#' Thus a larger right tail on the histogram are more important variables.
#'
#' @param SG SGGP object
#' @param ylog Should the y axis be put on a log scale?
#'
#' @return Histogram plot made using ggplot2
#' @export
#'
#' @examples
#' # All dimensions should look similar
#' d <- 8
#' SG = SGGPcreate(d,201)
#' SGGPhist(SG)
#' SGGPhist(SG, ylog=F)
#' 
#' # The first dimension is more active and will have greater depth
#' SG <- SGGPreate(d=5, batchsize=10)
#' SG <- SGGPappend(logtheta=c(-2,2,2,2,2), SG=SG, batchsize=100)
#' #SGhist(SG)
SGGPhist <- function(SGGP, ylog=T) {
  p <- ggplot2::ggplot(reshape2::melt(data.frame(SGGP$designindex), id.vars=NULL),
                       ggplot2::aes_string(x='value')) + 
          ggplot2::geom_histogram(binwidth = 1) + ggplot2::facet_grid(variable ~ .)
  if (ylog) {
    p <- p + ggplot2::scale_y_log10() #limits=c(.9999, NA))
  }
  p
}

#' Plot validation prediction errors
#'
#' @param SG SGGP object
#' @param y Measurements at SG$design
#' @param Xval X validation data
#' @param Yval Y validation data
#' @param ypred (optional) Predictions of SG at Xval
#'
#' @return None, makes a plot
#' @export
#'
#' @examples
#' SG <- SGcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2+rnorm(1,0,.01)}
#' y <- apply(SG$design, 1, f1)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' SG <- logthetaVALID(SG=SG, y=y, xval=Xval, yval=Yval)
#' SGGPvalplot(SG=SG, y=y, xval=Xval, yval=Yval)
SGGPvalplot <- function(SG, y, xval, yval, ypred) {
  if (missing(xpred)) {xpred <- SGGPpred(xp=Xval, SG=SG, y=y)}
  errmax <- max(sqrt(ypred$var), abs(Ypred$mean - Yp))
  plot(ypred$mean-Yval, sqrt(ypred$var), xlim=errmax*c(-1,1), ylim=c(0,errmax))#;abline(a=0,b=1,col=2)
  polygon(1.1*errmax*c(0,-2,2),1.1*errmax*c(0,1,1), col=3, density=10, angle=135)
  polygon(1.1*errmax*c(0,-1,1),1.1*errmax*c(0,1,1), col=2, density=30)
  points(SG$mean-Yp, sqrt(SG$var), xlim=errmax*c(-1,1), ylim=c(0,errmax))
}