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
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::scale_x_continuous(breaks = 1:SGGP$d)  +
    ggplot2::scale_y_continuous(breaks = 1:SGGP$d) # labels=c() to set names
}


#' Histogram of measurements at each design depth of each input dimension
#' 
#' A greater design depth signifies a more important dimension.
#' Thus a larger right tail on the histogram are more important variables.
#'
#' @param SGGP SGGP object
#' @param ylog Should the y axis be put on a log scale?
#'
#' @return Histogram plot made using ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' # All dimensions should look similar
#' d <- 8
#' SG = SGGPcreate(d,201)
#' SGGPhist(SG)
#' SGGPhist(SG, ylog=FALSE)
#' 
#' # The first dimension is more active and will have greater depth
#' SG <- SGGPreate(d=5, batchsize=10)
#' SG <- SGGPappend(SGGP=SG, batchsize=100)
#' SGGPhist(SG)
#' }
SGGPhist <- function(SGGP, ylog=TRUE) {
  p <- ggplot2::ggplot(reshape2::melt(data.frame(SGGP$uo), id.vars=NULL),
                       ggplot2::aes_string(x='value'))
  # Tried a power transformation, but I can't get breaks to work as expected
  # p <- p + ggplot2::coord_trans(y=scales::trans_new(name="test", 
  #                                                   transform=function(x) x^.1,
  #                                                   inverse=function(x) x^10,
  #                                                   breaks=function(...) c(0,.3,.5,1),
  #                                                   minor_breaks=c(0,1,2)
  #                                                   )
  #                               )
  # p <- p + ggplot2::coord_trans(y=scales::boxcox_trans(p=.2)  )
  p <- p +ggplot2::geom_histogram(binwidth = 1) + ggplot2::facet_grid(variable ~ .)
  if (ylog) {
    p <- p + ggplot2::scale_y_log10() #limits=c(.9999, NA))
  }
  p <- p + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) # rotates labels from vert to hor
  
  p
}

#' Plot validation prediction errors
#'
#' @param SGGP SGGP object that has been fitted
#' @param Xval X validation data
#' @param Yval Y validation data
#' @param plot_with Should the plot be made with "base" or "ggplot2"?
#' @param d If output is multivariate, which column to use
#'
#' @return None, makes a plot
#' @export
#' @importFrom graphics plot points polygon
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2}
#' y <- apply(SG$design, 1, f1)
#' SG <- SGGPfit(SG, y)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' SGGPvalplot(SGGP=SG, Xval=Xval, Yval=Yval)
SGGPvalplot <- function(SGGP, Xval, Yval, plot_with="ggplot2", d=NULL) {
  ypred <- SGGPpred(xp=Xval, SGGP=SGGP)
  if (!is.null(d)) {
    ypred <- list(mean=ypred$mean[,d], var=ypred$var[,d])
    Yval <- Yval[,d]
  }
  errmax <- max(sqrt(ypred$var), abs(ypred$mean - Yval))
  if (plot_with == "base") {
    plot(ypred$mean-Yval, sqrt(ypred$var), xlim=errmax*c(-1,1), ylim=c(0,errmax))#;abline(a=0,b=1,col=2)
    polygon(1.1*errmax*c(0,-2,2),1.1*errmax*c(0,1,1), col=3, density=10, angle=135)
    polygon(1.1*errmax*c(0,-1,1),1.1*errmax*c(0,1,1), col=2, density=30)
    points(ypred$mean-Yval, sqrt(ypred$var), xlim=errmax*c(-1,1), ylim=c(0,errmax))
  } else {
    
    tdf <- data.frame(err=ypred$mean-Yval, psd=sqrt(ypred$var))
    # ggplot(tdf, aes(x=err, y=psd)) + geom_point()
    values <- data.frame(id=factor(c(1, 2)), value=factor(c('095%','68%')))
    positions <- data.frame(id=rep(values$id, each=3),
                            x=1.1*c(0,errmax*2,-errmax*2, 0,errmax,-errmax),
                            y=1.1*c(0,errmax,errmax,0,errmax,errmax))
    # Currently we need to manually merge the two together
    datapoly <- merge(values, positions, by = c("id"))
    
    # ggplot(datapoly, aes(x = x, y = y)) +
    # geom_polygon(aes(fill = value, group = id))
    # ggplot(tdf, aes(x=err, y=psd)) + geom_polygon(aes(fill = value, group = id, x=x, y=y), datapoly, alpha=.2) + geom_point() +
    # xlab("Predicted - Actual") + ylab("Predicted error") + coord_cartesian(xlim=c(-errmax,errmax), ylim=c(0,errmax))
    ggplot2::ggplot(tdf, ggplot2::aes_string(x='err', y='psd')) + 
      ggplot2::geom_polygon(ggplot2::aes_string(fill = 'value', group = 'id', x='x', y='y'), datapoly, alpha=.2) + 
      ggplot2::geom_point() +
      ggplot2::xlab("Predicted - Actual") + ggplot2::ylab("Predicted error") + 
      ggplot2::coord_cartesian(xlim=c(-errmax,errmax), ylim=c(0,errmax))
    
  }
}

#' Calculate stats for SGGP prediction on validation data
#'
#' @param SGGP SGGP object
#' @param Xval X validation matrix
#' @param Yval Y validation data
#' @param d Dimension if multivariate output
#'
#' @return data frame
#' @export
#' @importFrom stats pnorm dnorm
#' @references Gneiting, Tilmann, and Adrian E. Raftery.
#' "Strictly proper scoring rules, prediction, and estimation."
#' Journal of the American Statistical Association 102.477 (2007): 359-378.
#'
#' @examples
#' SG <- SGGPcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2}
#' y <- apply(SG$design, 1, f1)
#' SG <- SGGPfit(SG, y)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' SGGPvalstats(SGGP=SG, Xval=Xval, Yval=Yval)
SGGPvalstats <- function(SGGP, Xval, Yval, d=NULL) {
  
  ypred <- SGGPpred(xp=Xval, SGGP=SGGP)
  if (!is.null(d)) {
    ypred <- list(mean=ypred$mean[,d], var=ypred$var[,d])
    Yval <- Yval[,d]
  }
  
  m <- ypred$mean
  v <- pmax(ypred$var, 0)
  s <- sqrt(v)
  z <- (Yval - m) / s
  RMSE <- sqrt(mean((ypred$mean - Yval)^2))
  score <- mean((Yval-ypred$mean)^2/ypred$var+log(ypred$var))
  CRPscore <- - mean(s * (1/sqrt(pi) - 2*dnorm(z) - z * (2*pnorm(z) - 1)))
  coverage <- mean((Yval<= ypred$mean+1.96*sqrt(ypred$var)) & 
                     (Yval>= ypred$mean-1.96*sqrt(ypred$var)))
  
  data.frame(RMSE=RMSE, score=score, CRPscore=CRPscore, coverage=coverage)
  
  
}

#' Plot correlation samples
#' 
#' Plot samples for a given correlation function and parameters.
#' Useful for getting an idea of what the correlation parameters mean
#' in terms of smoothness.
#'
#' @param Corr Correlation function or SGGP object.
#' If SGGP object, it will make plots for thetaMAP,
#' the max a posteriori theta.
#' @param theta Parameters for Corr
#' @param numlines Number of sample paths to draw
#' @param plot_with Should "base" or "ggplot2" be used to make the plot?
#' @param zero Should the sample paths start at y=0?
#'
#' @return Plot
#' @export
#' @importFrom graphics par
#'
#' @examples
#' SGGP_internal_CorrPlot()
#' SGGP_internal_CorrPlot(theta=c(-2,-1,0,1))
#' 
#' SG <- SGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+sin(2*pi*x[2]*3)}
#' y <- apply(SG$design, 1, f)
#' SG <- SGGPfit(SG, Y=y)
#' SGGP_internal_CorrPlot(SG)
SGGP_internal_CorrPlot <- function(Corr=SGGP_internal_CorrMatGaussian, theta=NULL,
                                   numlines=20, plot_with="ggplot",
                                   zero=TRUE) {
  # Points along x axis
  n <- 100
  xl <- seq(0,1,l=n)
  
  if (inherits(Corr, "SGGP")) {
    if (!is.null(theta)) {theta <- Corr$thetaMAP}
    Corr <- Corr$CorrMat
  }
  
  nparam <- Corr(return_numpara=TRUE)
  if (is.null(theta)) {theta <- rep(0, nparam)}
  
  ncorr <- length(theta) / nparam
  for (i in 1:ncorr) {
    # Can change px.mean and px.cov to be conditional on other data
    px.mean <- rep(0,n)
    px.cov <- Corr(xl, xl, theta=theta[1:nparam + nparam*(i-1)])
    # Generate sample paths
    samplepaths <- try(newy <- MASS::mvrnorm(n=numlines, mu=px.mean, Sigma=px.cov))
    if (inherits(samplepaths, "try-error")) {
      message("Adding nugget to cool1Dplot")
      nug <- 1e-6
      Sigma.try2 <- try(
        newy <- MASS::mvrnorm(n=numlines, mu=px.mean,
                              Sigma=px.cov + diag(nug, nrow(px.cov)))
      )
      if (inherits(Sigma.try2, "try-error")) {
        stop("Can't do cool1Dplot")
      }
    }
    
    # Set so all start at 0.
    if (zero) {
      samplepathsplot <- sweep(samplepaths, 1, samplepaths[,1])
    }
    
    if (plot_with == "base") {
      # Put plots in column if more than one
      if (i==1 && ncorr>1) {
        orig.mfrow <- par()$mfrow
        par(mfrow=c(ncorr, 1))
      }
      plot(xl, samplepathsplot[1,], type='l',
           ylim=c(min(samplepathsplot), max(samplepathsplot)),
           ylab="Sample path", xlab="x")
      for (i in 2:numlines) {
        points(xl, samplepathsplot[i,], type='l', col=i)
      }
    } else { # Use ggplot2
      # ggplot2::ggplot(cbind(reshape2::melt(data.frame(t(samplepathsplot)), id.vars=c()),
      #                       x=rep(xl, numlines)), 
      #                 ggplot2::aes_string(x="x", y="value", color="variable")) + 
      #   ggplot2::geom_line() + ggplot2::theme(legend.position="none")
      newdf <- cbind(reshape2::melt(data.frame(t(samplepathsplot)), id.vars=c()),
                    x=rep(xl, numlines), d=i)
      if (i==1) {ggdf <- newdf}
      else {ggdf <- rbind(ggdf, newdf)}
    }
  }
  if (plot_with == "base") {
    # Reset graphical parameters
    par(mfrow=orig.mfrow)
  } else { # Use ggplot2
    # Return plot
    p <- ggplot2::ggplot(ggdf, 
                    ggplot2::aes_string(x="x", y="value", color="variable")) + 
      ggplot2::geom_line() + ggplot2::theme(legend.position="none")
    if (i > 1) {
      p <- p + ggplot2::facet_grid(d ~ .)
    }
    p
  }
}
