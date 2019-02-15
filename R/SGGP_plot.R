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
#' \dontrun{
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
#' # SG <- SGGPfit(SG, Y=apply(SG$design, 1, f))
#' SGGPheat(SG)
#' }
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
#' SG <- SGGPcreate(d=5, batchsize=10)
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


#' SGGP block plot
#' 
#' Plot the 2D projections of the blocks of an SGGP object.
#'
#' @param SGGP SGGP object
#' @param singleplot If only two dimensions, should a single plot be made?
#'
#' @return ggplot2 plot
#' @export
#'
#' @examples
#' \dontrun{
#' # The first and fourth dimensions are most active and will have greater depth
#' ss <- SGGPcreate(d=5, batchsize=50)
#' f <- function(x) {cos(2*pi*x[1]*3) + x[3]*exp(4*x[4])}
#' ss <- SGGPfit(ss, Y=apply(ss$design, 1, f))
#' ss <- SGGPappend(SGGP=ss, batchsize=100)
#' SGGPblockplot(ss)
#' 
#' mat <- matrix(c(1,1,1,2,2,1,2,2,1,3), ncol=2, byrow=TRUE)
#' SGGPblockplot(mat)
#' }
SGGPblockplot <- function(SGGP, singleplot=TRUE) {
  
  if (inherits(SGGP, "SGGP")) {
    d <- SGGP$d
    uo <- SGGP$uo[1:SGGP$uoCOUNT,]
  } else if (is.matrix(SGGP)) {
    d <- ncol(SGGP)
    uo <- SGGP
  } else {stop("blockplot only works on SGGP or matrix")}
  if (d>2) {singleplot <- FALSE} # Can only use singleplot in 2D
  
  alldf <- NULL
  # ds <- c(1,2)
  d1set <- if (singleplot) {1} else {1:d}
  for (d1 in d1set) {
    d2set <- setdiff(1:d, d1)
    for (d2 in d2set) {
      ds <- c(d1, d2)
      uods <- uo[,ds]
      uodsdf <- data.frame(uods)
      uods.unique <- plyr::ddply(uodsdf, c("X1", "X2"), function(x) data.frame(count=nrow(x)))
      
      gdf <- uods.unique
      gdf$xmin <- gdf$X1-1
      gdf$xmax <- gdf$X1
      gdf$ymin <- gdf$X2-1
      gdf$ymax <- gdf$X2
      gdf$d1 <- d1
      gdf$d2 <- d2
      alldf <- if (is.null(alldf)) gdf else rbind(alldf, gdf)
    }
  }
  p <- ggplot2::ggplot(alldf) +
    ggplot2::geom_rect(ggplot2::aes_string(xmin="xmin", xmax="xmax",
                                           ymin="ymin", ymax="ymax",
                                           fill="count"),
                       color="black")+
    ggplot2::scale_fill_gradient(low = "yellow", high = "red")
  
  if (!singleplot) {p <- p + ggplot2::facet_grid(d1 ~ d2)}
  else {p <- p+ggplot2::xlab("X1") + ggplot2::ylab("X2")}
  p
}



#' Plot validation prediction errors
#'
#' @param predmean Predicted mean
#' @param predvar Predicted variance
#' @param Yval Y validation data
#' @param plot_with Should the plot be made with "base" or "ggplot2"?
#' @param d If output is multivariate, which column to use. Will do all if
#' left as NULL.
#'
#' @return None, makes a plot
#' @export
#' @importFrom graphics plot points polygon
#'
#' @examples
#' x <- matrix(runif(100*3), ncol=3)
#' f1 <- function(x){x[1]+x[2]^2}
#' y <- apply(x, 1, f1)
#' # Create a linear model on the data
#' mod <- lm(y ~ ., data.frame(x))
#' # Predict at validation data
#' Xval <- matrix(runif(3*100), ncol=3)
#' mod.pred <- predict.lm(mod, data.frame(Xval), se.fit=TRUE)
#' # Compare to true results
#' Yval <- apply(Xval, 1, f1)
#' valplot(mod.pred$fit, mod.pred$se.fit^2, Yval=Yval)
valplot <- function(predmean, predvar, Yval, plot_with="ggplot2", d=NULL) {
  
  if (!is.null(d)) {
    predmean <- predmean[,d]
    predvar  <- predvar[,d]
    Yval <- Yval[,d]
  }
  errmax <- max(sqrt(predvar), abs(predmean - Yval))
  if (plot_with == "base") {
    plot(predmean-Yval, sqrt(predvar), xlim=errmax*c(-1,1), ylim=c(0,errmax))#;abline(a=0,b=1,col=2)
    polygon(1.1*errmax*c(0,-2,2),1.1*errmax*c(0,1,1), col=3, density=10, angle=135)
    polygon(1.1*errmax*c(0,-1,1),1.1*errmax*c(0,1,1), col=2, density=30)
    points(predmean-Yval, sqrt(predvar), xlim=errmax*c(-1,1), ylim=c(0,errmax))
  } else {
    errs <- unname(predmean-Yval)
    psds <- unname(sqrt(predvar))
    if (!is.matrix(errs) || ncol(errs)==1) { # vector, single output
      tdf <- data.frame(err=errs, psd=psds)
    } else { # multiple outputs, need to melt
      tdf <- cbind(reshape2::melt(errs), reshape2::melt(psds))[,c(3,6,2)]
      tdf <- data.frame(tdf)
      names(tdf) <- c("err", "psd", "outputdim")
    }
    # ggplot(tdf, aes(x=err, y=psd)) + geom_point()
    values <- data.frame(id=factor(c(2,1)), Error=factor(c('68%','95%')))
    positions <- data.frame(id=rep(values$id, each=3),
                            x=1.1*c(0,errmax*2,-errmax*2, 0,errmax,-errmax),
                            y=1.1*c(0,errmax,errmax,0,errmax,errmax))
    # Currently we need to manually merge the two together
    datapoly <- merge(values, positions, by = c("id"))
    
    # ggplot(datapoly, aes(x = x, y = y)) +
    # geom_polygon(aes(fill = value, group = id))
    # ggplot(tdf, aes(x=err, y=psd)) + geom_polygon(aes(fill = value, group = id, x=x, y=y), 
    #              datapoly, alpha=.2) + geom_point() +
    # xlab("Predicted - Actual") + ylab("Predicted error") + 
    #   coord_cartesian(xlim=c(-errmax,errmax), ylim=c(0,errmax))
    p <- ggplot2::ggplot(tdf, ggplot2::aes_string(x='err', y='psd')) + 
      ggplot2::geom_polygon(ggplot2::aes_string(fill = 'Error', group = 'id', x='x', y='y'),
                            datapoly[datapoly$id==2,], alpha=.2) + # Separate these so it does the right order
      ggplot2::geom_polygon(ggplot2::aes_string(fill = 'Error', group = 'id', x='x', y='y'),
                            datapoly[datapoly$id==1,], alpha=.2) + 
      ggplot2::geom_point() +
      ggplot2::xlab("Predicted - Actual") + ggplot2::ylab("Predicted error") + 
      ggplot2::coord_cartesian(xlim=c(-errmax,errmax), ylim=c(0,max(psds))) #errmax))
    if (is.matrix(errs) && ncol(errs) > 1) {
      p <- p + ggplot2::facet_wrap("outputdim")
    }
    p
  }
}



#' Plot validation prediction errors for SGGP object
#'
#' @param SGGP SGGP object that has been fitted
#' @param Xval X validation data
#' @param Yval Y validation data
#' @param plot_with Should the plot be made with "base" or "ggplot2"?
#' @param d If output is multivariate, which column to use. Will do all if
#' left as NULL.
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
  ypred <- SGGPpred(SGGP=SGGP, xp=Xval)
  valplot(ypred$mean, ypred$var, Yval, plot_with=plot_with, d=d)
}


#' Calculate stats for prediction on validation data
#'
#' @param predmean Predicted mean
#' @param predvar Predicted variance
#' @param Yval Y validation data
#' @param bydim If multiple outputs, should it be done separately by dimension?
#'
#' @return data frame
#' @export
#' @importFrom stats pnorm dnorm cor
#' @references Gneiting, Tilmann, and Adrian E. Raftery.
#' "Strictly proper scoring rules, prediction, and estimation."
#' Journal of the American Statistical Association 102.477 (2007): 359-378.
#'
#' @examples
#' valstats(c(0,1,2), c(.01,.01,.01), c(0,1.1,1.9))
#' valstats(cbind(c(0,1,2), c(1,2,3)),
#'          cbind(c(.01,.01,.01),c(.1,.1,.1)),
#'          cbind(c(0,1.1,1.9),c(1,2,3)))
#' valstats(cbind(c(0,1,2), c(8,12,34)),
#'          cbind(c(.01,.01,.01),c(1.1,.81,1.1)),
#'          cbind(c(0,1.1,1.9),c(10,20,30)), bydim=FALSE)
#' valstats(cbind(c(.8,1.2,3.4), c(8,12,34)),
#'          cbind(c(.01,.01,.01),c(1.1,.81,1.1)),
#'          cbind(c(1,2,3),c(10,20,30)), bydim=FALSE)
valstats <- function(predmean, predvar, Yval, bydim=TRUE, metrics) {
  if ((is.matrix(predmean) && 
       (any(dim(predmean)!=dim(predvar)) || any(dim(predmean)!=dim(Yval)))
  ) ||
  (length(predmean)!=length(predvar) || length(predmean)!=length(Yval))
  ) {
    stop("Shape of predmean, predvar, and Yval must match in valstats")
  }
  if (bydim && is.matrix(predmean) && ncol(predmean) > 1) {
    return(do.call("rbind",
                   lapply(1:ncol(predmean),
                          function(i) {
                            valstats(predmean=predmean[,i], predvar=predvar[,i], Yval=Yval[,i])
                          })))
  }
  
  m <- predmean
  v <- pmax(predvar, 0)
  s <- sqrt(v)
  z <- (Yval - m) / s
  RMSE <- sqrt(mean((predmean - Yval)^2))
  score <- mean((Yval-predmean)^2/predvar+log(predvar))
  CRPscore <- - mean(s * (1/sqrt(pi) - 2*dnorm(z) - z * (2*pnorm(z) - 1)))
  coverage <- mean((Yval<= predmean+1.96*sqrt(predvar)) & 
                     (Yval>= predmean-1.96*sqrt(predvar)))
  corr <- cor(c(predmean), c(Yval))
  R2 <- 1 - (sum((Yval - predmean)^2) / sum((Yval - mean(Yval))^2))
  # Return df with values
  out <- data.frame(RMSE=RMSE, score=score, CRPscore=CRPscore, coverage=coverage,
                    corr=corr, R2=R2)
  
  # Add normalized RMSE, dividing MSE in each dimension by variance
  if (is.matrix(predmean) && ncol(predmean) > 1) {
    # out$RMSEnorm <- sqrt(mean(sweep((predmean - Yval)^2, 2, apply(Yval, 2, var), "/")))
    out$RMSEnorm <- sqrt(mean(colMeans((predmean - Yval)^2) / apply(Yval, 2, var)))
  }
  
  # Check if user wants more metrics
  if (!missing(metrics)) {
    if (is.function(metrics)) {
      out$metric <- metrics(predmean, predvar, Yval)
    } else if (is.list(metrics)) {
      metrics.names <- names(metrics)
      if (is.null(metrics.names)) {
        metrics.names <- paste0("metrics", 1:length(metrics))
      }
      for (iii in 1:length(metrics)) {
        out[[metrics.names[iii]]] <- metrics[[iii]](predmean, predvar, Yval)
      }
    }
  }
  
  out
}

#' Calculate stats for SGGP prediction on validation data
#'
#' @param SGGP SGGP object
#' @param Xval X validation matrix
#' @param Yval Y validation data
#' @param bydim If multiple outputs, should it be done separately by dimension?
#' @param fullBayesian Should prediction be done fully Bayesian? Much slower.
#' Averages over theta samples instead of using thetaMAP.
#'
#' @return data frame
#' @export
#'
#' @examples
#' \dontrun{
#' SG <- SGGPcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2}
#' y <- apply(SG$design, 1, f1)
#' SG <- SGGPfit(SG, y)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- apply(Xval, 1, f1)
#' SGGPvalstats(SGGP=SG, Xval=Xval, Yval=Yval)
#' SGGPvalstats(SGGP=SG, Xval=Xval, Yval=Yval, fullBayesian=TRUE)
#' 
#' # Multiple outputs
#' SG <- SGGPcreate(d=3, batchsize=100)
#' f1 <- function(x){x[1]+x[2]^2}
#' f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
#' y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
#' y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
#' y <- cbind(y1, y2)
#' SG <- SGGPfit(SG, Y=y)
#' Xval <- matrix(runif(3*100), ncol=3)
#' Yval <- cbind(apply(Xval, 1, f1),
#'               apply(Xval, 1, f2))
#' SGGPvalstats(SG, Xval, Yval)
#' SGGPvalstats(SG, Xval, Yval, bydim=FALSE)
#' }
SGGPvalstats <- function(SGGP, Xval, Yval, bydim=TRUE, fullBayesian=FALSE) {
  # Make predictions
  ypred <- SGGPpred(SGGP=SGGP, xp=Xval, fullBayesian=fullBayesian)
  # Use valstats to get df with values
  valstats(ypred$mean, ypred$var, Yval=Yval, bydim=bydim)
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
#' @param outdims Which output dimensions should be used?
#'
#' @return Plot
#' @export
#' @importFrom graphics par
#'
#' @examples
#' \dontrun{
#' SGGPcorrplot()
#' SGGPcorrplot(theta=c(-2,-1,0,1))
#' 
#' SG <- SGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+sin(2*pi*x[2]*3)}
#' y <- apply(SG$design, 1, f)
#' SG <- SGGPfit(SG, Y=y)
#' SGGPcorrplot(SG)
#' }
SGGPcorrplot <- function(Corr=SGGP_internal_CorrMatGaussian, theta=NULL,
                         numlines=20, plot_with="ggplot",
                         outdims=NULL,
                         zero=TRUE) {
  # Points along x axis
  n <- 100
  xl <- seq(0,1,l=n)
  
  if (inherits(Corr, "SGGP")) {
    if (is.null(theta)) {theta <- Corr$thetaMAP}
    Corr <- Corr$CorrMat
  }
  
  nparam <- Corr(return_numpara=TRUE)
  if (is.null(theta)) {theta <- rep(0, nparam)}
  
  thetadims <- if (is.matrix(theta)) {ncol(theta)} else {1}
  if (is.null(outdims)) {
    outdims <- 1:thetadims
  } else {
    theta <- theta[,outdims, drop=FALSE]
  }
  numoutdims <- length(outdims)
  
  thetadims <- if (is.matrix(theta)) {ncol(theta)} else {1}
  
  # ncorr <- length(theta) / nparam
  numindims <- length(theta) / nparam / numoutdims
  indims <- 1:numindims
  ggdf <- NULL
  for (indim in indims) {
    for (outdim in outdims) {
      theta.thisloop <- if (thetadims>1) {theta[1:nparam + nparam*(indim-1), outdim]} 
      else {theta[1:nparam + nparam*(indim-1)]}
      # Can change px.mean and px.cov to be conditional on other data
      px.mean <- rep(0,n)
      px.cov <- Corr(xl, xl, theta=theta.thisloop)
      # Generate sample paths
      samplepaths <- newy <- MASS::mvrnorm(n=numlines, mu=px.mean, Sigma=px.cov)
      
      # Set so all start at 0.
      if (zero) {
        samplepathsplot <- sweep(samplepaths, 1, samplepaths[,1])
      }
      
      if (plot_with == "base") {
        # Put plots in column if more than one
        if (indim==1 && outdim==1 && numindims+numoutdims>2) {
          orig.mfrow <- par()$mfrow
          par(mfrow=c(numindims, numoutdims))
        }
        plot(xl, samplepathsplot[1,], type='l',
             ylim=c(min(samplepathsplot), max(samplepathsplot)),
             ylab="Sample path", xlab="x")
        for (numline in 2:numlines) {
          points(xl, samplepathsplot[numline,], type='l', col=numline)
        }
      } else { # Use ggplot2
        # ggplot2::ggplot(cbind(reshape2::melt(data.frame(t(samplepathsplot)), id.vars=c()),
        #                       x=rep(xl, numlines)), 
        #                 ggplot2::aes_string(x="x", y="value", color="variable")) + 
        #   ggplot2::geom_line() + ggplot2::theme(legend.position="none")
        newdf <- cbind(reshape2::melt(data.frame(t(samplepathsplot)), id.vars=c()),
                       x=rep(xl, numlines), d=paste0("X",indim), outdim=paste0("Y",outdim))
        if (is.null(ggdf)) {ggdf <- newdf}
        else {ggdf <- rbind(ggdf, newdf)}
      }
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
    if (numindims > 1 && numoutdims==1) {
      p <- p + ggplot2::facet_grid(d ~ .)
    } else if (numindims > 1 && numoutdims > 1) {
      p <- p + ggplot2::facet_grid(d ~ outdim)
    }
    p
  }
}


#' SGGP projection plot
#' 
#' Show prediction plots when projected down to one dimension.
#' Most useful when setting all values to 0.5 because it will
#' have the most points.
#'
#' @param SGGP  SGGP object
#' @param proj Point to project onto
#' @param np Number of points to use along each dimension
#' @param color Color to make error region
#' @param outdims If multiple outputs, which of them should be plotted?
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' d <- 5
#' f1 <- function(x){x[1]+x[2]^2 + cos(x[3]^2*2*pi*4) - 3.3}
#' s1 <- SGGPcreate(d, 200)
#' s1 <- SGGPfit(s1, apply(s1$design, 1, f1))
#' #s1 <- SGGPappend(s1, 200)
#' #s1 <- SGGPfit(s1, apply(s1$design, 1, f1))
#' SGGPprojectionplot(s1)
#' SGGPprojectionplot(s1, 0.)
#' SGGPprojectionplot(s1, s1$design[nrow(s1$design),])
#' }
SGGPprojectionplot <- function(SGGP, proj=.5, np=300, color="pink", outdims) {
  if (!is.null(SGGP$design_unevaluated)) {stop("SGGP must be updated with all data")}
  if (length(proj) == 1) {proj <- rep(proj, SGGP$d)}
  if (length(proj) != SGGP$d) {stop("proj should be of length SGGP$d or 1")}
  d <- SGGP$d
  
  tdfall <- NULL
  pointdfall <- NULL
  
  Y <- as.matrix(SGGP$Y)
  if (missing(outdims)) {
    numoutdims <- ncol(Y)
    outdims <- 1:numoutdims
  } else {
    numoutdims <- length(outdims)
  }
  
  for (d5 in 1:d) {
    xl <- seq(0,1,l=np)
    m <- matrix(proj,np,d, byrow=T)
    m[,d5] <- xl
    # Won't work if it used PCA or shared parameters
    
    p5 <- try(SGGPpred(SGGP, m, outdims=outdims), silent = TRUE)
    if (inherits(p5, "try-error")) {
      p5 <- SGGPpred(SGGP, m)
    }
    # p5 %>% str
    poly <- cbind(c(rep(xl,each=2))[-c(1,2*np)])
    # If multiple outputs, get all of them
    for (outdim in outdims) {
      if (numoutdims > 1) {
        tdf <- data.frame(mean=p5$mean[,outdim], var=p5$var[,outdim])
      } else { # Single out dim
        tdf <- as.data.frame(p5)
      }
      tdf$outdim <- outdim
      tdf$var <- pmax(0, tdf$var) # No negative values
      tdf$sd <- sqrt(tdf$var)
      tdf$meanp2sd <- tdf$mean + 2*tdf$sd
      tdf$meanm2sd <- tdf$mean - 2*tdf$sd
      tdf$x <- xl
      tdf$d <- d5
      tdfall <- rbind(tdfall, tdf)
    }
    
    w2.5 <- apply(SGGP$design[,-d5, drop=FALSE], 1, function(x) all(abs(x - proj[-d5]) < 1e-8))
    x2.5 <- SGGP$design[w2.5,, drop=FALSE]
    y2.5 <- Y[w2.5,, drop=FALSE]
    # plot(x2.5[,d5], y2.5)
    pointdf <- NULL
    if (length(y2.5) > 0) {
      for (outdim in outdims) {
        pointdf <- rbind(pointdf, data.frame(x=x2.5[,d5], y=y2.5[,outdim], d=d5, outdim=outdim))
      }
    }
    pointdfall <- rbind(pointdfall, pointdf)
  }
  
  tdfall$d <- paste0("X", tdfall$d)
  tdfall$outdim <- paste0("Y", tdfall$outdim)
  if (!is.null(pointdfall)) {
    pointdfall$d <- paste0("X", pointdfall$d)
    pointdfall$outdim <- paste0("Y", pointdfall$outdim)
  }
  p <- ggplot2::ggplot(tdfall, ggplot2::aes_string(x='x')) + 
    ggplot2::geom_ribbon(ggplot2::aes_string(ymin='meanm2sd', ymax='meanp2sd'), color=color, fill=color) +
    ggplot2::geom_line(ggplot2::aes_string(y='mean'))
  if (!is.null(pointdfall)) {
    p <- p + ggplot2::geom_point(ggplot2::aes_string(x='x', y='y'), data=pointdfall)
  }
  if (numoutdims == 1) {
    p <- p + ggplot2::facet_grid(d ~ .)
  } else {
    p <- p +ggplot2::facet_grid(outdim ~ d, scales="free_y")
  }
  p
}


#' Plot something similar to a semivariogram
#' 
#' It's not actually a variogram or semivariogram.
#' It shows how the correlation function falls off as distance increases.
#'
#' @param SGGP SGGP object
#' @param facet How should the plots be faceted? If 1, in a row,
#' if 2, in a column, if 3, wrapped around.
#' @param outdims Which output dimensions should be shown.
#'
#' @return ggplot2 object
#' @export
#'
#' @examples 
#' SG <- SGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+x[3]^.4*sin(2*pi*x[2]^2*3) + .1*exp(3*x[3])}
#' y <- apply(SG$design, 1, f)
#' SG <- SGGPfit(SG, Y=y)
#' SGGPvariogram(SG)
SGGPvariogram <- function(SGGP, facet=1, outdims=NULL) {
  vdf <- NULL
  xl <- seq(0,1,l=101)
  if (is.null(outdims)) {
    outdims <- if (is.matrix(SGGP$Y)) {1:ncol(SGGP$Y)} else {1}
  }
  for (di in 1:SGGP$d) {
    for (outdim in outdims) {
      theta.thisiter <- if (is.matrix(SGGP$thetaMAP)) {SGGP$thetaMAP[1:SGGP$numpara + (di-1)*SGGP$numpara, outdim]
      } else {
        SGGP$thetaMAP[1:SGGP$numpara + (di-1)*SGGP$numpara]
      }
      tmp <- c(SGGP$CorrMat(c(0), xl, theta.thisiter))
      tmpdf <- data.frame(x=xl, y=tmp, d=paste0("X",di), outdim=paste0("Y", outdim))
      vdf <- if (is.null(vdf)) tmpdf else rbind(vdf, tmpdf)
    }
  }
  p <- ggplot2::ggplot(vdf, ggplot2::aes_string(x='x', y='y')) + ggplot2::geom_line() + ggplot2::ylim(c(0,1))
  # p <- p + facet_grid(d ~ .)
  if (facet==1) {
    p <- p + ggplot2::facet_grid(outdim ~ d)
    p <- p + ggplot2::scale_x_continuous(breaks=c(0,1))
  } else if (facet==2) {
    p <- p + ggplot2::facet_grid(d ~ outdim)
  } else if (facet==3) {
    p <- p + ggplot2::facet_wrap(d ~ outdim)
  } else {
    stop("facet is not 1, 2, or 3")
  }
  p
}


#' Plot theta samples
#'
#' @param SGGP SGGP object
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' gs <- SGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+x[3]^.4*sin(2*pi*x[2]^2*3) + .1*exp(3*x[3])}
#' y <- apply(gs$design, 1, f)
#' gs <- SGGPfit(gs, Y=y)
#' SGGPplottheta(gs)
SGGPplottheta <- function(SGGP) {
  # stripchart(data.frame(t(SGGP$thetaPostSamples)), xlim=c(-1,1))
  # stripchart(data.frame(t(SGGP$thetaMAP)), add=T, col=2, pch=17)
  if (is.matrix(SGGP$thetaMAP)) {
    tsamp <- NULL
    tmap <- NULL
    for (i in 1:dim(SGGP$thetaPostSamples)[3]) {
      tsamp <- rbind(tsamp,
                  reshape2::melt(SGGP$thetaPostSamples[,,i]), pdim=i)
      tmap <- rbind(tmap,
                  data.frame(value=SGGP$thetaMAP[,i], Var1=1:nrow(SGGP$thetaMAP), pdim=i))
    }
  } else { # single output
    tsamp <- reshape2::melt(SGGP$thetaPostSamples)
    tmap <- data.frame(value=SGGP$thetaMAP, Var1=1:length(SGGP$thetaMAP))
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data=tmap, mapping=ggplot2::aes_string("value", "Var1"), color="green", size=5) + 
    ggplot2::geom_point(data=tsamp, mapping=ggplot2::aes_string("value", "Var1"))
  if (is.matrix(SGGP$thetaMAP)) {
    p <- p + ggplot2::facet_grid(. ~ pdim)
  }
  p
}


#' Plot negative log posterior likelihood of samples
#'
#' @param SGGP SGGP object
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' gs <- SGGPcreate(d=3, batchsize=100)
#' f <- function(x){x[1]^1.2+x[3]^.4*sin(2*pi*x[2]^2*3) + .1*exp(3*x[3])}
#' y <- apply(gs$design, 1, f)
#' gs <- SGGPfit(gs, Y=y)
#' SGGPplotsamplesneglogpost(gs)
SGGPplotsamplesneglogpost <- function(SGGP) {
  neglogpost_thetaMAP <- SGGP_internal_neglogpost(SGGP$thetaMAP, SGGP, SGGP$y)
  neglogpost_samples <- apply(SGGP$thetaPostSamples, 2, SGGP_internal_neglogpost, SGGP=SGGP, y=SGGP$y)
  ggplot2::ggplot() + 
    ggplot2::geom_histogram(ggplot2::aes_string(x='x'), data.frame(x=neglogpost_samples), bins=30) + 
    ggplot2::geom_vline(xintercept = neglogpost_thetaMAP, color="blue", size=2) +
    ggplot2::xlab("neglogpost") + ggplot2::ggtitle("neglogpost of theta samples (blue is MAP)")
}
