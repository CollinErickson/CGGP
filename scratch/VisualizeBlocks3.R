require(ggplot2); require(CGGP); require(dplyr)

plot_index <- function(i1=1, i2=1) {
  cg <- CGGPcreate(2, 5)
  cg$uo[1:cg$uoCOUNT,]
  # i1 <- 3
  # i2 <- 2
  Xi1 <- cg$xb[1:cumsum(cg$sizes)[i1]]
  Xi2 <- cg$xb[1:cumsum(cg$sizes)[i2]]
  tdf <- data.frame(x=.5, y=.5)
  tdf <- expand.grid(Xi1, Xi2)
  colnames(tdf) <- c("x1", "x2")
  ggplot() + geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), data.frame(xmin=0, xmax=1, ymin=0, ymax=1), fill='yellow', color="black") + 
    geom_point(aes(x=x1, y=x2), tdf, size=2) + geom_point(aes(x=x1, y=x2), data.frame(x1=Xi1, x2=-.1), size=2.5)+geom_point(aes(x=x1, y=x2), data.frame(x1=-.1, x2=Xi2), size=2.5) +
    xlab(expression(x[1])) + ylab(expression(x[2])) + ggtitle(paste0("i = (",i1,", ",i2,")")) + coord_fixed()
}
plot_index()
plot_index(2,3)

# maybe_save from VisualizeBlocks.R
maybe_save("IndexExample11", plot_index(1,1), height=3, width=3)
maybe_save("IndexExample12", plot_index(1,2), height=3, width=3)
maybe_save("IndexExample21", plot_index(2,1), height=3, width=3)



plot_indexset <- function(i1=1, i2=1) {
  message("plot_indexset is basically the same think as blocks_points with b_plot=F")
  cg <- CGGPcreate(2, 5)
  cg$uo[1:cg$uoCOUNT,]
  ddf <- NULL
  if (length(i1) != length(i2)) {stop("i1 != length i2")}
  for (ii in 1:length(i1)) {
    t1 <- i1[ii]
    t2 <- i2[ii]
    Xi1 <- cg$xb[1:cumsum(cg$sizes)[t1]]
    Xi2 <- cg$xb[1:cumsum(cg$sizes)[t2]]
    tdf <- data.frame(x=.5, y=.5)
    tdf <- expand.grid(Xi1, Xi2)
    colnames(tdf) <- c("x1", "x2")
    ddf <- rbind(ddf, tdf)
  }
  gt <- "I = {"
  for (ii in 1:length(i1)) {if (ii>1) {gt <- paste0(gt, ", ")};gt <- paste0(gt, "(", i1[ii], ", ", i2[ii], ")")}
  gt <- paste0(gt, "}")
  ggplot() + geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), data.frame(xmin=0, xmax=1, ymin=0, ymax=1), fill='yellow', color="black", alpha=0) + 
    geom_point(aes(x=x1, y=x2), ddf, size=2) + #geom_point(aes(x=x1, y=x2), data.frame(x1=Xi1, x2=-.1), size=2.5)+geom_point(aes(x=x1, y=x2), data.frame(x1=-.1, x2=Xi2), size=2.5) +
    xlab(expression(x[1])) + ylab(expression(x[2])) + ggtitle(gt) + coord_fixed()
}
plot_indexset()
plot_indexset(c(1,2), c(1,1))
plot_indexset(c(1,2,1), c(1,1,2))

maybe_save("IndexSetExample_11_21_12", plot_indexset(c(1,2,1), c(1,1,2)), height=3, width=3)
