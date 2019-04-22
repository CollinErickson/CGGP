# createbox
x1 <- c(0,.5,1)
x2 <- c(.5)
pts <- expand.grid(x1, x2)
box1 <- 1
mat <- matrix(c(1,1,1,2,1,3,2,1,1,1,2,2,3,1,4,1), ncol=2, byrow=T)
require(ggplot2); require(CGGP)
CGGPplotblocks(mat) + coord_fixed()

blocks_points <- function(SG, mat, scale=.9, includebelow=FALSE, b_plot=TRUE, add_beneath_points) {
  use_xb <- (SG$xb-.5)*scale + .5
  # mat <- matrix(c(1,1,1,2,1,3,2,1,1,1), ncol=2, byrow=T)
  # a refers to points in 0 to 1
  # b refers to points location in the block shown. So the point in 0 to 1
  #   plus the block depth in that direction
  aptsall <- NULL
  bptsall <- NULL
  blocksall <- NULL
  for (i in 1:nrow(mat)) {
    a1 <- mat[i,1]
    a2 <- mat[i,2]
    cssizes <- c(0, cumsum(SG$sizes))
    if (includebelow) {
      
      a1inds <- 1:sum(SG$sizes[1:a1]) # + cssizes[a1]
      a2inds <- 1:sum(SG$sizes[1:a2]) # + cssizes[a2]
    } else {
      
      a1inds <- 1:SG$sizes[a1] + cssizes[a1]
      a2inds <- 1:SG$sizes[a2] + cssizes[a2]
    }
    a1pts <- use_xb[a1inds]
    b1pts <- a1pts + a1-1
    a2pts <- use_xb[a2inds]
    b2pts <- a2pts + a2-1
    apts <- expand.grid(X1=a1pts, X2=a2pts)
    bpts <- expand.grid(X1=b1pts, X2=b2pts)
    aptsall <- rbind(aptsall, apts)
    bptsall <- rbind(bptsall, bpts)
    blocksall <- rbind(blocksall, data.frame(xmin=a1-1, xmax=a1,
                                             ymin=a2-1, ymax=a2))
  }
  
  if (!missing(add_beneath_points)) {
    p <- add_beneath_points
  } else {
    p <- ggplot()
  }
  if (b_plot) {
    p <- p + xlim(0,max(mat)) + ylim(0,max(mat)) +
      geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), data=blocksall, color='black', fill='yellow') + 
      geom_point(data=bptsall, mapping=aes(X1,X2)) + xlab(expression(i[1])) + ylab(expression(i[2]))
  } else {
    p <- p + coord_fixed(xlim=c(0,1), ylim=c(0,1)) +#xlim(0,1) + ylim(0,1) +
      geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), data=data.frame(xmin=0, xmax=1, ymin=0, ymax=1), color='black', fill='palegreen1', alpha=0) + 
      geom_point(data=(aptsall-.5)/scale+.5, mapping=aes(X1,X2), size=2) + xlab(expression(x[1])) + ylab(expression(x[2]))
  }
  p
}

eg4 <- expand.grid(1:4,1:4)
SG <- CGGPcreate(d=2, batchsize=100)
# Full factorial
blocks_points(SG, eg4) + coord_fixed()
blocks_points(SG, eg4, b_plot=F) + coord_fixed()
# Sparse grid
blocks_points(SG, eg4[rowSums(eg4)<6,]) + coord_fixed()
blocks_points(SG, eg4[rowSums(eg4)<6,], b_plot=F) + coord_fixed()
# Adaptive grid 1, the cross
blocks_points(SG, eg4[apply(eg4, 1, function(x)any(x==1)),]) + coord_fixed()
blocks_points(SG, eg4[apply(eg4, 1, function(x)any(x==1)),], b_plot=F) + coord_fixed()
# Adaptive grid 2, focus on D1
blocks_points(SG, eg4[c(1,2,3,4,5,6,9, 7),]) + coord_fixed()
blocks_points(SG, eg4[c(1,2,3,4,5,6,9, 7),], b_plot=F) + coord_fixed()
# Full factorial with all below
blocks_points(SG, eg4, includebelow = T) + coord_fixed()
# Sparse grid with all below
blocks_points(SG, eg4[rowSums(eg4)<6,], includebelow = T) + coord_fixed()


# SAVE IMAGES
SAVEPLOT <- T
maybe_save <- function(filepath, p,folderpath="./scratch/thesis/", device='eps', width=4, height=4) {
  if (SAVEPLOT) {
    ggsave(paste0(folderpath, "/", filepath, ".", device), p, device=device, width=width, height=height, units="in")
  } else {p}
}
# Full factorial
maybe_save("CGGP2DexFF", blocks_points(SG, eg4) + coord_fixed())
maybe_save("CGGP2DexFF_a", blocks_points(SG, eg4, b_plot=F) + coord_fixed())
# Sparse grid
maybe_save("CGGP2DexSG", blocks_points(SG, eg4[rowSums(eg4)<6,]) + coord_fixed())
maybe_save("CGGP2DexSG_a", blocks_points(SG, eg4[rowSums(eg4)<6,], b_plot=F) + coord_fixed())
# Adaptive grid 1, the cross
maybe_save("CGGP2DexAGcross", blocks_points(SG, eg4[apply(eg4, 1, function(x)any(x==1)),]) + coord_fixed())
maybe_save("CGGP2DexAGcross_a", blocks_points(SG, eg4[apply(eg4, 1, function(x)any(x==1)),], b_plot=F) + coord_fixed())
# Adaptive grid 2, focus on D1
maybe_save("CGGP2DexAGd1", blocks_points(SG, eg4[c(1,2,3,4,5,6,9, 7),]) + coord_fixed())
maybe_save("CGGP2DexAGd1_a", blocks_points(SG, eg4[c(1,2,3,4,5,6,9, 7),], b_plot=F) + coord_fixed())


# Demonstration of selecting points
d <- 2
f <- function(x) {(1+.1*x[2]^1.1)*(1-.5*x[1]) + (.1*x[2]+.5)*cos(2*pi*x[1]^1.3)}
f <- function(x) {(1+.1*x[2])*(1-.5*x[1]) + (.1*x[2]+.5)*cos(2*pi*x[1])}
ContourFunctions::cf(f)
c1 <- CGGPcreate(d, 5+12)
c1 <- CGGPfit(c1, apply(c1$design, 1, f))
blocks_points(c1, c1$uo[1:c1$uoCOUNT,])
c1 <- CGGPfit(c1, apply(c1$design, 1, f))
c1 <- CGGPappend(c1, 32)
blocks_points(c1, c1$uo[1:c1$uoCOUNT,])
c1 <- CGGPfit(c1, apply(c1$design, 1, f))
c1 <- CGGPappend(c1, 64)
blocks_points(c1, c1$uo[1:c1$uoCOUNT,])
blocks_points(c1, c1$uo[1:c1$uoCOUNT,], b_plot=F)
epsil <- .0; t1 <- expand.grid(seq(0-epsil,1+epsil,l=101), seq(0,1,l=101))
y1 <- apply(t1, 1, f)
t2 <- cbind(t1, y=y1)
ggplot() + geom_contour(aes(x=Var1, y=Var2, z=y), t2)
ggplot() + geom_contour(aes(x=Var1, y=Var2, z=y), t2)+ geom_raster(aes(x=Var1, y=Var2, fill = y), t2) +  geom_contour(colour = "white")
## Wrong order
blocks_points(c1, c1$uo[1:c1$uoCOUNT,], b_plot=F) + geom_contour(aes(x=Var1, y=Var2, z=y), t2)+ geom_raster(aes(x=Var1, y=Var2, fill = y), t2, alpha=.6) +  geom_contour(colour = "white")
ggplot() + geom_contour(aes(x=Var1, y=Var2, z=y), t2)+ geom_raster(aes(x=Var1, y=Var2, fill = y), t2) +  geom_contour(colour = "white") + blocks_points(c1, c1$uo[1:c1$uoCOUNT,], b_plot=F)
# Here's a good one
blocks_points(c1, c1$uo[1:c1$uoCOUNT,], b_plot=F, add_beneath_points = ggplot() + geom_raster(aes(x=Var1, y=Var2, fill = y), t2)+
                scale_fill_gradientn(colours=c("#639fff","#FFFFFFFF","#ff5959")) +  geom_contour(aes(x=Var1, y=Var2, z=y), t2, colour = "white"))
# Same but with no lines
blocks_points(c1, c1$uo[1:c1$uoCOUNT,], b_plot=F, add_beneath_points = ggplot() + geom_raster(aes(x=Var1, y=Var2, fill = y), t2)+
                scale_fill_gradientn(colours=c("#639fff","#FFFFFFFF","#ff5959")))

# Save these plots
SAVEPLOT <- T
c1 <- CGGPcreate(d, 5+12)
c1 <- CGGPfit(c1, apply(c1$design, 1, f))
maybe_save("CGGPsequentialdemoa", blocks_points(c1, c1$uo[1:c1$uoCOUNT,]))
c1 <- CGGPfit(c1, apply(c1$design, 1, f))
c1 <- CGGPappend(c1, 32)
maybe_save("CGGPsequentialdemob", blocks_points(c1, c1$uo[1:c1$uoCOUNT,]))
c1 <- CGGPfit(c1, apply(c1$design, 1, f))
c1 <- CGGPappend(c1, 64)
maybe_save("CGGPsequentialdemoc", blocks_points(c1, c1$uo[1:c1$uoCOUNT,]))
maybe_save("CGGPsequentialdemod", blocks_points(c1, c1$uo[1:c1$uoCOUNT,], b_plot=F))
# Save d with contour plot underneath
maybe_save("CGGPsequentialdemod_cont", blocks_points(c1, c1$uo[1:c1$uoCOUNT,], b_plot=F, add_beneath_points = ggplot() + geom_raster(aes(x=Var1, y=Var2, fill = y), t2)+
                                                       scale_fill_gradientn(colours=c("#639fff","#FFFFFFFF","#ff5959"))))
