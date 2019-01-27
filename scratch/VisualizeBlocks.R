# createbox
x1 <- c(0,.5,1)
x2 <- c(.5)
pts <- expand.grid(x1, x2)
box1 <- 1
mat <- matrix(c(1,1,1,2,1,3,2,1,1,1,2,2,3,1,4,1), ncol=2, byrow=T)
SGGPblockplot(mat) + coord_fixed()

blocks_points <- function(SG, mat, scale=.9, includebelow=FALSE) {
  use_xb <- (SG$xb-.5)*scale + .5
  # mat <- matrix(c(1,1,1,2,1,3,2,1,1,1), ncol=2, byrow=T)
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
    bptsall <- rbind(bptsall, bpts)
    blocksall <- rbind(blocksall, data.frame(xmin=a1-1, xmax=a1,
                                             ymin=a2-1, ymax=a2))
  }
  
  ggplot() + xlim(0,max(mat)) + ylim(0,max(mat)) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), data=blocksall, color='black', fill='pink') + 
    geom_point(data=bptsall, mapping=aes(X1,X2))
}

eg4 <- expand.grid(1:4,1:4)
blocks_points(SG, eg4)
blocks_points(SG, eg4[rowSums(eg4)<6,])
blocks_points(SG, eg4, includebelow = T)
blocks_points(SG, eg4[rowSums(eg4)<6,], includebelow = T)
