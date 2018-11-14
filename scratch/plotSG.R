# I'm trying to create ways to visualize the sparse grid designs.

ggplot2::ggplot(reshape2::melt(data.frame(SG$designindex)), ggplot2::aes(x=value)) + ggplot2::geom_histogram(binwidth = 1) + ggplot2::facet_grid(variable ~ .)
ggplot2::ggplot(reshape2::melt(data.frame(SG$designindex)), ggplot2::aes(x=value)) + ggplot2::geom_histogram(binwidth = 1) + ggplot2::facet_grid(variable ~ .) + ggplot2::scale_y_log10()

heatmatrix <- matrix(NA, SG$d, SG$d)
for (i in 1:SG$d) {
  heatmatrix[i,i] <- max(SG$designindex[,i])
}
for (i in 1:(SG$d-1)) {
  for (j in (i+1):SG$d) {
    heatmatrix[i,j] <- heatmatrix[j,i] <- max(apply(SG$designindex[,c(i,j)], 1, min))
  }
}
heatmatrix
image(heatmatrix)
heatmap(heatmatrix, Rowv=NA, Colv=NA, symm=T)
lattice::levelplot(heatmatrix)