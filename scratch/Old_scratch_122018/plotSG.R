# I'm trying to create ways to visualize the sparse grid designs.

ggplot2::ggplot(reshape2::melt(data.frame(SG$designindex)), ggplot2::aes(x=value)) + ggplot2::geom_histogram(binwidth = 1) + ggplot2::facet_grid(variable ~ .)
ggplot2::ggplot(reshape2::melt(data.frame(SG$designindex)), ggplot2::aes(x=value)) + ggplot2::geom_histogram(binwidth = 1) + ggplot2::facet_grid(variable ~ .) + ggplot2::scale_y_log10()

heatmatrix <- matrix(NaN, SG$d, SG$d)
skinny <- NULL
for (i in 1:SG$d) {
  heatmatrix[i,i] <- max(SG$designindex[,i])
  skinny <- rbind(skinny, c(i, i, max(SG$designindex[,i])))
}
for (i in 1:(SG$d-1)) {
  for (j in (i+1):SG$d) {
    heatmatrix[i,j] <- heatmatrix[j,i] <- max(apply(SG$designindex[,c(i,j)], 1, min))
    skinny <- rbind(skinny,
                    c(i, j, max(apply(SG$designindex[,c(i,j)], 1, min))),
                    c(j, i, max(apply(SG$designindex[,c(i,j)], 1, min)))
                    )
  }
}
heatmatrix
image(heatmatrix)
heatmap(heatmatrix, Rowv=NA, Colv=NA, symm=T)
lattice::levelplot(heatmatrix)
# Try to make heatmap look better
heatmap(heatmatrix, Rowv=NA, symm=T, 
        # RowSideColors = terrain.colors(max(diag(heatmatrix)))[diag(heatmatrix)], 
        RowSideColors = heat.colors(max(diag(heatmatrix)))[(max(diag(heatmatrix))):1][diag(heatmatrix)],
        col=heat.colors(32)[32:1])
text(x=.5,y=.5, labels=c('.5'))
text(x=1,y=.5, labels=c('1,.5'))
text(x=0,y=.5, labels=c('0,.5'))
text((0:4)/5, y=0, labels=1:5)
# Can't get text to show up in spots easily

# Try ggplot2 https://stackoverflow.com/questions/14290364/heatmap-with-values-ggplot2
skdf <- data.frame(skinny)
names(skdf) <- c('Var1', 'Var2', 'value')
ggplot(skdf, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "red") 
# Looks good